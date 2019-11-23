!-----------------------------------------------------------------------
!     train.x - train (fit) atomic energy neural network potential
!-----------------------------------------------------------------------
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2019 Nongnuch Artrith and Alexander Urban
!+
!+ This Source Code Form is subject to the terms of the Mozilla Public
!+ License, v. 2.0. If a copy of the MPL was not distributed with this
!+ file, You can obtain one at http://mozilla.org/MPL/2.0/.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!+ Mozilla Public License, v. 2.0, for more details.
!+ ---------------------------------------------------------------------
!+ If you make use of AENET for your publication, please cite:
!+ [1] N. Artrith and A. Urban, Comput. Mater. Sci. 114 (2016) 135-150.
!+ [2] J. Behler and M. Parrinello, Phys. Rev. Lett. 98 (2007) 146401.
!+
!+ If you used the Chebyshev descriptor, please cite:
!+ [3] N. Artrith, A. Urban, and G. Ceder, PRB 96 (2017) 014112.
!-----------------------------------------------------------------------
! 2011-10-19 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

program train

  use aeio,        only: aeio_readline,          &
                         aeio_header,            &
                         aeio_timestamp,         &
                         aeio_print_copyright,   &
                         PATHLEN, TYPELEN

  use feedforward, only: Network,                &
                         new_Network,            &
                         del_Network,            &
                         save_Network,           &
                         load_Network,           &
                         ff_change_activation,   &
                         ff_random_init_weights, &
                         ff_update_weights,      &
                         ff_get_nweights,        &
                         ff_print_info,          &
                         ff_eval,                &
                         ff_deriv,               &
                         ff_wderiv

  use geometry,    only: geo_itype_of_name

  use input,       only: InputData,              &
                         read_InpTrain,          &
                         del_InputData,          &
                         inp_read_networks

  use io,          only: io_adjustl,             &
                         io_center,              &
                         io_lower,               &
                         io_readval,             &
                         io_readnext,            &
                         io_unit

  use optimize,    only: opt_init_training,      &
                         opt_final,              &
                         opt_before_batch,       &
                         opt_after_sample,       &
                         opt_after_batch,        &
                         opt_print_info,         &
                         opt_schedule_epoch,     &
                         opt_schedule,           &
                         opt_samplesize_local,   &
                         opt_batchsize_local,    &
                         opt_nbatch

  use parallel,    only: pp_init,                &
                         pp_final,               &
                         pp_print_info,          &
                         pp_bcast,               &
                         pp_recv,                &
                         pp_send,                &
                         pp_sum,                 &
                         pp_sum2d,               &
                         pp_bcast_Network,       &
                         pp_bcast_InputData,     &
                         pp_bcast_TrnSet_info,   &
                         pp_sum_weights,         &
                         ppMaster,               &
                         ppRank,                 &
                         ppSize

  use random,      only: random_init,            &
                         random_reinit,          &
                         random_final,           &
                         random_save_state,      &
                         random_load_state

  use sfsetup,     only: Setup,                  &
                         load_Setup,             &
                         save_Setup,             &
                         skip_Setup,             &
                         del_Setup

  use timing,      only: tng_init,               &
                         tng_final,              &
                         tng_timing,             &
                         tng_timing2,            &
                         tng_timing3,            &
                         tng_dump

  use trainset,    only: TrnSet,                  &
                         new_TrnSet,              &
                         open_TrnSet,             &
                         rewind_TrnSet,           &
                         close_TrnSet,            &
                         save_TrnSet_info,        &
                         ts_count_atoms,          &
                         ts_load_Setups,          &
                         ts_print_info,           &
                         ts_read_atom_info,       &
                         ts_read_sf_info,         &
                         ts_read_sf_values,       &
                         ts_read_structure_info,  &
                         ts_write_atom_info,      &
                         ts_write_sf_info,        &
                         ts_write_structure_info, &
                         ts_skip_atoms

  implicit none

  !--------------------------------------------------------------------!
  ! RNG_STATE_FILE    filename to save state of the random number      !
  !                   generator, needed for restarting the testing set !
  ! STOP_FILE         name of the file that causes training to stop    !
  !--------------------------------------------------------------------!

  character(len=*), parameter :: RNG_STATE_FILE = 'train.rngstate'
  character(len=*), parameter :: STOP_FILE = 'STOP'

  !--------------------------------------------------------------------!
  ! A '*' in front of the variable name means that it is a broadcasted !
  ! variable and has the same value on each process.  A '+' means that !
  ! an array is allocated on all parallel processes, but does not      !
  ! necessarily have the same contents.                                !
  !                                                                    !
  !----------------------------- general ------------------------------!
  ! inFile          name of the input file                             !
  !*inp             input data object                                  !
  !                                                                    !
  !--------------------------- training set ---------------------------!
  ! restarted       .true. if the testing set is restarted             !
  ! ts              training set; instance of TrnSet                   !
  ! E_av, E_min, E_max  average, minimum and maximum cohesive energy   !
  !                 per atom in the training set structures            !
  !*Escale          the cohesive energy will be normalized to the      !
  !                 interval [-1:1] before the training                !
  !                 E_new = Escale*E                                   !
  !*nTrain          number of structures actually used for training    !
  !*nTest           size of the testing set = ts%nStrucs - nTrain      !
  !*nAtomsTrain     total number of atoms in the training set          !
  !*nAtomsTest      total number of atoms in the testing set           !
  !*isTest(i)       .true., if the i-th structure belongs to the       !
  !                 testing set                                        !
  !+trnEnergies(i)  reference energy of the i-th structure             !
  !+trnErrors(i)    current error of the i-th structure                !
  !                                                                    !
  !------------------------- training method --------------------------!
  !*iepoch          current epoch/training iteration                   !
  !*conv            .true., if training has converged                  !
  !*batchsize       size of the sliding window batch, if any           !
  !*nextbatch       .true., if the window shall be updated             !
  ! ibatch          process local batch counter                        !
  ! batchiter       iterations for the same batch so far               !
  !                                                                    !
  !------------------- structural fingerprint basis -------------------!
  ! stp(itype)      basis function set-up for atom type itype          !
  !*nsf_max         max. number of basis functions over all set-ups    !
  !                                                                    !
  !----------------------------- networks -----------------------------!
  !*net(itype)      neural network for central atom type itype         !
  !*nw_max          max. number of network weights over all networks   !
  !*nw_tot          total number (sum) of network weights              !
  ! iw, nw          generic weight counters                            !
  !+Dw(k,l)         = sum_i Dw_i(k) ; for atom type l                  !
  !+ann_values(i)   value of the i-th ANN node                         !
  !+ann_derivs(i)   derivative of the i-th ANN node value              !
  !+ann_jacobian(i) derivative of the ANN output wrt. the i-th weight  !
  !                                                                    !
  !----------------------------- parallel -----------------------------!
  !+ts_trn, ts_tst  process local train and test data sets             !
  !*stopnow         if .true., stop immediately                        !
  !--------------------------------------------------------------------!

  character(len=100)                               :: inFile
  type(InputData)                                  :: inp

  logical                                          :: restarted
  type(TrnSet)                                     :: ts
  integer                                          :: nTrain
  integer                                          :: nTest
  integer                                          :: nAtomsTrain
  integer                                          :: nAtomsTest
  logical,           dimension(:),     allocatable :: isTest

  integer                                          :: iepoch
  logical                                          :: conv
  integer                                          :: batchsize
  logical                                          :: do_nextbatch
  integer                                          :: ibatch
  double precision, dimension(:),      allocatable :: trnEnergies
  double precision, dimension(:),      allocatable :: trnErrors

  logical                                          :: do_deriv

  type(Setup),       dimension(:),     allocatable :: stp
  integer                                          :: nsf_max

  type(Network),     dimension(:),     allocatable :: net
  integer                                          :: nnodes_max
  integer                                          :: nw_max
  integer                                          :: nw_tot
  double precision,  dimension(:,:),   allocatable :: Dw
  double precision,  dimension(:),     allocatable :: ann_values
  double precision,  dimension(:),     allocatable :: ann_derivs
  double precision,  dimension(:),     allocatable :: ann_jacobian

  type(TrnSet)                                     :: ts_trn, ts_tst
  integer                                          :: itrn
  logical                                          :: stopnow

  double precision                                 :: dE
  double precision                                 :: MAE_trn, MAE_tst
  double precision                                 :: SSE_trn, SSE_tst
  double precision                                 :: RMSE_trn, RMSE_tst

  integer                                          :: u_dbg, u_tng

  integer                                          :: irec, isample

  ! timing registers
  integer, parameter :: R_TRN = 1, R_TST = 2

  external :: DSYR2K ! BLAS

  !-------------------------- initialization --------------------------!

  call initialize(inFile)

  stopnow = .false.
  if (ppMaster) then
     inp = read_InpTrain(inFile)
     if (.not. inp%init) then
        stopnow = .true.
     else
        ! load training set
        ts = open_TrnSet(inp%trn_file, maxenergy=inp%trn_maxenergy)
        ! load basis function set ups
        allocate(stp(ts%nTypes))
        call ts_load_Setups(ts,stp)
        ! synchronize file types and read NN file names
        inp%nTypes = ts%nTypes
        allocate(inp%typeName(inp%nTypes))
        inp%typeName(:) = ts%typeName(:)
        call inp_read_networks(inp, file=trim(inp%file), readarch=.true.)
        allocate(net(ts%nTypes))
        ! set up NNs
        call init_networks(inp, stp, net)
        ! max basis functions and weights
        nsf_max = max_num_sf()
        call get_num_nodes(ts, nnodes_max)
        call get_num_weights(ts, nw_max, nw_tot)
     end if
  end if
  call pp_bcast(stopnow)
  if (stopnow) then
     call finalize()
     stop
  end if
  call pp_bcast_InputData(inp)
  call pp_bcast_TrnSet_info(ts)
  call broadcast_networks()
  call pp_bcast(nsf_max)
  call pp_bcast(nnodes_max)
  call pp_bcast(nw_max)
  call pp_bcast(nw_tot)

  ! allocate memory for ANN values, derivatives, and Jacobian
  allocate(ann_values(nnodes_max), &
           ann_derivs(nnodes_max), &
           ann_jacobian(nw_max))

  if (inp%do_timing .and. ppMaster) then
     u_tng = io_unit()
     call tng_init(unit=u_tng, file='train.time', registers=2)
     write(*,*) 'Timing info will be written to: train.time'
     write(*,*)
  end if
  if (inp%do_debug) then
     u_dbg = io_unit()
     open(u_dbg, file='train.debug'//trim(io_adjustl(ppRank)), &
          status='replace', action='write')
  end if

  if (ppMaster) call save_all_networks(verbose=.true.)

  !------------- set-up of the training and testing sets --------------!

  allocate(isTest(ts%nStrucs))
  if (ppMaster) then
     ! devide structures into test and train set:
     call decide_testing_set(ts%nStrucs, inp%trn_testset, nTrain, nTest, &
                             isTest, restarted)
     ! count total number of atoms in train and test structures:
     call rewind_TrnSet(ts)
     call ts_count_atoms(ts, isTest, nAtomsTrain, nAtomsTest)
     call ts_print_info(ts)
     call print_training_info()
  end if
  call pp_bcast(isTest, ts%nStrucs)
  call pp_bcast(nAtomsTrain)
  call pp_bcast(nTrain)
  call pp_bcast(nAtomsTest)
  call pp_bcast(nTest)

  if (inp%do_timing .and. ppMaster) &
       call tng_timing('Training set initialized.')

  ! distribute the training and testing structures:
  call distribute_trnfile(isTest, nsf_max, ts, ts_trn, ts_tst)

  ! allocate arrays for process local energies and errors
  allocate(trnEnergies(ts_trn%nStrucs), trnErrors(ts_trn%nStrucs))
  ! the reference energies don't change. Store once and for all.
  call store_training_set_energies(ts_trn, trnEnergies)

  if (inp%do_timing .and. ppMaster) &
       call tng_timing('Structures distributed over processes')

  !------------------ training method initialization ------------------!

  call opt_init_training(inp%trn_method, inp%trn_param, inp%trn_sampling, &
                         nw_tot, nw_max, nTrain, ts_trn%nStrucs, &
                         ts_trn%nTypes)

  allocate(Dw(nw_max,ts%nTypes))

  if (inp%do_timing .and. ppMaster) then
     call tng_timing('Optimization method (' // trim(inp%trn_method) &
                     // ') initialized.')
  end if

  !--------------- initial status (before any training) ---------------!

  if (ppMaster) &
     call print_training_header(inp%trn_steps, inp%trn_methodName)

  call eval_entire_trainset(ts_trn, nTrain, nsf_max, nw_max, net, &
                            MAE_trn, SSE_trn, errors=trnErrors)
  call eval_entire_trainset(ts_tst, nTest, nsf_max, nw_max, net, &
                            MAE_tst, SSE_tst)

  if (ppMaster) call print_energies(0, &
       MAE_trn/ts_trn%scale, sqrt(2.0d0*SSE_trn/dble(nTrain))/ts_trn%scale, &
       MAE_tst/ts_trn%scale, sqrt(2.0d0*SSE_tst/dble(nTest))/ts_tst%scale)

  if (inp%do_timing .and. ppMaster) &
       call tng_timing('Initial energies evaluated (before training)')

  !----------------------------- training -----------------------------!

  conv = .false.
  epochs : do iepoch = 1, inp%trn_steps

     if (conv) then
        if (ppMaster) then
           write(*,*) 'The optimization has converged. Training stopped.'
           write(*,*)
        end if
        exit epochs
     end if

     call opt_schedule_epoch(ts_trn%nStrucs, trnEnergies, trnErrors)

     if (inp%do_timing .and. ppMaster) &
          call tng_timing('Starting epoch ' // io_adjustl(iepoch))

     ibatch       = 1
     do_nextbatch = .true.
     do_deriv     = .true.

     batches : do while (ibatch <= opt_nbatch)

        call opt_before_batch()

        !--------------!
        ! training set !
        !--------------!

        irec = (ibatch - 1)*opt_batchsize_local
        training : do itrn = 1, opt_batchsize_local

           ! go to next record in training set file
           ! after the last record, start from beginning
           irec = mod(irec, opt_samplesize_local) + 1
           isample = opt_schedule(irec)

           if (isample /= (ts_trn%iStruc + 1)) then
              call rewind_TrnSet(ts_trn, rec=isample)
           end if

           if (do_deriv) then
              ! get error function and gradient
              call eval_next_structure(ts_trn, nsf_max, nw_max, net, dE, Dw=Dw)
           else
              ! only compute error
              call eval_next_structure(ts_trn, nsf_max, nw_max, net, dE)
           end if

           call opt_after_sample(net, ts, dE, Dw)

        end do training

        call opt_after_batch(net, ts, do_deriv, do_nextbatch, conv)

        if (inp%do_timing .and. ppMaster) then
           call tng_timing2('Done with training iteration.')
           call tng_timing3(register=R_TRN)
        end if

        if (do_nextbatch) ibatch = ibatch + 1

     end do batches

     !-----------------------------!
     ! current training set errors !
     !-----------------------------!

     ! FIXME: evaluation here not needed when using batch training
     call eval_entire_trainset(ts_trn, nTrain, nsf_max, nw_max, net, &
                               MAE_trn, SSE_trn, errors=trnErrors)

     if (inp%do_timing .and. ppMaster) then
        call tng_timing2('Done with training set energies.')
        call tng_timing3(register=R_TRN)
     end if

     !-------------------------!
     ! current test set errors !
     !-------------------------!

     call eval_entire_trainset(ts_tst, nTest, nsf_max, nw_max, net, &
                               MAE_tst, SSE_tst)

     if (inp%do_timing .and. ppMaster) then
        call tng_timing2('Done with testing set energies.')
        call tng_timing3(register=R_TST)
     end if

     if (ppMaster) then
        RMSE_trn = sqrt(2.0d0*SSE_trn/dble(nTrain))
        RMSE_tst = sqrt(2.0d0*SSE_tst/dble(nTest))
        call print_energies(iepoch, &
             MAE_trn/ts_trn%scale, RMSE_trn/ts_trn%scale, &
             MAE_tst/ts_trn%scale, RMSE_tst/ts_tst%scale  )
        call save_all_networks(iter=iepoch)
     end if

     ! synchronize networks (to avoid numerical problems)
     if (ppSize>1) call broadcast_networks()

     if (inp%do_timing .and. ppMaster) &
          call tng_timing2('Done with epoch '//trim(io_adjustl(iepoch)))

     if (check_stopfile()) exit epochs

  end do epochs

  if (ppMaster) then
     write(*,*)
     write(*,*) 'Training finished.'
     write(*,*)
  end if

  if (inp%do_timing .and. ppMaster) then
     call tng_timing('Training finished.')
     call tng_dump(R_TRN, 'Total time spent for training set.')
     call tng_dump(R_TST, 'Total time spent for testing set.')
  end if

  !---------------- save final energies, if requested -----------------!

  if (inp%do_save_energies) then
     call save_all_energies(&
          ts_trn, ts_tst, nTrain, nTest, nsf_max, nw_max, net)
     if (inp%do_timing .and. ppMaster) &
          call tng_timing('Stored energies of all structures.')
  end if

  !--------------------------- finalization ---------------------------!

  call finalize()


contains !=============================================================!


  subroutine initialize(inFile)

    implicit none

    character(len=*), intent(out) :: inFile

    logical :: fexists
    integer :: nargs
    logical :: stopnow

    call pp_init()

    stopnow = .false.
    if (ppMaster) then
       call aeio_header("Training process started.", char='=')
       call aeio_header(aeio_timestamp(), char=' ')
       write(*,*)

       call aeio_print_copyright('2015-2018', 'Nongnuch Artrith and Alexander Urban')

       nargs = command_argument_count()
       if (nargs < 1) then
          write(0,*) "Error: No input file provided."
          call print_usage()
          stopnow = .true.
       end if

       call get_command_argument(1, value=inFile)
       inquire(file=trim(inFile), exist=fexists)
       if (.not. fexists) then
          write(0,*) "Error: File not found: ", trim(inFile)
          call print_usage()
          stopnow = .true.
       end if

       call random_init()
    end if

    call pp_bcast(stopnow)
    if (stopnow) then
       call finalize()
       stop
    end if

    call pp_print_info()

  end subroutine initialize

  !--------------------------------------------------------------------!

  subroutine init_networks(inp, stp, net)

    implicit none

    type(InputData),                intent(in)    :: inp
    type(Setup),      dimension(:), intent(in)    :: stp
    type(Network),    dimension(:), intent(inout) :: net

    integer           :: i, ipos, itype
    integer           :: il, nl, nlmax
    logical           :: fexists
    character(len=20) :: nodes

    integer,           dimension(:), allocatable :: arch
    character(len=10), dimension(:), allocatable :: ftype

    nlmax = 10
    allocate(arch(nlmax), ftype(nlmax))

    call aeio_header("Networks")
    write(*,*)

    do itype = 1, inp%nTypes
       inquire(file=trim(inp%netFile(itype)), exist=fexists)

       if (fexists) then
          ! restart the network:
          write(*,*) 'Restarting the ', trim(inp%typeName(itype)), &
                     ' network from file : ', trim(inp%netFile(itype))
          net(itype) = load_Network(trim(inp%netFile(itype)))
       else
          ! set up a new network:
          write(*,*) 'Creating a new ', trim(inp%typeName(itype)), ' network'
          ipos = 1
          call io_readnext(inp%netArch(itype), ipos, nl)
          if (nl+2 > nlmax) then
             ! reallocate arrays, if more memory is needed
             nlmax = nl + 2
             if (allocated(arch)) deallocate(arch, ftype)
             allocate(arch(nlmax), ftype(2:nlmax))
          end if
          arch(1)     = stp(itype)%nsf
          arch(nl+2)  = 1
          ftype(nl+2) = 'linear'
          do il = 1, nl
             call io_readnext(inp%netArch(itype), ipos, nodes)
             i = scan(nodes, ':')
             read(nodes(1:i-1), *)          arch(il+1)
             read(nodes(i+1:len(nodes)), *) ftype(il+1)
          end do
          net(itype) = new_Network(arch(1:nl+2))
          do il = 2, nl+2
             call ff_change_activation(net(itype), il, ftype(il))
          end do
          call ff_random_init_weights(net(itype))
       end if ! restart

       ! write out network information:
       write(*,*)
       call ff_print_info(net(itype))

    end do ! itype

  end subroutine init_networks

  !--------------------------------------------------------------------!

  subroutine finalize()

    implicit none

    integer :: itype

    if (ppMaster) then
       ! write current networks to files:
       call save_all_networks(verbose=.true.)
    end if

    ! deallocate memory needed for optimization:
    call opt_final()

    if (allocated(trnEnergies)) deallocate(trnEnergies, trnErrors)

    if (allocated(ann_values)) deallocate(ann_values, ann_derivs, ann_jacobian)

    if (allocated(stp)) then
       do itype = 1, ts%nTypes
          call del_Setup(stp(itype))
       end do
       deallocate(stp)
    end if

    if (allocated(net)) then
       do itype = 1, ts%nTypes
          call del_Network(net(itype))
       end do
       deallocate(net)
    end if

    ! close training sets
    if (ts%init)     call close_TrnSet(ts)
    if (ts_trn%init) call close_TrnSet(ts_trn, status='delete')
    if (ts_tst%init) call close_TrnSet(ts_tst, status='delete')

    if (allocated(isTest))  deallocate(isTest)
    if (allocated(Dw))      deallocate(Dw)

    if (ppMaster) then
       call aeio_header(aeio_timestamp(), char=' ')
       call aeio_header("Neural Network training done.", char='=')
       call random_final()
    end if

    if (inp%do_debug) close(u_dbg)

    if (inp%do_timing .and. ppMaster) call tng_final()

    call pp_final()

  end subroutine finalize

  !--------------------------------------------------------------------!

  function check_stopfile() result(fexists)

    implicit none

    logical :: fexists

    if (ppMaster) then
       inquire(file=STOP_FILE, exist=fexists)
       if (fexists) write(*,*) "File `STOP' detected.  Training halted."
    end if
    call pp_bcast(fexists)

  end function check_stopfile

  !--------------------------------------------------------------------!

  subroutine print_usage()

    implicit none

    write(*,*)
    write(*,*) "train.x -- Train an atomic energy NN."
    write(*,'(1x,70("-"))')
    write(*,*) 'Usage: train.x <input-file>'
    write(*,*)
    write(*,*) 'See the documentation or the source code for a description of the '
    write(*,*) 'input file format.'
    write(*,*)

  end subroutine print_usage

  !--------------------------------------------------------------------!
  !             target energies of all training structures             !
  !--------------------------------------------------------------------!

  subroutine store_training_set_energies(ts, energies)

    implicit none

    type(TrnSet),                   intent(inout) :: ts
    double precision, dimension(:), intent(out)   :: energies

    character(len=1024) :: filename
    integer             :: nAtoms
    integer             :: nTypes
    double precision    :: E_coh
    integer             :: i

    call rewind_TrnSet(ts)
    do i = 1, ts%nStrucs
       call ts_read_structure_info(ts, filename, nAtoms, nTypes, E_coh)
       energies(i) = E_coh
       call ts_skip_atoms(ts, nAtoms)
    end do
    call rewind_TrnSet(ts)

  end subroutine store_training_set_energies

  !--------------------------------------------------------------------!
  !                     actual network evaluation                      !
  !--------------------------------------------------------------------!

  !--------------------- evaluate whole data set ----------------------!

  subroutine eval_entire_trainset(ts, N, nsf_max, nw_max, net, &
                                  MAE, SSE, errors, save_energies)

    implicit none

    !------------------------------------------------------------------!
    ! ts          training set instance to be read from                !
    ! N           total number of structures for normalization         !
    !             (not necessarily equal to ts%nStrucs)                !
    ! nsf_max     max. number of basis functions                       !
    ! nw_max      max. number of network weights                       !
    ! net(i,j)    NN for atom types i and j                            !
    ! MAE         mean absolute error                                  !
    ! SSE         sum of squared errors                                !
    !------------------------------------------------------------------!

    type(TrnSet),                             intent(inout) :: ts
    integer,                                  intent(in)    :: N
    integer,                                  intent(in)    :: nsf_max
    integer,                                  intent(in)    :: nw_max
    type(Network),    dimension(:),           intent(inout) :: net
    double precision,                         intent(out)   :: MAE
    double precision,                         intent(out)   :: SSE
    double precision, dimension(:), optional, intent(out)   :: errors
    integer,                        optional, intent(in)    :: save_energies


    double precision :: dE, error
    integer :: i

    MAE = 0.0d0
    SSE = 0.0d0

    call rewind_TrnSet(ts)
    do i = 1, ts%nStrucs
       if (present(save_energies)) then
          call eval_next_structure(ts, nsf_max, nw_max, net, dE, &
                                   save_energies=save_energies)
       else
          call eval_next_structure(ts, nsf_max, nw_max, net, dE)
       end if
       error = 0.5d0*dE*dE
       if (present(errors)) errors(i) = error
       MAE = MAE + abs(dE)/dble(N)
       SSE = SSE + error
    end do

    if (ppSize > 1) then
       call pp_sum(MAE)
       call pp_sum(SSE)
    end if

  end subroutine eval_entire_trainset

  !-------------------- evaluate single structure ---------------------!

  subroutine eval_next_structure(ts, nsf_max, nw_max, net, dE, dF, Dw, &
                                 save_energies)
    implicit none

    !------------------------------------------------------------------!
    ! ts          training set instance to be read rom                 !
    ! nsf_max     max. number of basis functions                       !
    ! net(i,j)    NN for atom types i and j                            !
    ! dE          energy difference per atom: dE=(E_coh-E_nn)/nAtoms   !
    ! dF(i,j)     difference in the i-th force component of the j-th   !
    !             atom: dF(i,j)=F_nn(i,j)-F(i,j)  [optional]           !
    !             Note: force evaluation not yet implemented here !!   !
    ! Dw          Jacobian with weight derivatives [optional]          !
    !------------------------------------------------------------------!

    type(TrnSet),                                 intent(inout) :: ts
    integer,                                      intent(in)    :: nsf_max
    integer,                                      intent(in)    :: nw_max
    type(Network),    dimension(:),               intent(inout) :: net
    double precision,                             intent(out)   :: dE
    double precision, dimension(:,:),   optional, intent(out)   :: dF
    double precision, dimension(:,:),   optional, intent(out)   :: Dw
    integer,                            optional, intent(in)    :: save_energies

    !------------------------------------------------------------------!
    ! iatom, nAtoms  counter and total number for atoms in structure   !
    ! itype, nTypes  counter and total number for species in structure !
    ! itype          the atomic species of the current network         !
    ! E_coh          cohesive (target) energy                          !
    ! forCart(1:3)   Cartesian (target) force of current atom          !
    ! cooCart(1:3)   Cartesian coordinates of current atom             !
    ! nsf, nw        number of basis functions and weights             !
    ! sfval(i)       value of the i-th basis function                  !
    ! E_i            predicted atomic energy of the current atom       !
    ! E_nn           total NN predicted energy = sum_i E_i             !
    ! F_i(i)         derivative of E_i wrt. i-th basis function        !
    ! Dw_i(i)        derivative of E_i wrt. i-th weight                !
    !------------------------------------------------------------------!

    integer                                :: iatom, nAtoms
    integer                                :: nTypes
    integer                                :: itype
    character(len=1024)                    :: filename
    double precision                       :: E_coh
    double precision, dimension(3)         :: forCart, cooCart
    integer                                :: nsf, nw
    double precision, dimension(nsf_max)   :: sfval
    double precision                       :: E_nn, E_i
    double precision, dimension(nsf_max)   :: F_i
    double precision, dimension(nw_max)    :: Dw_i
    integer,          dimension(ts%nTypes) :: natoms_type

    if (present(save_energies)) natoms_type(1:ts%nTypes) = 0

    E_nn = 0.0d0
    if (present(Dw)) Dw(:,:) = 0.0d0
    call ts_read_structure_info(ts, filename, nAtoms, nTypes, E_coh)
    do iatom = 1, nAtoms
       call ts_read_atom_info(ts, itype, cooCart, forCart)
       call ts_read_sf_info(ts, nsf)
       call ts_read_sf_values(ts, nsf, sfval(1:nsf))
       if (present(Dw)) then
          nw = ff_get_nweights(net(itype))
          call eval_net(net(itype), nsf, nw, sfval(1:nsf), &
                        E_i, F_i(1:nsf), Dw_i(1:nw))
          Dw(1:nw, itype) = Dw(1:nw, itype) - Dw_i(1:nw)/dble(nAtoms)
       else
          call eval_net2(net(itype), nsf, sfval(1:nsf), E_i, F_i(1:nsf))
       end if
       E_nn = E_nn + E_i
       if (present(save_energies)) natoms_type(itype) = natoms_type(itype) + 1
    end do ! iatom central atom

    dE = (E_coh - E_nn)/dble(nAtoms)

    if (present(save_energies)) then
       ! save energies to file for comparison
       E_coh = E_coh/ts%scale + dble(nAtoms)*ts%shift
       E_nn  = E_nn/ts%scale + dble(nAtoms)*ts%shift
       write(save_energies,'(1x,2(ES14.6,1x),I5,1x,4(ES14.6,1x))', &
            advance='no') E_coh, E_nn, nAtoms, E_coh/dble(nAtoms), &
            E_nn/dble(nAtoms), dE/ts%scale, dE
       ! write out number of atoms for each species
       do itype = 1, ts%nTypes
          write(save_energies, '(1x,I4)', advance='no') natoms_type(itype)
       end do
       write(save_energies, '(1x,A)') trim(filename)
    end if
    if (present(dF)) dF(:,:) = 0.0d0

  end subroutine eval_next_structure

  !--------------------------------------------------------------------!

  subroutine eval_net(net, nsf, nw, sfval, E, F, Dw)

    implicit none

    type(Network),                    intent(inout) :: net
    integer,                          intent(in)    :: nsf
    integer,                          intent(in)    :: nw
    double precision, dimension(nsf), intent(in)    :: sfval
    double precision,                 intent(out)   :: E
    double precision, dimension(nsf), intent(out)   :: F
    double precision, dimension(nw),  intent(out)   :: Dw

    double precision, dimension(1) :: Ebuff

    ann_values(:) = 0.0d0
    ann_derivs(:) = 0.0d0
    ann_jacobian(:) = 0.0d0
    call ff_eval(net, nsf, sfval(1:nsf), 1, ann_values, ann_derivs, Ebuff)
    call ff_deriv(net, nsf, 1, ann_derivs, ann_jacobian, F(1:nsf))
    call ff_wderiv(net, nw, 1, ann_values, ann_derivs, ann_jacobian, Dw)

    E = Ebuff(1)

  end subroutine eval_net

  !--------------------------------------------------------------------!

  subroutine eval_net2(net, nsf, sfval, E, F)

    implicit none

    type(Network),                    intent(inout) :: net
    integer,                          intent(in)    :: nsf
    double precision, dimension(nsf), intent(in)    :: sfval
    double precision,                 intent(out)   :: E
    double precision, dimension(nsf), intent(out)   :: F

    double precision, dimension(1) :: Ebuff

    ann_values(:) = 0.0d0
    ann_derivs(:) = 0.0d0
    call ff_eval(net, nsf, sfval(1:nsf), 1, ann_values, ann_derivs, Ebuff)
    call ff_deriv(net, nsf, 1, ann_derivs, ann_jacobian, F(1:nsf))

    E = Ebuff(1)

  end subroutine eval_net2

  !--------------------------------------------------------------------!
  !                           weight updates                           !
  !--------------------------------------------------------------------!

  subroutine update_weights(nTypes, Dw, net)

    implicit none

    integer,                                    intent(in)    :: nTypes
    double precision, dimension(:,:), optional, intent(in)    :: Dw
    type(Network),    dimension(:),             intent(inout) :: net

    integer :: itype
    integer :: nw

    do itype = 1, nTypes
       nw = ff_get_nweights(net(itype))
       call ff_update_weights(net(itype), nw, Dw(1:nw,itype))
    end do

  end subroutine update_weights

  !--------------------------------------------------------------------!
  !      all processes (in a parallel run) must know the networks      !
  !--------------------------------------------------------------------!

  subroutine broadcast_networks()

    implicit none

    integer :: itype

    call pp_bcast(ts%nTypes)
    if (.not. allocated(net)) allocate(net(ts%nTypes))
    do itype = 1, ts%nTypes
       call pp_bcast_Network(net(itype))
    end do

  end subroutine broadcast_networks

  !--------------------------------------------------------------------!
  !            divide training / testing set over processes            !
  !--------------------------------------------------------------------!

  subroutine distribute_trnfile(isTest, nsf_max, ts, ts_trn, ts_tst)

    implicit none

    logical, dimension(:), intent(in)    :: isTest
    integer,               intent(in)    :: nsf_max
    type(TrnSet),          intent(inout) :: ts
    type(TrnSet), target,  intent(out)   :: ts_trn
    type(TrnSet), target,  intent(out)   :: ts_tst

    integer                                 :: nTrn, nTst
    character(len=1024)                     :: procTrnFile, procTstFile
    integer                                 :: iproc
    integer                                 :: ifile
    character(len=1024)                     :: filename
    integer                                 :: iatom, natoms
    integer                                 :: itype, ntypes
    integer                                 :: nsf
    double precision, dimension(3)          :: coo, for
    double precision, dimension(nsf_max)    :: sfval
    !$$ double precision, dimension(3,nsf_max)  :: sfderiv
    double precision                        :: energy
    integer,          dimension(0:ppSize-1) :: nAtProc
    type(TrnSet), pointer                   :: ts_p
    double precision                        :: scale, shift

    ! determine how many structures each process will receive
    nTrn = 0
    nTst = 0
    do ifile = 1, ts%nStrucs
       iproc = mod(ifile, ppSize)
       if (iproc == ppRank) then
          if (isTest(ifile)) then
             nTst = nTst + 1
          else
             nTrn = nTrn + 1
          end if
       end if
    end do

    if (ppMaster) then
       call rewind_TrnSet(ts)
       scale = ts%scale
       shift = ts%shift
    end if
    call pp_bcast(scale)
    call pp_bcast(shift)

    ! file names for process local data set files
    procTrnFile = 'TRAIN.'//trim(io_adjustl(ppRank))
    procTstFile = 'TEST.'//trim(io_adjustl(ppRank))

    ts_trn = new_TrnSet(ts%nTypes, ts%typeName, ts%E_atom, nTrn, &
                        procTrnFile, scale=scale, shift=shift)
    ts_tst = new_TrnSet(ts%nTypes, ts%typeName, ts%E_atom, nTst, &
                        procTstFile, scale=scale, shift=shift)

    nAtProc(:) = 0
    structures : do ifile = 1, ts%nStrucs

       ! target process:
       iproc = mod(ifile, ppSize)

       ! target file:
       if (isTest(ifile)) then
          ts_p => ts_tst
       else
          ts_p => ts_trn
       end if

       !-----------------------------------------------!
       ! read record and send it to the target process !
       !-----------------------------------------------!

       if (ppMaster) then
          call ts_read_structure_info(ts, filename, natoms, ntypes, energy)
          if (iproc == ppRank) then
             call ts_write_structure_info(ts_p, filename, natoms, ntypes, energy)
          else
             call pp_send(filename, dest=iproc)
             call pp_send(natoms,   dest=iproc)
             call pp_send(ntypes,   dest=iproc)
             call pp_send(energy,   dest=iproc)
          end if
          nAtProc(iproc) = nAtProc(iproc) + natoms
          do iatom = 1, natoms
             call ts_read_atom_info(ts, itype, coo, for)
             if (iproc == ppRank) then
                call ts_write_atom_info(ts_p, itype, coo, for)
             else
                call pp_send(itype, dest=iproc)
                call pp_send(coo, 3, dest=iproc)
                call pp_send(for, 3, dest=iproc)
             end if
             call ts_read_sf_info(ts, nsf)
             call ts_read_sf_values(ts, nsf, sfval(1:nsf))
             if (iproc == ppRank) then
                call ts_write_sf_info(ts_p, nsf, sfval(1:nsf))
                !$$ call ts_write_sf_info(ts_p, itype2, nsf, sfval(1:nsf), &
                !$$                       sfderiv(1:3,1:nsf))
             else
                call pp_send(nsf, dest=iproc)
                call pp_send(sfval(1:nsf), nsf, dest=iproc)
                !$$ do isf = 1, nsf
                !$$    call pp_send(sfderiv(1:3,isf), 3, dest=iproc)
                !$$ end do
             end if
          end do ! iatom
       end if ! ppMaster

       !-----------------------------------------------!
       ! receive record and write it to the right file !
       !-----------------------------------------------!

       if ((iproc == ppRank) .and. (.not. ppMaster)) then
          call pp_recv(filename)
          call pp_recv(natoms)
          call pp_recv(ntypes)
          call pp_recv(energy)
          call ts_write_structure_info(ts_p, filename, natoms, ntypes, energy)
          do iatom = 1, natoms
             call pp_recv(itype)
             call pp_recv(coo, 3)
             call pp_recv(for, 3)
             call ts_write_atom_info(ts_p, itype, coo, for)
             call pp_recv(nsf)
             call pp_recv(sfval(1:nsf), nsf)
             !$$ do isf = 1, nsf
             !$$    call pp_recv(sfderiv(1:3,isf), 3)
             !$$ end do
             call ts_write_sf_info(ts_p, nsf, sfval(1:nsf))
             !$$ call ts_write_sf_info(ts_p, itype2, nsf, sfval(1:nsf), &
             !$$                       sfderiv(1:3,1:nsf))
          end do ! iatom
       end if

    end do structures

    call close_TrnSet(ts_trn)
    call close_TrnSet(ts_tst)
    ts_trn = open_TrnSet(procTrnFile)
    ts_tst = open_TrnSet(procTstFile)

    if (ppMaster) then
       ! write out info about static work load balance:
       if (ppSize > 1) then
          call aeio_header("Static work load balance")
          write(*,*)
          write(*,*) 'Number of atoms on the different processes'
          write(*,*)
          do iproc = 0, ppSize-1
             write(*,'(1x,I4,2x,":",2x,I8)') iproc, nAtProc(iproc)
          end do
          write(*,*)
       end if
    end if

  end subroutine distribute_trnfile

  !--------------------------------------------------------------------!
  !               print info about the training settings               !
  !--------------------------------------------------------------------!

  subroutine print_training_info()

    implicit none

    integer :: ifile, itest

    call aeio_header('Training details')
    write(*,*)

    write(*,*) 'Training method         : ', trim(inp%trn_methodName)
    select case(trim(inp%trn_method))
    case('ekf')
       write(*,*) '  Forgetting factor     : ', inp%trn_param(1)
       write(*,*) '  Ini. state covariance : ', inp%trn_param(2)
       write(*,*) '  Measuring noise       : ', inp%trn_param(3)
       write(*,*) '  Process noise         : ', inp%trn_param(4)
       write(*,*) '  Adaptive EKF threshold: ', inp%trn_param(5)
    case('lm')
       write(*,*) '  Batch size            : ', int(inp%trn_param(1))
       write(*,*) '  Batch iterations      : ', int(inp%trn_param(3))
       write(*,*) '  Initial learning rate : ', inp%trn_param(2)
       write(*,*) '  Target error          : ', inp%trn_param(4)
    case('online_sd')
       write(*,*) '  Learning rate         : ', inp%trn_param(2)
       write(*,*) '  Momentum rate         : ', inp%trn_param(1)
    end select
    write(*,*)
    write(*,*) 'Number of iterations    : ', trim(io_adjustl(inp%trn_steps))
    write(*,*)
    write(*,*) 'Training structures     : ', trim(io_adjustl(nTrain))
    write(*,*) 'Testing  structures     : ', trim(io_adjustl(nTest))
    write(*,*)

    if (restarted) then
       write(*,'(1x,"Testing set (restarted from previous run): ")')
    else
       write(*,'(1x,"Testing set : ")')
    end if
    itest = 0
    do ifile = 1, ts%nStrucs
       if (isTest(ifile)) then
          if (mod(itest,8) == 0) then
             write(*,*)
          end if
          itest = itest + 1
          write(*, '(I8,1x)', advance='no') ifile
       end if
    end do
    write(*,*)
    write(*,*)

  end subroutine print_training_info

  !--------------------------------------------------------------------!
  !                print header info for training steps                !
  !--------------------------------------------------------------------!

  subroutine print_training_header(nIters, method)

    implicit none

    integer,          intent(in) :: nIters
    character(len=*), intent(in) :: method

    call aeio_header("Training process")
    write(*,*)
    write(*,*) 'Weight optimization for ' // trim(io_adjustl(nIters)) &
         // ' epochs using the ' // trim(adjustl(method)) &
         // ' method.'
    write(*,*)

    call opt_print_info()

    write(*,*)
    write(*,'(8x,A30,2x,A30)') &
         '|------------TRAIN-----------|', &
         '|------------TEST------------|'
    write(*,'(1x,A5,2x,A14,2x,A14,2x,A14,2x,A14)') &
         'epoch', 'MAE', '<RMSE>', 'MAE', '<RMSE>'

  end subroutine print_training_header

  !------------------ energies at current iteration -------------------!

  subroutine print_energies(istep, MAE_trn, RMSE_trn, MAE_tst, RMSE_tst)

    implicit none

    integer,          intent(in) :: istep
    double precision, intent(in) :: MAE_trn
    double precision, intent(in) :: RMSE_trn
    double precision, intent(in) :: MAE_tst
    double precision, intent(in) :: RMSE_tst

    write(*,'(1x,I5,2x,ES14.6,2x,ES14.6,2x,ES14.6,2x,ES14.6," <")') &
          istep, MAE_trn, RMSE_trn, MAE_tst, RMSE_tst

  end subroutine print_energies

  !--------------------------------------------------------------------!
  !  save all networks to files (incl. structural fingerprint set-up)  !
  !--------------------------------------------------------------------!

  subroutine save_all_networks(verbose, iter)

    implicit none

    logical, optional, intent(in) :: verbose
    integer, optional, intent(in) :: iter

    integer :: itype
    integer :: u_sav
    character(len=10) :: enum

    if (.not. allocated(net)) return
    if (.not. allocated(stp)) return

    if (present(verbose)) then
       call aeio_header("Storing networks")
       write(*,*)
    end if

    if (present(iter)) then
       write(enum,'("-",I0.5)') iter
    else
       enum = " "
    end if

    u_sav = io_unit()
    do itype = 1, ts%nTypes
       if (present(verbose)) then
          write(*,*) 'Saving the ', trim(ts%typeName(itype)), &
               ' network to file : ', trim(inp%netFile(itype)) // trim(enum)
       end if
       open(u_sav, file=trim(inp%netFile(itype))//trim(enum), status='replace', &
            action='write', form='unformatted')
       call save_Network(net(itype), unit=u_sav)
       call save_Setup(stp(itype), unit=u_sav)
       call save_TrnSet_info(ts, unit=u_sav)
       close(u_sav)
    end do
    if (present(verbose)) write(*,*)

  end subroutine save_all_networks

  !--------------------------------------------------------------------!
  !                    save all structural energies                    !
  !--------------------------------------------------------------------!

  subroutine save_all_energies(ts_trn, ts_tst, nTrain, nTest, nsf_max, &
                               nw_max, net)

    implicit none

    type(TrnSet),                intent(inout) :: ts_trn
    type(TrnSet),                intent(inout) :: ts_tst
    integer,                     intent(in)    :: nTrain
    integer,                     intent(in)    :: nTest
    integer,                     intent(in)    :: nsf_max
    integer,                     intent(in)    :: nw_max
    type(Network), dimension(:), intent(inout) :: net

    character(len=PATHLEN) :: fname
    character(len=64)      :: frmt, str
    character(len=1024)    :: header
    integer                :: u_sav, itype
    double precision       :: MAE_trn, SSE_trn, RMSE_trn
    double precision       :: MAE_tst, SSE_tst, RMSE_tst

    if (ppMaster) then
       call aeio_header("Storing final energies")
       write(*,*)
       write(*,*) 'Energies of training structures : energies.train.PROCESS'
       write(*,*) 'Energies of testing structures  : energies.test.PROCESS'
       write(*,*) '(Manually concatenate the files from different processes.)'
       write(*,*)
    end if

    u_sav = io_unit()

    ! file header with column description
    header = "  Ref(eV)        ANN(eV)      #atoms  Ref(eV/atom)" &
             // "   ANN(eV/atom) Ref-ANN(eV/atom)    Cost-Func"
    frmt = '(2x,"#",A' // trim(io_adjustl(TYPELEN)) // ')'
    write(str, frmt) ts_trn%typeName(1)
    header = trim(header) // " " // trim(str)
    do itype = 2, ts_trn%nTypes
       write(str, frmt) ts_trn%typeName(itype)
       header = trim(header) // trim(str)
    end do
    header = trim(header) // "    Path-of-input-file"

    ! training set
    fname = 'energies.train.' // trim(io_adjustl(ppRank))
    open(u_sav, file=trim(fname), status='replace', action='write')
    write(u_sav, '(A)') trim(header)
    call eval_entire_trainset(ts_trn, nTrain, nsf_max, nw_max, net, &
                              MAE_trn, SSE_trn, save_energies=u_sav)
    close(u_sav)

    ! testing set
    fname = 'energies.test.' // trim(io_adjustl(ppRank))
    open(u_sav, file=trim(fname), status='replace', action='write')
    write(u_sav, '(A)') trim(header)
    call eval_entire_trainset(ts_tst, nTest, nsf_max, nw_max, net, &
                              MAE_tst, SSE_tst, save_energies=u_sav)
    close(u_sav)

    if (ppMaster) then
       MAE_trn  = 1000.0d0*MAE_trn/ts_trn%scale
       RMSE_trn = 1000.0d0*sqrt(2.0d0*SSE_trn/dble(nTrain))/ts_trn%scale
       MAE_tst  = 1000.0d0*MAE_tst/ts_tst%scale
       RMSE_tst = 1000.0d0*sqrt(2.0d0*SSE_tst/dble(nTest))/ts_tst%scale

       write(*,'(1x,"Final MAE of training set  = ",F8.1," meV/atom")') MAE_trn
       write(*,'(1x,"Final MAE of testing set   = ",F8.1," meV/atom")') MAE_tst
       write(*,*)
       write(*,'(1x,"Final RMSE of training set = ",F8.1," meV/atom")') RMSE_trn
       write(*,'(1x,"Final RMSE of testing set  = ",F8.1," meV/atom")') RMSE_tst
       write(*,*)
    end if

  end subroutine save_all_energies

  !--------------------------------------------------------------------!
  !                 maximum number of basis functions                  !
  !--------------------------------------------------------------------!

  function max_num_sf() result(nsf)

    implicit none

    integer :: nsf
    integer :: itype

    nsf = 0
    do itype = 1, ts%nTypes
       nsf = max(nsf, stp(itype)%nsf)
    end do

  end function max_num_sf

  !--------------------------------------------------------------------!
  !                 maximum number of network nodes                    !
  !--------------------------------------------------------------------!

  subroutine get_num_nodes(ts, nnodes_max)

    implicit none

    type(TrnSet), intent(in)  :: ts
    integer,      intent(out) :: nnodes_max

    integer :: itype

    nnodes_max = 0
    do itype = 1, ts%nTypes
       nnodes_max = max(nnodes_max, net(itype)%nvalues)
    end do

  end subroutine get_num_nodes

  !--------------------------------------------------------------------!
  !                 maximum number of network weights                  !
  !--------------------------------------------------------------------!

  subroutine get_num_weights(ts, nw_max, nw_tot)

    implicit none

    type(TrnSet), intent(in)  :: ts
    integer,      intent(out) :: nw_max, nw_tot

    integer :: nw, itype

    nw_max = 0
    nw_tot = 0
    do itype = 1, ts%nTypes
       nw = ff_get_nweights(net(itype))
       nw_max = max(nw_max, nw)
       nw_tot = nw_tot + nw
    end do

  end subroutine get_num_weights

  !--------------------------------------------------------------------!
  !           decide, which structures enter the testing set           !
  !--------------------------------------------------------------------!

  subroutine decide_testing_set(nFiles, testpercent, nTrain, nTest, &
                                isTest, restarted)

    implicit none

    integer,                    intent(in)  :: nFiles
    double precision,           intent(in)  :: testpercent
    integer,                    intent(out) :: nTrain
    integer,                    intent(out) :: nTest
    logical, dimension(nFiles), intent(out) :: isTest
    logical,                    intent(out) :: restarted

    logical          :: fexists
    double precision :: r
    integer          :: ifile, itrain, itest

    inquire(file=RNG_STATE_FILE, exist=fexists)
    if (fexists) then
       restarted = .true.
       call random_load_state(file=RNG_STATE_FILE, unit=io_unit())
    else
       restarted = .false.
       call random_save_state(file=RNG_STATE_FILE, unit=io_unit())
    end if

    nTrain = ceiling((1.0d0 - testpercent/100.0d0)*dble(nFiles))
    nTest  = nFiles - nTrain

    itrain = 0
    itest  = 0
    do ifile = 1, nFiles
       if (itest == nTest) then
          itrain = itrain + 1
          isTest(ifile) = .false.
       else if (itrain == nTrain) then
          itest = itest + 1
          isTest(ifile) = .true.
       else
          call random_number(r)
          if (r <= testpercent/100.0d0) then
             itest = itest + 1
             isTest(ifile) = .true.
          else
             itrain = itrain + 1
             isTest(ifile) = .false.
          end if
       end if
    end do

    ! re-initialize random number generator, in case it is used for
    ! something else afterwards
    call random_reinit()

  end subroutine decide_testing_set

end program train
