!-----------------------------------------------------------------------
!     optimize.f90 - an interface module to optimization libraries
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
! 2011-11-18 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

module optimize

  use feedforward, only: Network,         &
                         ff_get_nweights, &
                         ff_update_weights

  use io,          only: io_adjustl,      &
                         io_lower,        &
                         io_unit

  use parallel,    only: ppSize,          &
                         ppRank,          &
                         ppMaster,        &
                         pp_sum,          &
                         pp_sum2d,        &
                         pp_bcast,        &
                         pp_bcast_Network

  use random,      only: random_reinit

  use sortlib,     only: argsort

  use trainset,    only: TrnSet

  implicit none
  save

  public  :: opt_init,             &
             opt_init_training,    &
             opt_final,            &
             opt_save_state,       &
             opt_load_state,       &
             opt_before_batch,     &
             opt_after_sample,     &
             opt_after_batch,      &
             opt_optimize_coords,  &
             opt_optimize,         &
             opt_print_info,       &
             opt_schedule_epoch

  private :: opt_defaults,          &
             setulb,                &
             opt_init_sd,           &
             opt_final_sd,          &
             opt_after_sample_sd,   &
             opt_after_batch_sd,    &
             opt_init_ekf,          &
             opt_final_ekf,         &
             opt_save_state_ekf,    &
             opt_load_state_ekf,    &
             opt_after_sample_ekf,  &
             ekf_weight_update,     &
             opt_init_lm,           &
             opt_final_lm,          &
             opt_before_batch_lm,   &
             opt_after_sample_lm,   &
             opt_after_batch_lm,    &
             lm_jacobian,           &
             lm_inner_loop,         &
             opt_init_bfgs,         &
             opt_final_bfgs,        &
             opt_save_state_bfgs,   &
             opt_load_state_bfgs,   &
             opt_after_batch_bfgs,  &
             opt_bfgs_weights,      &
             opt_bfgs_coords,       &
             opt_assert_moduleinit, &
             opt_update_weights

  !---------------------------- constants -----------------------------!
  ! OPT_STATE_FILE    name of the binary restart file                  !
  !--------------------------------------------------------------------!

  character(len=*), parameter :: OPT_STATE_FILE = 'train.restart'

  !----------------------------- general ------------------------------!
  ! method       string ID of the optimization method used             !
  ! memsize      estimate number of allocated words for optimizer      !
  ! ftol, gtol   convergence criteria for function and radient         !
  ! batchsize    size of training batch (# of samples per iteration)   !
  ! batchsize_local size of the process local batch                    !
  ! nbatch       number of batches per epoch                           !
  ! batchiter    current bach iteration                                !
  ! nbatchiters  total number iterations per batch                     !
  ! nw_tot       total number of weights in weight optimization        !
  ! nw_max       max. number of weights per atomic NN                  !
  ! ntypes       number of different atomic species                    !
  ! SSE          Sum of Squared Errors                                 !
  ! sampling_type 'random', 'sequential', or 'weighted'                !
  ! samplesize   Number of samples for each epoch                      !
  ! samplesize_local Process local number of samples for each epoch    !
  ! schedule     Indices of the samples for one epoch                  !
  ! idx          temporary sort index                                  !
  !--------------------------------------------------------------------!

  character(len=50), public :: opt_method
  integer,           public :: opt_memsize
  double precision,  public :: opt_ftol
  double precision,  public :: opt_gtol
  integer,           public :: opt_batchsize
  integer,           public :: opt_batchsize_local
  integer,           public :: opt_nbatch
  integer,           public :: opt_batchiter
  integer,           public :: opt_nbatchiters
  integer,           public :: opt_nw_tot
  integer,           public :: opt_nw_max
  integer,           public :: opt_ntypes
  double precision,  public :: opt_SSE
  character(len=10), public :: opt_sampling_type
  integer,           public :: opt_samplesize
  integer,           public :: opt_samplesize_local

  integer, dimension(:), allocatable, public  :: opt_schedule

  integer, dimension(:), allocatable, private :: opt_idx

  ! optional work space
  double precision, dimension(:,:),   allocatable, private :: A, B, C

  !--------------- steepest descent method (online_sd) ----------------!
  ! momentum    weight of the previous iteration                       !
  ! learnrate   learning rate = amount of weight gradient to use       !
  ! Dw_prev     derivatives (Jacobian) of previous iteration           !
  ! Dw_sum      cumulative sum of derivatives (batch Jacobian)         !
  !--------------------------------------------------------------------!

  double precision,                              private :: sd_momentum
  double precision,                              private :: sd_learnrate
  double precision, dimension(:,:), allocatable, private :: sd_Dw_prev
  double precision, dimension(:,:), allocatable, private :: sd_Dw_sum

  !---------------- Adaptive Moment Estimation (Adam) -----------------!
  ! learnrate   initial learning rate                                  !
  ! b1          decay rate for gradients                               !
  ! b2          decay rate for squared gradients                       !
  ! eps         small parameter to prevent division by zero            !
  ! samplesize  number of randomly selected samples per epoch          !
  ! batchsize   size of the mini-batch for Jacobian calculation        !
  ! m           exponentially decaying average of the gradient         !
  ! v           exponentially decaying average of the squared gradient !
  ! m2 and v2   same, but renormalized not to tend towards zero        !
  ! t           step counter                                           !
  ! Dw_sum      cumulative sum of derivatives (batch Jacobian)         !
  ! Dw_up       weight update                                          !
  !--------------------------------------------------------------------!

  double precision,                              private :: adam_learnrate
  double precision,                              private :: adam_b1
  double precision,                              private :: adam_b2
  double precision,                              private :: adam_eps
  integer,                                       private :: adam_samplesize
  integer,                                       private :: adam_batchsize
  double precision, dimension(:,:), allocatable, private :: adam_m
  double precision, dimension(:,:), allocatable, private :: adam_v
  double precision, dimension(:,:), allocatable, private :: adam_m2
  double precision, dimension(:,:), allocatable, private :: adam_v2
  integer,                                       private :: adam_t
  double precision, dimension(:,:), allocatable, private :: adam_Dw_sum
  double precision, dimension(:,:), allocatable, private :: adam_Dw_up

  !------------------- parameters for the L-BFGS-B --------------------!
  ! Dw_sum      sum of derivatives (Jacobian of batch)                 !
  !--------------------------------------------------------------------!
  ! x(i)        i-th parameters to be optimized                        !
  ! g(i)        function gradient gradient dy/dx(i)                    !
  ! n           total number of parameters                             !
  ! m           LM-BFGS memory parameter; select 3 < m < 20            !
  ! iprint      output control option; 0 = no output                   !
  ! factr       factr*EPS is the convergence criterion;                !
  ! pgtol       convergence criterion for the gradient                 !
  ! task        on first entry = 'START'; values: 'FG', 'NEW_X', ...   !
  ! l(i)        lower bound of parameter i                             !
  ! u(i)        upper bound of parameter i                             !
  ! nbd(i)      bounds for param i -- 0=no; 1=lower; 2=both; 3=upper   !
  ! wa          double precision work array                            !
  ! iwa         integer work array                                     !
  ! csave       character working array                                !
  ! lsave       logical working array                                  !
  ! isave       integer working array                                  !
  ! dsave       double precision working array                         !
  !--------------------------------------------------------------------!

  double precision, dimension(:,:), allocatable, private :: bfgs_Dw_sum

  integer,                                       private :: bfgs_n
  integer,                                       private :: bfgs_m
  integer,                                       private :: bfgs_iprint
  double precision,                              private :: bfgs_factr
  double precision,                              private :: bfgs_pgtol
  character(len=60),                             private :: bfgs_task
  double precision, dimension(:),   allocatable, private :: bfgs_x
  double precision, dimension(:),   allocatable, private :: bfgs_g
  double precision, dimension(:),   allocatable, private :: bfgs_l
  double precision, dimension(:),   allocatable, private :: bfgs_u
  integer,          dimension(:),   allocatable, private :: bfgs_nbd
  double precision, dimension(:),   allocatable, private :: bfgs_wa
  integer,          dimension(:),   allocatable, private :: bfgs_iwa
  character(len=60),                             private :: bfgs_csave
  logical,          dimension(4),                private :: bfgs_lsave
  integer,          dimension(44),               private :: bfgs_isave
  double precision, dimension(29),               private :: bfgs_dsave

  !------------------ Levenberg-Marquardt algorithm -------------------!
  ! learnrate     initial learning rate (only used in first step)      !
  ! adjust        factor to adjust the dynamic learning rate           !
  ! cycle         current inner cycle of the LM algorithm              !
  ! stat          current LM status (evaluate just energy or Jacobian) !
  ! SSE_prev      Sum of Squared Errors of previous LM iteration       !
  ! iJw           index for Jw                                         !
  ! W_save(i)     i-th saved weight from previous iteration            !
  ! Jw(:,i)       Jacobian for i-th batch                              !
  ! dEv(i)        cumulative error of i-th batch                       !
  ! H(i,j)        approximate Hessian matrix                           !
  ! grad(i)       gradient wrt. i-th weight                            !
  ! Wup(i)        update for i-th weight                               !
  !--------------------------------------------------------------------!

  double precision,                              private :: lm_learnrate
  double precision,                              private :: lm_adjust
  integer,                                       private :: lm_cycle
  character(len=6),                              private :: lm_stat
  double precision,                              private :: lm_SSE_prev
  integer,                                       private :: lm_iJw
  double precision, dimension(:),   allocatable, private :: lm_W_save
  double precision, dimension(:,:), allocatable, private :: lm_Jw
  double precision, dimension(:),   allocatable, private :: lm_dEv
  double precision, dimension(:,:), allocatable, private :: lm_H
  double precision, dimension(:),   allocatable, private :: lm_grad
  double precision, dimension(:),   allocatable, private :: lm_Wup

  !---------------------- Extended Kalman Filter ----------------------!
  ! lambda      initial value of the time dependent forgetting factor  !
  ! lambda0     constant decay rate of lambda                          !
  ! P_initial   P (see below) will be set to 1/this times the identity !
  ! mnoise      amplitude of the measuring noise                       !
  ! pnoise      amplitude of the processing noise                      !
  ! adaptive    RMSE threshold for adaptive filtering                  !
  ! nwg         number of weight groups                                !
  ! wg(i)       size of weight group i                                 !
  ! wg_max      max size of weight group                               !
  ! J           Jacobian                                               !
  ! P           state/error covariance matrix                          !
  ! KT          transposed Kalman gain matrix                          !
  ! E           error vector                                           !
  ! Wup         weight update                                          !
  !--------------------------------------------------------------------!

  double precision,                                private :: ekf_lambda
  double precision,                                private :: ekf_lambda0
  double precision,                                private :: ekf_P_initial
  double precision,                                private :: ekf_mnoise
  double precision,                                private :: ekf_pnoise
  double precision,                                private :: ekf_adaptive
  integer,                                         private :: ekf_nwg
  integer,          dimension(:),     allocatable, private :: ekf_wg
  integer,                                         private :: ekf_wg_max
  double precision, dimension(:,:),   allocatable, private :: ekf_J
  double precision, dimension(:,:,:), allocatable, private :: ekf_P
  double precision, dimension(:,:),   allocatable, private :: ekf_KT
  double precision, dimension(:),     allocatable, private :: ekf_E
  double precision, dimension(:),     allocatable, private :: ekf_Wup

  !--------------------------------------------------------------------!
  !                             interfaces                             !
  !--------------------------------------------------------------------!

  interface ! from the L-BFGS-B library
     subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, wa, &
                       iwa, task, iprint, csave, lsave, isave, dsave)
       implicit none
       integer,                         intent(in)    :: n
       integer,                         intent(in)    :: m
       double precision, dimension(n),  intent(inout) :: x
       double precision, dimension(n),  intent(in)    :: l
       double precision, dimension(n),  intent(in)    :: u
       integer,          dimension(n),  intent(in)    :: nbd
       double precision,                intent(inout) :: f
       double precision, dimension(n),  intent(inout) :: g
       double precision,                intent(in)    :: factr
       double precision,                intent(in)    :: pgtol
       double precision, dimension(*),  intent(inout) :: wa
       integer,          dimension(*),  intent(inout) :: iwa
       character(len=60),               intent(inout) :: task
       integer,                         intent(in)    :: iprint
       character(len=60),               intent(inout) :: csave
       logical,          dimension(4),  intent(inout) :: lsave
       integer,          dimension(44), intent(inout) :: isave
       double precision, dimension(29), intent(inout) :: dsave
     end subroutine setulb
  end interface

  interface opt_optimize
     module procedure opt_optimize_coords
  end interface opt_optimize

  !--------------------------------------------------------------------!

  logical, private :: isInit = .false.

contains

  subroutine opt_defaults()

    implicit none

    opt_method = ''
    opt_memsize = 0
    opt_ftol = 0.0d0
    opt_gtol = 1.0d-6
    opt_batchsize = 0
    opt_batchsize_local = 0
    opt_nw_tot = 0
    opt_nw_max = 0
    opt_sampling_type = 'sequential'
    opt_samplesize = 0

  end  subroutine opt_defaults

  !--------------------------------------------------------------------!

  subroutine opt_init(method, nparam, ftol, gtol)

    implicit none

    character(len=*),                         intent(in) :: method
    integer,                                  intent(in) :: nparam
    double precision,               optional, intent(in) :: ftol
    double precision,               optional, intent(in) :: gtol

    if (isInit) then
       write(0,*) "Warning: repeated initialization of module `optimize'."
       return
    end if

    isInit = .true.

    call opt_defaults()
    opt_method = trim(adjustl(io_lower(method)))
    if (present(ftol)) opt_ftol = ftol
    if (present(gtol)) opt_gtol = gtol

    select case(trim(opt_method))
    case('bfgs')
       call opt_init_bfgs(nparam, opt_ftol, opt_gtol, opt_memsize)
    case default
       write(0,*) "Error: unknown optimization method: " // trim(opt_method)
    end select

  end subroutine opt_init

  !--------------------------------------------------------------------!

  subroutine opt_init_training(method, methodparam, sampling_type, &
                               nw_tot, nw_max, ntrain, ntrain_local, ntypes, &
                               ftol, gtol)

    implicit none

    character(len=*),               intent(in)  :: method
    double precision, dimension(:), intent(in)  :: methodparam
    character(len=*),               intent(in)  :: sampling_type
    integer,                        intent(in)  :: nw_tot
    integer,                        intent(in)  :: nw_max
    integer,                        intent(in)  :: ntrain
    integer,                        intent(in)  :: ntrain_local
    integer,                        intent(in)  :: nTypes
    double precision, optional,     intent(in)  :: ftol
    double precision, optional,     intent(in)  :: gtol

    integer :: localbatch, batchsize
    integer :: i

    if (isInit) then
       write(0,*) "Warning: repeated initialization of module `optimize'."
       return
    end if

    isInit = .true.

    call opt_defaults()
    opt_method = trim(adjustl(io_lower(method)))
    if (present(ftol)) opt_ftol = ftol
    if (present(gtol)) opt_gtol = gtol
    opt_nw_tot = nw_tot
    opt_nw_max = nw_max
    opt_ntypes = ntypes

    opt_sampling_type = trim(sampling_type)

    localbatch = ntrain_local
    batchsize = ntrain
    opt_nbatchiters = 1
    opt_batchiter = 0

    opt_samplesize = ntrain
    opt_samplesize_local = ntrain_local

    select case(trim(opt_method))
    case('bfgs')
       if (ppMaster) then
          call opt_init_bfgs(nw_tot, opt_ftol, opt_gtol, opt_memsize)
       end if
       allocate(bfgs_Dw_sum(nw_max,ntypes))
       opt_memsize = opt_memsize + nw_max*nTypes*2
       if (trim(opt_sampling_type) /= "sequential") then
          write(0,*) "Warning: Sampling type set to 'sequential' for BFGS."
          opt_sampling_type = 'sequential'
       end if
    case('ekf')
       call opt_init_ekf(nw_tot, methodparam, opt_memsize)
    case('lm')
       call opt_init_lm(nw_tot, ntrain, methodparam, batchsize, &
                        localbatch, opt_nbatchiters, opt_memsize)
    case('online_sd')
       call opt_init_sd(nw_max, ntypes, methodparam, opt_memsize)
    case('adam')
       call opt_init_adam(nw_max, ntypes, ntrain_local, methodparam, &
                          opt_memsize)
       opt_samplesize_local = adam_samplesize
       opt_samplesize = adam_samplesize
       batchsize = adam_batchsize
       localbatch = adam_batchsize
       if (ppSize>1) then
          call pp_sum(opt_samplesize)
          call pp_sum(batchsize)
       end if
    case default
       write(0,*) "Error: unknown optimization method: " &
                 // trim(opt_method)
    end select

    opt_batchsize  = batchsize
    opt_batchsize_local = localbatch

    opt_nbatch = nint(dble(opt_samplesize)/dble(opt_batchsize))

    ! initialize sequential schedule; for other sampling types this
    ! will be overwritten later
    allocate(opt_schedule(opt_samplesize_local), opt_idx(opt_samplesize_local))
    do i = 1, opt_samplesize_local
       opt_schedule(i) = i
    end do

    ! RNG should have been initialized in train.f90
    call random_reinit()

  end subroutine opt_init_training

  !--------------------------------------------------------------------!

  subroutine opt_final()

    implicit none

    if (.not. isInit) return

    select case(trim(io_lower(opt_method)))
    case('bfgs')
       call opt_final_bfgs()
    case('ekf')
       call opt_final_ekf()
    case('lm')
       call opt_final_lm()
    case('online_sd')
       call opt_final_sd()
    case('adam')
       call opt_final_adam()
    end select

    if (allocated(opt_schedule)) deallocate(opt_schedule, opt_idx)

    opt_memsize = 0
    isInit = .false.

  end subroutine opt_final

  !--------------------------------------------------------------------!
  !              print info about the optimization method              !
  !--------------------------------------------------------------------!

  subroutine opt_print_info()

    implicit none

    if (.not. isInit) return

    write(*,*) "Sampling type               : " &
         // trim(adjustl(opt_sampling_type))

    select case(trim(io_lower(opt_method)))
    case('bfgs')
       continue
    case('ekf')
       continue
    case('lm')
       continue
    case('online_sd')
       continue
    case('adam')
       call opt_print_info_adam()
    end select

  end subroutine opt_print_info

  !--------------------------------------------------------------------!
  !                         sampling schedule                          !
  !                                                                    !
  ! Here, the samples for each epoch are determined.                   !
  !--------------------------------------------------------------------!

  subroutine opt_schedule_epoch(ntrain, energies, errors)

    implicit none

    integer,                             intent(in) :: ntrain
    double precision, dimension(ntrain), intent(in) :: energies
    double precision, dimension(ntrain), intent(in) :: errors

    integer            :: i
    double precision   :: r
    integer, parameter :: a = 4

    select case(trim(adjustl(io_lower(opt_sampling_type))))
    case('sequential')
       ! the default schedule is sequential, so nothing to do
       return
    case('random')
       ! draw samples randomly
       do i = 1, opt_samplesize_local
          call random_number(r)
          opt_schedule(i) = max(ceiling(r*dble(ntrain)), 1)
       end do
    case('weighted')
       ! samples with large errors are learned more frequently
       call argsort(abs(errors), opt_idx)
       do i = 1, opt_samplesize_local
          call random_number(r)
          r = r**a
          opt_schedule(i) = opt_idx(&
               ntrain - max(ceiling(r*dble(ntrain)), 1) + 1)
       end do
    case('energy')
       ! samples with low energies are learned more frequently
       call argsort(energies, opt_idx)
       do i = 1, opt_samplesize_local
          call random_number(r)
          r = r**a
          opt_schedule(i) = opt_idx(max(ceiling(r*dble(ntrain)), 1))
       end do
    case default
       write(0,*) 'Error: Unknown sampling type: ' // trim(opt_sampling_type)
       stop
    end select

    ! Sanity check
    if (any(opt_schedule > opt_samplesize_local)) then
       write(0,*) "Error: invalid schedule: ", opt_schedule(i)
       stop
    end if

  end subroutine opt_schedule_epoch

  !--------------------------------------------------------------------!
  !                         save/restore state                         !
  !--------------------------------------------------------------------!

  subroutine opt_save_state(file, unit)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer :: u

    call opt_assert_moduleinit('save_state')

    if (.not. ppMaster) return

    if (.not. (present(file) .or. present(unit))) then
       write(0,*) "Warning: neither file name nor unit given in `opt_save_state'."
       write(0,*) "         Nothing will be saved."
       return
    end if

    if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='replace', action='write', &
            form='unformatted')
    else
       u = unit
    end if

    write(u) opt_method

    select case(trim(opt_method))
    case('bfgs')
       call opt_save_state_bfgs(u)
    case('ekf')
       call opt_save_state_ekf(u)
    case default
       write(0,*) "Error: this should never have happened! (1)"
       stop
    end select

    if (present(file)) close(u)

  end subroutine opt_save_state

  !--------------------------------------------------------------------!

  subroutine opt_load_state(file, unit)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    character(len=50) :: method
    integer           :: u

    call opt_assert_moduleinit('load_state')

    if (.not. ppMaster) return

    if (.not. (present(file) .or. present(unit))) then
       write(0,*) "Warning: neither file name nor unit given in `opt_load_state'."
       write(0,*) "         Nothing will be loaded."
       return
    end if

    if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='old', action='read', &
            form='unformatted')
       write(*,*) 'Attempting restart of optimization algorithm from file: ' &
                  // trim(file)
    else
       u = unit
    end if

    read(u) method
    if (trim(method) /= trim(opt_method)) then
       write(*,*) 'Optimization method has changed.  Not restarting.'
    else
       select case(trim(opt_method))
       case('bfgs')
          call opt_load_state_bfgs(u)
       case('ekf')
          call opt_load_state_ekf(u)
       case default
          write(0,*) "Error: this should never have happened! (1)"
          stop
       end select
    end if
    write(*,*)

    if (present(file)) close(u)

  end subroutine opt_load_state

  !--------------------------------------------------------------------!
  !      reset state before evaluation of batch (weight training)      !
  !--------------------------------------------------------------------!

  subroutine opt_before_batch()

    implicit none

    call opt_assert_moduleinit('before_batch')

    opt_SSE = 0.0d0

    select case(trim(opt_method))
    case('bfgs')
       bfgs_Dw_sum(:,:) = 0.0d0
    case('ekf')
       continue
    case('lm')
       call opt_before_batch_lm()
    case('online_sd')
       sd_Dw_sum(:,:) = 0.0d0
    case('adam')
       adam_Dw_sum(:,:) = 0.0d0
    end select

  end subroutine opt_before_batch

  !--------------------------------------------------------------------!
  !                       during batch training                        !
  !--------------------------------------------------------------------!

  subroutine opt_after_sample(net, ts, dE, Dw)

    implicit none

    type(Network),    dimension(:),   intent(inout) :: net
    type(TrnSet),                     intent(in)    :: ts
    double precision,                 intent(in)    :: dE
    double precision, dimension(:,:), intent(in)    :: Dw

    call opt_assert_moduleinit('after_sample')

    opt_SSE = opt_SSE + 0.5d0*dE*dE

    select case(trim(opt_method))
    case('bfgs')
       bfgs_Dw_sum(:,:)  = bfgs_Dw_sum(:,:) + dE*Dw(:,:)
    case('ekf')
       call opt_after_sample_ekf(net, ts, dE, Dw)
    case('lm')
       call opt_after_sample_lm(net, ts, dE, Dw)
    case('online_sd')
       call opt_after_sample_sd(net, ts, dE, Dw)
    case('adam')
       call opt_after_sample_adam(dE, Dw)
    end select

  end subroutine opt_after_sample

  !--------------------------------------------------------------------!

  subroutine opt_after_batch(net, ts, do_deriv, do_nextbatch, conv)

    implicit none

    type(Network),    dimension(:),   intent(inout) :: net
    type(TrnSet),                     intent(in)    :: ts
    logical,                          intent(out)   :: do_deriv
    logical,                          intent(out)   :: do_nextbatch
    logical,                          intent(out)   :: conv

    call opt_assert_moduleinit('after_batch')

    ! combine sum of squared errors from all processes
    if (ppSize > 1) call pp_sum(opt_SSE)

    do_deriv = .true.

    select case(trim(opt_method))
    case('bfgs')
       call opt_after_batch_bfgs(net, ts, conv)
       opt_batchiter = opt_batchiter + 1
    case('ekf')
       opt_batchiter = opt_batchiter + 1
       if (ppMaster) call opt_save_state(file=OPT_STATE_FILE)
    case('lm')
       call opt_after_batch_lm(net, ts, do_deriv, conv)
    case('online_sd')
       call opt_after_batch_sd(net, ts, conv)
       opt_batchiter = opt_batchiter + 1
    case('adam')
       call opt_after_batch_adam(net, ts, conv)
       opt_batchiter = opt_batchiter + 1
    case default
       write(0,*) "Error: this should never have happened! (2)"
       stop
    end select

    if (opt_batchiter >= opt_nbatchiters) then
       do_nextbatch = .true.
       opt_batchiter = 1
    end if

  end subroutine opt_after_batch

  !--------------------------------------------------------------------!
  !                       geometry optimization                        !
  !--------------------------------------------------------------------!

  subroutine opt_optimize_coords(E, n, X, F, conv, dmax)

    implicit none

    double precision,                         intent(inout) :: E
    integer,                                  intent(in)    :: n
    double precision, dimension(3,n),         intent(inout) :: X
    double precision, dimension(3,n),         intent(inout) :: F
    logical,                                  intent(out)   :: conv
    double precision, dimension(3), optional, intent(in)    :: dmax

    call opt_assert_moduleinit('optimize_coords')

    select case(trim(opt_method))
    case('bfgs')
       if (present(dmax)) then
          call opt_bfgs_coords(E, n, X, F, conv, dmax=dmax)
       else
          call opt_bfgs_coords(E, n, X, F, conv)
       end if
    case default
       write(0,*) "Error: this should never have happened! (3)"
       stop
    end select

  end subroutine opt_optimize_coords

  !====================================================================!
  !                                                                    !
  !          online steepest descent (error backpropagation)           !
  !                                                                    !
  !====================================================================!

  subroutine opt_init_sd(nw_max, ntypes, methodparam, memsize)

    implicit none

    integer,                        intent(in)  :: nw_max
    integer,                        intent(in)  :: nTypes
    double precision, dimension(:), intent(in)  :: methodparam
    integer,                        intent(out) :: memsize

    if (size(methodparam) /= 2) then
       write(0,*) "Error: wrong number of parameters for " &
               // "the SD method in `opt_init_sd'."
       stop
    end if

    sd_momentum  = methodparam(1)
    sd_learnrate = methodparam(2)

    allocate(sd_Dw_prev(nw_max, ntypes), sd_Dw_sum(nw_max, ntypes))
    sd_Dw_prev(1:nw_max, 1:ntypes) = 0.0d0
    sd_Dw_sum(1:nw_max, 1:ntypes) = 0.0d0
    memsize = 2*nw_max*ntypes*2

  end subroutine opt_init_sd

  !--------------------------------------------------------------------!

  subroutine opt_final_sd()

    implicit none

    if (allocated(sd_Dw_prev)) deallocate(sd_Dw_prev, sd_Dw_sum)

  end subroutine opt_final_sd

  !--------------------------------------------------------------------!

  subroutine opt_after_sample_sd(net, ts, dE, Dw)

    implicit none

    type(Network),    dimension(:),   intent(inout) :: net
    type(TrnSet),                     intent(in)    :: ts
    double precision,                 intent(in)    :: dE
    double precision, dimension(:,:), intent(in)    :: Dw

    sd_Dw_prev(:,:) = sd_momentum*sd_Dw_prev(:,:) &
                    - sd_learnrate*dE*Dw(:,:)

    call opt_update_weights(ts%nTypes, sd_Dw_prev, net)
    sd_Dw_sum(:,:)  = sd_Dw_sum(:,:) + sd_Dw_prev(:,:)

  end subroutine opt_after_sample_sd

  !--------------------------------------------------------------------!

  subroutine opt_after_batch_sd(net, ts, conv)

    implicit none

    type(Network), dimension(:), intent(inout) :: net
    type(TrnSet),                intent(in)    :: ts
    logical,                     intent(out)   :: conv

    if (ppSize>1) then
       ! subtract accumulated weight updates (reset network):
       call opt_update_weights(ts%nTypes, -sd_Dw_sum, net)
       ! gather weights from all processes:
       call pp_sum2d(sd_Dw_sum, opt_nw_max, ts%nTypes)
       ! apply combined weight update:
       call opt_update_weights(ts%nTypes, sd_Dw_sum, net)
    end if

    ! FIXME implement convergence check
    conv = .false.

  end subroutine opt_after_batch_sd

  !====================================================================!
  !                                                                    !
  !              Adaptive Moment Estimation (ADAM) Method              !
  !                                                                    !
  !====================================================================!

  subroutine opt_init_adam(nw_max, ntypes, ntrain_local, methodparam, &
                           memsize)

    implicit none

    integer,                        intent(in)  :: nw_max
    integer,                        intent(in)  :: ntypes
    integer,                        intent(in)  :: ntrain_local
    double precision, dimension(:), intent(in)  :: methodparam
    integer,                        intent(out) :: memsize

    if (size(methodparam) /= 6) then
       write(0,*) "Error: wrong number of parameters for " &
               // "the Adam method in `opt_init_adam'."
       stop
    end if

    adam_learnrate  = methodparam(1)
    adam_b1         = methodparam(2)
    adam_b2         = methodparam(3)
    adam_eps        = methodparam(4)

    ! Process local batch size
    adam_batchsize = max(floor(methodparam(6)/dble(ppSize)), 1)
    adam_batchsize = min(adam_batchsize, ntrain_local)

    ! Process local sample size
    ! make sure that the sample size is a multiple of the batch size
    ! and less than the total number of samples available to the process
    adam_samplesize = max(floor(methodparam(5)/dble(ppSize)), 1)
    adam_samplesize = max(floor(dble(adam_samplesize) &
                        / dble(adam_batchsize)), 1) * adam_batchsize
    adam_samplesize = min(adam_samplesize, ntrain_local)

    allocate(adam_Dw_sum(nw_max, ntypes), adam_Dw_up(nw_max, ntypes), &
             adam_m(nw_max, ntypes), adam_v(nw_max, ntypes), &
             adam_m2(nw_max, ntypes), adam_v2(nw_max, ntypes))
    adam_Dw_sum(1:nw_max, 1:ntypes) = 0.0d0
    adam_m(1:nw_max, 1:ntypes) = 0.0d0
    adam_v(1:nw_max, 1:ntypes) = 0.0d0
    adam_t = 1

    memsize = 6*nw_max*ntypes*2

  end subroutine opt_init_adam

  !--------------------------------------------------------------------!

  subroutine opt_final_adam()

    implicit none

    if (allocated(adam_Dw_sum)) then
       deallocate(adam_Dw_sum, &
                  adam_Dw_up,  &
                  adam_m,      &
                  adam_v,      &
                  adam_m2,     &
                  adam_v2)
    end if

  end subroutine opt_final_adam

  !--------------------------------------------------------------------!

  subroutine opt_print_info_adam()

    implicit none

    write(*,'(1x,"Learning rate               : ",ES9.3)') adam_learnrate
    write(*,*) "Gradient decay rate (b1)    : " // trim(io_adjustl(adam_b1, 4))
    write(*,*) "Sq. gradient decay rate (b1): " // trim(io_adjustl(adam_b2, 4))
    write(*,'(1x,"Epsilon                     : ",ES9.3)') adam_eps
    write(*,*) "Sample size                 : " &
         // trim(io_adjustl(opt_samplesize)) // " (" &
         // trim(io_adjustl(opt_samplesize_local)) // ")"
    write(*,*) "Mini-batch size             : " &
         // trim(io_adjustl(opt_batchsize)) // " (" &
         // trim(io_adjustl(opt_batchsize_local)) // ")"

  end subroutine opt_print_info_adam

  !--------------------------------------------------------------------!

  subroutine opt_after_sample_adam(dE, Dw)

    implicit none

    double precision,                 intent(in)    :: dE
    double precision, dimension(:,:), intent(in)    :: Dw

    ! cumulate mini-batch Jacobian
    adam_Dw_sum(:,:)  = adam_Dw_sum(:,:) + dE*Dw(:,:)

  end subroutine opt_after_sample_adam

  !--------------------------------------------------------------------!

  subroutine opt_after_batch_adam(net, ts, conv)

    implicit none

    type(Network), dimension(:), intent(inout) :: net
    type(TrnSet),                intent(in)    :: ts
    logical,                     intent(out)   :: conv

    if (ppSize>1) then
       ! gather Jacobian from all processes:
       call pp_sum2d(adam_Dw_sum, opt_nw_max, ts%nTypes)
    end if

    ! exponentially decaying averages of the gradient and the squared gradient
    adam_m = adam_b1*adam_m + (1.0d0 - adam_b1)*adam_Dw_sum(:,:)
    adam_v = adam_b2*adam_v + (1.0d0 - adam_b2)*adam_Dw_sum(:,:)**2

    ! renormalize to avoid bias towards zero
    adam_m2 = 1.0d0/(1.0d0 - adam_b1**adam_t)*adam_m
    adam_v2 = 1.0d0/(1.0d0 - adam_b2**adam_t)*adam_v

    ! Adam weight update
    adam_Dw_up(:,:) = -adam_learnrate/(sqrt(adam_v2)+adam_eps)*adam_m2
    call opt_update_weights(ts%nTypes, adam_Dw_up, net)

    ! reset the cumulative Jacobian
    adam_Dw_sum(:,:)  = 0.0d0

    ! increment step
    adam_t = adam_t + 1

    ! FIXME implement convergence check
    conv = .false.

  end subroutine opt_after_batch_adam

  !====================================================================!
  !                                                                    !
  !                       Extended Kalman Filter                       !
  !                                                                    !
  !====================================================================!

  subroutine opt_init_ekf(nw_tot, methodparam, memsize)

    implicit none

    integer,                        intent(in) :: nw_tot
    double precision, dimension(:), intent(in) :: methodparam
    integer,                        intent(out) :: memsize

    integer :: iw, iwg
    logical :: fexists

    if (size(methodparam) /= 6) then
       write(0,*) "Error: wrong number of parameters for " &
               // "the EKF method in `opt_init_ekf'."
       stop
    end if

    ekf_lambda     = methodparam(1)
    ekf_lambda0    = methodparam(2)
    ekf_P_initial  = methodparam(3)
    ekf_pnoise     = methodparam(4)
    ekf_adaptive   = methodparam(5)
    ekf_wg_max     = nint(methodparam(6))

    ! FIXME implement checks to confirm reasonable parameters

    ekf_nwg = ceiling(dble(opt_nw_tot)/dble(ekf_wg_max))

    allocate(ekf_J(nw_tot,ppSize),                 &
             ekf_P(ekf_wg_max,ekf_wg_max,ekf_nwg), &
             ekf_E(ppSize),                        &
             ekf_Wup(nw_tot),                      &
             ekf_KT(ppSize,ekf_wg_max),            &
             A(ppSize,ppSize),                     &
             B(ekf_wg_max,ppSize),                 &
             C(ekf_wg_max,ekf_wg_max),             &
             ekf_wg(ekf_nwg))


    memsize = nw_tot*ppSize*2 + 2*ekf_wg_max*ekf_wg_max*ekf_nwg*2 &
            + ppSize*ekf_wg_max*2 + ppSize*2 + ppSize*ppSize*2    &
            + ekf_wg_max*ppSize*2 + ekf_nwg

    ! determine sizes of weight groups
    iw = nw_tot
    do iwg = 1, ekf_nwg
       ekf_wg(iwg) = min(ekf_wg_max, iw)
       iw = iw - ekf_wg(iwg)
    end do
    if (.not. iw == 0) write(0,*) "Warning: error in weight groups !"

    ! initial state covariance matrix
    do iwg = 1, ekf_nwg
       do iw = 1, ekf_wg(iwg)
          ekf_P(iw,iw,iwg) = 1.0d0/ekf_P_initial
       end do
    end do

    if (ppMaster) then
       inquire(file=OPT_STATE_FILE, exist=fexists)
       if (fexists) call opt_load_state(file=OPT_STATE_FILE)
    end if

  end subroutine opt_init_ekf

  !--------------------------------------------------------------------!

  subroutine opt_final_ekf()

    implicit none

    if (allocated(ekf_J)) then
       deallocate(ekf_J, ekf_P, ekf_KT, ekf_E, ekf_Wup, ekf_wg, A, B, C)
    end if

  end subroutine opt_final_ekf

  !--------------------------------------------------------------------!

  subroutine opt_save_state_ekf(unit)

    implicit none

    integer, intent(in) :: unit

    write(unit) opt_nw_tot, ekf_wg_max, ekf_nwg
    write(unit) ekf_lambda
    write(unit) ekf_wg(1:ekf_nwg)
    write(unit) ekf_P(1:ekf_wg_max,1:ekf_wg_max,1:ekf_nwg)

  end subroutine opt_save_state_ekf

  !--------------------------------------------------------------------!

  subroutine opt_load_state_ekf(unit)

    implicit none

    integer, intent(in) :: unit

    integer :: nw_tot, wg_max, nwg

    read(unit) nw_tot, wg_max, nwg

    if ((nw_tot /= opt_nw_tot) .or. (wg_max /= ekf_wg_max) &
         .or. (nwg /= ekf_nwg)) then
       write(*,*) "Incompatible restart file.  Not restarting."
    else
       write(*,*) "Restarting Extended Kalman Filter."
       read(unit) ekf_lambda
       read(unit) ekf_wg(1:ekf_nwg)
       read(unit) ekf_P(1:ekf_wg_max,1:ekf_wg_max,1:ekf_nwg)
    end if

  end subroutine opt_load_state_ekf

  !--------------------------------------------------------------------!

  subroutine opt_after_sample_ekf(net, ts, dE, Dw)

    implicit none

    type(Network),    dimension(:),   intent(inout) :: net
    type(TrnSet),                     intent(in)    :: ts
    double precision,                 intent(in)    :: dE
    double precision, dimension(:,:), intent(in)    :: Dw

    integer :: iw, nw, it

    ekf_J(:,:) = 0.0d0
    ekf_E(:) = 0.0d0

    ! each process stores its Jacobian ...
    iw = 1
    do it = 1, ts%nTypes
       nw = net(it)%Wsize
       ekf_J(iw:iw+nw-1,ppRank+1) = -Dw(1:nw,it)
       iw = iw + nw
    end do

    ! ... and the corresponding error
    ekf_E(ppRank+1) = dE

    ! combine the information from all processes
    if (ppSize>1) then
       call pp_sum2d(ekf_J, opt_nw_tot, ppSize)
       call pp_sum(ekf_E, ppSize)
    end if

    ! compute and communicate weight updates
    if (ppMaster) call ekf_weight_update()
    if (ppSize>1) call pp_bcast(ekf_Wup, opt_nw_tot)

    ! apply weight updates
    iw = 1
    do it = 1, ts%nTypes
       nw = net(it)%Wsize
       net(it)%W(1:nw) = net(it)%W(1:nw) + ekf_Wup(iw:iw+nw-1)
       iw = iw + nw
    end do

    ekf_lambda = ekf_lambda0*(ekf_lambda - 1.0d0) + 1.0d0

  end subroutine opt_after_sample_ekf

  !--------------------------------------------------------------------!

  subroutine ekf_weight_update()

    implicit none

    integer :: ip, iw, nw, i1, i2, iwg, info

    external :: DSYR2K, DGEMV ! BLAS
    external :: DPOTRF, DPOTRS ! LAPACK

    ! loop over weight groups
    i1 = 1
    do iwg = 1, ekf_nwg
       nw = ekf_wg(iwg)
       i2 = i1 + nw - 1

       !---------------------------------------------------------------!
       ! Kalman gain:                                                  !
       !                                                               !
       ! K = P.J.(J^T.P.J + R)^-1                   (nw x ppSize)      !
       !   = B.(J^T.B + R)^-1  with  B = P.J        (nw x ppSize)      !
       !   = B.A^-1            with  A = J^T.B + R  (ppSize x ppSize)  !
       !                                                               !
       ! P : state covariance matrix (nw x nw_tot)                     !
       ! J : Jacobian (nw x ppSize)                                    !
       ! R : measurement error covariance matrix (ppSize x ppSize)     !
       !---------------------------------------------------------------!

       B(1:nw,1:ppSize) = matmul(ekf_P(1:nw,1:nw,iwg),ekf_J(i1:i2,1:ppSize))

       ! A = 0.5*(J^T.B + B^T.J)
       ! note that we only get the upper triangular A
       call DSYR2K('U', 'T', ppSize, nw, 0.5d0, &
                   ekf_J(i1:i2,1:ppSize), nw, B(1:nw,1:ppSize), &
                   nw, 0.0d0, A, ppSize)

       ! add noise R = l*I to complete A
       ! (if the noise is too small, the Cholesky decomp. might fail)
       do ip = 1, ppSize
          A(ip,ip) = A(ip,ip) + ekf_lambda
          ! alternatively: add fixed amount of measurement noise
          ! A(ip,ip) = A(ip,ip) + ekf_mnoise
       end do

       !---------------------------------------------------------------!
       ! K = B.A^-1 <=> A^T.K^T = B^T   with   A^T = A                 !
       ! solve for K^T!                                                !
       !---------------------------------------------------------------!

!$  if (ppSize > 1) then
       ! Cholesky decomposition of A
       call DPOTRF('U', ppSize, A, ppSize, info)
       if (info /= 0) then
          write(0,*) "Error: Cholesky decomposition failed."
          stop
       end if

       ekf_KT(1:ppSize,1:nw) = transpose(B(1:nw,1:ppSize))
       call DPOTRS('U', ppSize, nw, A, ppSize, ekf_KT(1:ppSize,1:nw), ppSize, info)

       ! weight update = KT^T . E
       ! ekf_Wup(i1:i2) = matmul(transpose(ekf_KT(1:ppSize,1:nw)), ekf_E)
       call DGEMV('T', ppSize, nw, 1.0d0, ekf_KT(1:ppSize,1:nw), &
                  ppSize, ekf_E, 1, 0.0d0, ekf_Wup(i1:i2), 1)
!$  else
!$     KT_ekf = transpose(B/A(1,1))
!$     Wup = B(:,1)/A(1,1)*E_ekf(1)
!$  end if

       !---------------------------------------------------------------!
       ! update state covariance matrix P                              !
       ! P(k+1) = P(k) - K.J^T.P(k) + Q                                !
       !        = P(k) - C^T.P(k) + Q  with  C = J.K^T  (nw x nw)      !
       ! since P is symmetric                                          !
       ! P(k+1) = (I - C^T).P(k) = P(k).(I - C) = - P(k).(C - I)       !
       !---------------------------------------------------------------!

       C(1:nw,1:nw) = matmul(transpose(ekf_KT(1:ppSize,1:nw)), &
                             transpose(ekf_J(i1:i2,1:ppSize)))

       ekf_P(1:nw,1:nw,iwg) = ekf_P(1:nw,1:nw,iwg) &
                            - matmul(C(1:nw,1:nw), ekf_P(1:nw,1:nw,iwg))

       ekf_P(1:nw,1:nw,iwg) = ekf_P(1:nw,1:nw,iwg)/ekf_lambda

!!$       C(1:nw,1:nw) = matmul(ekf_J(i1:i2,1:ppSize), ekf_KT(1:ppSize,1:nw))
!!$
!!$       ! ekf_P = ekf_P - matmul(transpose(C),ekf_P)
!!$       do iw = 1, nw
!!$          C(iw,iw) = C(iw,iw) - 1.0d0
!!$       end do
!!$       ekf_P(1:nw,1:nw,iwg) = -matmul(ekf_P(1:nw,1:nw,iwg),C(1:nw,1:nw))/ekf_lambda

       ! add processing noise
       if (ekf_pnoise > 0.0d0) then
          do iw = 1, nw
             ekf_P(iw,iw,iwg) = ekf_P(iw,iw,iwg) + ekf_pnoise
          end do
       end if

       i1 = i2 + 1
    end do ! weight group

  end subroutine ekf_weight_update

  !====================================================================!
  !                                                                    !
  !                   Levenberg-Marquardt algorithm                    !
  !                                                                    !
  !====================================================================!

  subroutine opt_init_lm(nw_tot, ntrain, methodparam, batchsize, &
                         localbatch, nbatchiters, memsize)

    implicit none

    integer,                        intent(in)  :: nw_tot
    integer,                        intent(in)  :: ntrain
    double precision, dimension(:), intent(in)  :: methodparam
    integer,                        intent(out) :: batchsize
    integer,                        intent(out) :: localbatch
    integer,                        intent(out) :: nbatchiters
    integer,                        intent(out) :: memsize

    if (size(methodparam) /= 5) then
       write(0,*) "Error: wrong number of parameters for LM method in `opt_init_lm'."
       stop
    end if

    batchsize    = min(nint(methodparam(1)), ntrain)
    lm_learnrate = methodparam(2)
    lm_adjust    = methodparam(5)
    lm_stat      = 'jacob'

    localbatch  = nint(batchsize/dble(ppSize))
    batchsize   = ppSize*localbatch
    nbatchiters = nint(methodparam(3))

    allocate(lm_Jw(nw_tot,batchsize), &
             lm_dEv(batchsize),       &
             lm_H(nw_tot, nw_tot),    &
             lm_grad(nw_tot),         &
             lm_Wup(nw_tot),          &
             lm_W_save(nw_tot),       &
             A(nw_tot, nw_tot))

    memsize = nw_tot*batchsize*2 + 2*nw_tot*nw_tot*2 + batchsize*2 &
            + 3*nw_tot*2

  end subroutine opt_init_lm

  !--------------------------------------------------------------------!

  subroutine opt_final_lm()

    implicit none

    if (allocated(lm_Jw))     deallocate(lm_Jw)
    if (allocated(lm_dEv))    deallocate(lm_dEv)
    if (allocated(lm_H))      deallocate(lm_H)
    if (allocated(lm_grad))   deallocate(lm_grad)
    if (allocated(lm_Wup))    deallocate(lm_Wup)
    if (allocated(lm_W_save)) deallocate(lm_W_save)
    if (allocated(A))         deallocate(A)

  end subroutine opt_final_lm

  !--------------------------------------------------------------------!

  subroutine opt_before_batch_lm()

    implicit none

    lm_Jw(:,:) = 0.0d0
    lm_dEv(:) = 0.0d0
    lm_iJw = ppRank*opt_batchsize_local + 1

  end subroutine opt_before_batch_lm

  !--------------------------------------------------------------------!

  subroutine opt_after_sample_lm(net, ts, dE, Dw)

    implicit none

    type(Network),    dimension(:),   intent(in) :: net
    type(TrnSet),                     intent(in) :: ts
    double precision,                 intent(in) :: dE
    double precision, dimension(:,:), intent(in) :: Dw

    integer :: iw, nw, it

    ! store one column of the Jacobian matrix
    iw = 1
    do it = 1, ts%nTypes
       nw = net(it)%Wsize
       lm_Jw(iw:iw+nw-1,lm_iJw) = Dw(1:nw,it)
       iw = iw + nw
    end do

    ! store error in error vector
    lm_dEv(lm_iJw) = dE

    lm_iJw = lm_iJw + 1

  end subroutine opt_after_sample_lm

  !--------------------------------------------------------------------!

  subroutine opt_after_batch_lm(net, ts, do_deriv, conv)

    implicit none

    type(Network), dimension(:), intent(inout) :: net
    type(TrnSet),                intent(in)    :: ts
    logical,                     intent(out)   :: do_deriv
    logical,                     intent(out)   :: conv

    ! FIXME implement convergence test
    conv = .false.

    select case(trim(lm_stat))
    case('jacob')
       call lm_jacobian(net, ts)
       do_deriv     = .false.
    case('error')
       call lm_inner_loop(net, ts, do_deriv)
    end select ! LM_stat

  end subroutine opt_after_batch_lm

  !--------------------------------------------------------------------!

  subroutine lm_jacobian(net, ts)

    implicit none

    type(Network), dimension(:), intent(inout) :: net
    type(TrnSet),                intent(in)    :: ts

    integer :: iw, nw, it, info

    ! FIXME: replace by module-wide interfaces
    external :: DPOTRF, DPOTRS

    ! combine Jacobian from different processes
    if (ppSize > 1) then
       call pp_sum2d(lm_Jw, opt_nw_tot, opt_batchsize)
    end if

    !---------------- Levenberg-Marquard weight update ----------------!
    !                                                                  !
    ! g : gradient | J : Jacobian | H : approx. Hessian | w : weights  !
    ! e : error vector                                                 !
    !                                                                  !
    ! w(k+1) = w(k) - H(k)^-1 . g(k)                                   !
    ! with  H(k) = [J(k)^T.J(k) + l*I]  and  g(k) = J(k)*e(k)          !
    ! <-->  H . [w(k) - w(k+1)] = g(k)                                 !
    ! with g(k)   = J(k) . e(k)                                        !
    !                                                                  !
    ! use LAPACK's DPOTRF & DPOTRS to Factorize and Solve for          !
    ! dw = [w(k) - w(k+1)] --> w(k+1) = w(k) - dw                      !
    !------------------------------------------------------------------!

    if (ppMaster) then
       lm_H(:,:) = matmul(lm_Jw,transpose(lm_Jw))
       do iw = 1, opt_nw_tot
          lm_H(iw,iw) = lm_H(iw,iw) + lm_learnrate
       end do
       lm_grad(:) = matmul(lm_Jw,lm_dEv)
       ! FIXME replace the following by parallel SCALAPACK calls
       A = lm_H
       lm_Wup = lm_grad
       call DPOTRF('U', opt_nw_tot, A, opt_nw_tot, info)
       if (info /= 0) then
          write(0,*) "Error: Cholesky decomposition failed in `lm_jacobian'."
       else
          call DPOTRS('U', opt_nw_tot, 1, A, opt_nw_tot, lm_Wup, &
                      opt_nw_tot, info)
          if (info /= 0) stop "Fatal LAPACK error in `lm_jacobian'"
       end if
    end if ! ppMaster
    if (ppSize>1) call pp_bcast(lm_Wup, opt_nw_tot)

    iw = 1
    do it = 1, ts%nTypes
       nw = net(it)%Wsize
       lm_W_save(iw:iw+nw-1) = net(it)%W(1:nw)
       net(it)%W(1:nw) = net(it)%W(1:nw) - lm_Wup(iw:iw+nw-1)
       iw = iw + nw
    end do

    lm_stat     = 'error'
    lm_cycle    = lm_cycle + 1
    lm_SSE_prev = opt_SSE

  end subroutine lm_jacobian

  !--------------------------------------------------------------------!

  subroutine lm_inner_loop(net, ts, do_deriv)

    implicit none

    type(Network), dimension(:), intent(inout) :: net
    type(TrnSet),                intent(in)    :: ts
    logical,                     intent(out)   :: do_deriv

    double precision :: learnrate_prev
    integer :: iw, nw, it, info

    if ((opt_SSE > lm_SSE_prev) .and. (lm_cycle <= 5)) then

       ! SSE got higher --> change learning rate and try again
       ! however, after 5 tries accept new weights anyway

       learnrate_prev = lm_learnrate
       lm_learnrate = lm_learnrate*lm_adjust
       if (ppMaster) then
          do iw = 1, opt_nw_tot
             lm_H(iw,iw) = lm_H(iw,iw) - learnrate_prev + lm_learnrate
          end do
          ! FIXME replace the following by parallel SCALAPACK calls
          A = lm_H
          lm_Wup = lm_grad ! lm_grad is available from previous cycle
          call DPOTRF('U', opt_nw_tot, A, opt_nw_tot, info)
          if (info /= 0) then
             write(0,*) "Error: Cholesky decomposition failed in `lm_inner_loop'."
          else
             call DPOTRS('U', opt_nw_tot, 1, A, opt_nw_tot, lm_Wup, &
                         opt_nw_tot, info)
             if (info /= 0) stop "Fatal LAPACK error in `lm_inner_loop'"
          end if
       end if
       if (ppSize>1) call pp_bcast(lm_Wup, opt_nw_tot)

       iw = 1
       do it = 1, ts%nTypes
          nw = net(it)%Wsize
          net(it)%W(1:nw) = lm_W_save(iw:iw+nw-1) - lm_Wup(iw:iw+nw-1)
          iw = iw + nw
       end do

       lm_stat  = 'error'
       lm_cycle = lm_cycle + 1
       do_deriv = .false.

    else

       ! only adjust learning rate, if SSE got lower
       if (opt_SSE <= lm_SSE_prev) lm_learnrate = lm_learnrate/lm_adjust

       lm_stat  = 'jacob'
       lm_cycle = 1
       do_deriv = .true.
       opt_batchiter = opt_batchiter + 1

    end if ! SSE got smaller or 5 cycles done

  end subroutine lm_inner_loop

  !====================================================================!
  !                                                                    !
  !                              L-BFGS-B                              !
  !                                                                    !
  !====================================================================!

  subroutine opt_init_bfgs(nparam, ftol, gtol, memsize)

    implicit none

    integer,          intent(in)  :: nparam
    double precision, intent(in)  :: ftol
    double precision, intent(in)  :: gtol
    integer,          intent(out) :: memsize

    logical :: fexists
    integer :: n

    bfgs_n = nparam
    bfgs_m = 20
    bfgs_factr = ftol
    bfgs_pgtol = gtol
    bfgs_task   = 'START'
    bfgs_iprint = -1
    n = (2*bfgs_m + 5)*bfgs_n + 11*bfgs_m*bfgs_m + 8*bfgs_m
    allocate(bfgs_x(bfgs_n), bfgs_g(bfgs_n), bfgs_l(bfgs_n), &
             bfgs_u(bfgs_n), bfgs_nbd(bfgs_n), bfgs_wa(n),   &
             bfgs_iwa(3*bfgs_n))
    bfgs_nbd(1:bfgs_n) = 0
    memsize = 4*bfgs_n*2 + bfgs_n + n*2 + 3*bfgs_n

    inquire(file=OPT_STATE_FILE, exist=fexists)
    if (fexists) call opt_load_state(file=OPT_STATE_FILE)

  end subroutine opt_init_bfgs

  !--------------------------------------------------------------------!

  subroutine opt_final_bfgs()

    implicit none

    if (allocated(bfgs_x))      deallocate(bfgs_x)
    if (allocated(bfgs_g))      deallocate(bfgs_g)
    if (allocated(bfgs_l))      deallocate(bfgs_l)
    if (allocated(bfgs_u))      deallocate(bfgs_u)
    if (allocated(bfgs_nbd))    deallocate(bfgs_nbd)
    if (allocated(bfgs_wa))     deallocate(bfgs_wa)
    if (allocated(bfgs_iwa))    deallocate(bfgs_iwa)
    if (allocated(bfgs_Dw_sum)) deallocate(bfgs_Dw_sum)

  end subroutine opt_final_bfgs

  !--------------------------------------------------------------------!

  subroutine opt_save_state_bfgs(unit)

    implicit none

    integer, intent(in) :: unit

    integer :: n1, n2

    n1 = (2*bfgs_m + 5)*bfgs_n + 11*bfgs_m*bfgs_m + 8*bfgs_m
    n2 = 3*bfgs_n
    write(unit) bfgs_m, bfgs_n
    write(unit) bfgs_factr, bfgs_pgtol  ! usually not restarted
    write(unit) bfgs_task(1:60)
    write(unit) bfgs_iprint
    write(unit) bfgs_l(1:bfgs_n), bfgs_u(1:bfgs_n), bfgs_nbd(1:bfgs_n), &
                bfgs_wa(1:n1), bfgs_iwa(1:n2), bfgs_csave(1:60),        &
                bfgs_lsave(1:4), bfgs_isave(1:44), bfgs_dsave(1:29)

  end subroutine opt_save_state_bfgs

  !--------------------------------------------------------------------!

  subroutine opt_load_state_bfgs(unit)

    implicit none

    integer, intent(in) :: unit

    double precision :: factr, pgtol
    integer :: m, n, n1, n2

    read(unit) m, n
    if ((m == bfgs_m) .and. (n == bfgs_n)) then
       write(*,*) 'Restarting L-BFGS.'
       n1 = (2*bfgs_m + 5)*bfgs_n + 11*bfgs_m*bfgs_m + 8*bfgs_m
       n2 = 3*bfgs_n
       read(unit) factr, pgtol
       read(unit) bfgs_task(1:60)
       read(unit) bfgs_iprint
       read(unit) bfgs_l(1:bfgs_n), bfgs_u(1:bfgs_n), bfgs_nbd(1:bfgs_n), &
               bfgs_wa(1:n1), bfgs_iwa(1:n2), bfgs_csave(1:60),  &
               bfgs_lsave(1:4), bfgs_isave(1:44), bfgs_dsave(1:29)
    else
       write(*,*) 'Incompatible restart file.  Not restarting L-BFGS.'
       ! skip remaining records in case the file is also
       ! used for something else
       read(unit) ! bfgs_factr, bfgs_pgtol
       read(unit) ! bfgs_task
       read(unit) ! bfgs_iprint
       read(unit) ! bfgs_l(:), ...
    end if

  end subroutine opt_load_state_bfgs

  !--------------------------------------------------------------------!

  subroutine opt_after_batch_bfgs(net, ts, conv)

    implicit none

    type(Network), dimension(:), intent(inout) :: net
    type(TrnSet),                intent(in)    :: ts
    logical,                     intent(out)   :: conv

    integer :: it

    if (ppSize > 1) then
       ! gather derivatives from all processes
       call pp_sum2d(bfgs_Dw_sum, opt_nw_max, ts%nTypes)
    end if
    if (ppMaster) then
       call opt_bfgs_weights(net, opt_SSE, bfgs_Dw_sum, conv)
       call opt_save_state(file=OPT_STATE_FILE)
    end if
    if (ppSize > 1) then
       call pp_bcast(conv)
       do it = 1, ts%nTypes
          call pp_bcast_Network(net(it))
       end do
    end if

  end subroutine opt_after_batch_bfgs

  !--------------------------------------------------------------------!

  subroutine opt_bfgs_weights(net, fval, Dw, conv)

    implicit none

    type(Network),    dimension(:),   intent(inout) :: net
    double precision,                 intent(inout) :: fval
    double precision, dimension(:,:), intent(inout) :: Dw
    logical,                          intent(out)   :: conv

    integer :: ix, itype, ntypes, nw

    conv   = .false.
    ntypes = size(net(:))

    ix = 1
    do itype = 1, ntypes
       nw = net(itype)%Wsize
       bfgs_x(ix:ix+nw-1) = net(itype)%W(1:nw)
       bfgs_g(ix:ix+nw-1) = Dw(1:nw, itype)
       ix = ix + nw
    end do

    call setulb(bfgs_n, bfgs_m, bfgs_x, bfgs_l, bfgs_u, bfgs_nbd,    &
         fval, bfgs_g, bfgs_factr, bfgs_pgtol, bfgs_wa, bfgs_iwa,    &
         bfgs_task, bfgs_iprint, bfgs_csave, bfgs_lsave, bfgs_isave, &
         bfgs_dsave)

    select case(io_lower(bfgs_task(1:2)))
    case('st') ! START --> initialization
       continue
    case('fg') ! FG --> new function value and gradient requested
    case('ne') ! NEW_X --> iteration done, next
       call setulb(bfgs_n, bfgs_m, bfgs_x, bfgs_l, bfgs_u, bfgs_nbd,    &
            fval, bfgs_g, bfgs_factr, bfgs_pgtol, bfgs_wa, bfgs_iwa,    &
            bfgs_task, bfgs_iprint, bfgs_csave, bfgs_lsave, bfgs_isave, &
            bfgs_dsave)
    case('co') ! CONV --> converged
       conv = .true.
    case('er') ! ERROR
       write(0,*) trim(adjustl(bfgs_task))
       stop
    end select

    ix = 1
    do itype = 1, ntypes
       nw = net(itype)%Wsize
       net(itype)%W(1:nw) = bfgs_x(ix:ix+nw-1)
       ix = ix + nw
    end do

  end subroutine opt_bfgs_weights

  !--------------------------------------------------------------------!

  subroutine opt_bfgs_coords(E, n, X, F, conv, dmax)

    implicit none

    double precision,                         intent(inout) :: E
    integer,                                  intent(in)    :: n
    double precision, dimension(3,n),         intent(inout) :: X
    double precision, dimension(3,n),         intent(inout) :: F
    logical,                                  intent(out)   :: conv
    double precision, dimension(3), optional, intent(in)    :: dmax

    double precision :: d
    integer :: i, j, ix

    conv   = .false.

    if (3*n /= bfgs_n) then
       write(0,*) "Error: wrong number of coordinates in `opt_optimize()'."
       stop
    end if

    ix = 0
    do j = 1, n
    do i = 1, 3
       ix = ix + 1
       bfgs_x(ix) =  X(i,j)
       bfgs_g(ix) = -F(i,j)
    end do
    end do

    call setulb(bfgs_n, bfgs_m, bfgs_x, bfgs_l, bfgs_u, bfgs_nbd,    &
         E, bfgs_g, bfgs_factr, bfgs_pgtol, bfgs_wa, bfgs_iwa,       &
         bfgs_task, bfgs_iprint, bfgs_csave, bfgs_lsave, bfgs_isave, &
         bfgs_dsave)

    select case(io_lower(bfgs_task(1:2)))
    case('st') ! START --> initialization
       continue
    case('fg') ! FG --> new function value and gradient requested
    case('ne') ! NEW_X --> iteration done, next
       call setulb(bfgs_n, bfgs_m, bfgs_x, bfgs_l, bfgs_u, bfgs_nbd,    &
            E, bfgs_g, bfgs_factr, bfgs_pgtol, bfgs_wa, bfgs_iwa,       &
            bfgs_task, bfgs_iprint, bfgs_csave, bfgs_lsave, bfgs_isave, &
            bfgs_dsave)
    case('co') ! CONV --> converged
       conv = .true.
    case('er') ! ERROR
       write(0,*) trim(adjustl(bfgs_task))
       stop
    end select

    ix = 0
    do j = 1, n
    do i = 1, 3
       ix = ix + 1
       if (present(dmax)) then
          d = bfgs_x(ix) - X(i,j)
          if (abs(d) <= dmax(i)) then
             X(i,j) = bfgs_x(ix)
          else
             X(i,j) = X(i,j) + d/abs(d)*dmax(i)
          end if
       else
          X(i,j) = bfgs_x(ix)
       end if
    end do
    end do

  end subroutine opt_bfgs_coords

  !--------------------------------------------------------------------!
  !                        auxiliary procedures                        !
  !--------------------------------------------------------------------!

  subroutine opt_assert_moduleinit(sub)
    implicit none
    character(len=*), optional, intent(in) :: sub
    if (.not. isInit) then
       if (present(sub)) then
          write(0,*) "Error: module `optimize' not initialized " &
                  // "in `" // trim(sub) // "'."
       else
          write(0,*) "Error: module `optimize' not initialized."
       end if
       if (ppSize>1) write(0,*) "       (rank = ", ppRank, ")"
       stop
    end if
  end subroutine opt_assert_moduleinit

  !--------------------------------------------------------------------!
  !                           weight updates                           !
  !--------------------------------------------------------------------!

  subroutine opt_update_weights(nTypes, Dw, net)

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

  end subroutine opt_update_weights

  !--------------------------------------------------------------------!
  !                          debugging tools                           !
  !--------------------------------------------------------------------!

  subroutine print_matrix(A)

    implicit none

    double precision, dimension(:,:), intent(in) :: A

    integer :: i, j

    do i = 1, size(A(:,1))
       do j = 1, size(A(1,:))
          write(*,'(ES11.4,2x)', advance='no') A(i,j)
       end do
       write(*,*)
    end do

  end subroutine print_matrix

end module optimize
