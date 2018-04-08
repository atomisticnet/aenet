!-----------------------------------------------------------------------
!                   input.f90 - input file handling
!-----------------------------------------------------------------------
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2018 Nongnuch Artrith and Alexander Urban
!+
!+ This Source Code Form is subject to the terms of the Mozilla Public
!+ License, v. 2.0. If a copy of the MPL was not distributed with this
!+ file, You can obtain one at http://mozilla.org/MPL/2.0/.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!+ Mozilla Public License, v. 2.0, for more details.
!-----------------------------------------------------------------------
! 2012-05-27 Nongnuch Artrith (NA), Alexander Urban (AU)
!-----------------------------------------------------------------------

module input

  use aeio,     only: aeio_header,             &
                      aeio_readline,           &
                      aeio_assert_file_exists, &
                      PATHLEN, LINELEN, TYPELEN

  use geometry, only: geo_itype_of_name

  use io,       only: io_adjustl,  &
                      io_lower,    &
                      io_readnext, &
                      io_readval,  &
                      io_unit

  implicit none
  private
  save

  public  :: read_InpGenerate,  &
             read_InpPredict,   &
             read_InpTrain,     &
             read_InpMC,        &
             del_InputData,     &
             inp_read_networks, &
             inp_print_info

  private :: inp_read_value,    &
             inp_read_value_c1, &
             inp_read_value_d1, &
             inp_read_value_i1

  !--------------------------------------------------------------------!
  !                      generic input data type                       !
  !--------------------------------------------------------------------!

  type, public :: InputData

     !-----------------------------------------------------------------!
     ! init             .true. when memory has been allocated          !
     ! file             name of the input file the data was read from  !
     ! outFileName      name of the output file                        !
     ! mode             run mode (predict, mc)                         !
     ! verbosity        verbosity level (0=low, 1=normal, 2=high)      !
     !                                                                 !
     ! nTypes           number of atomic species                       !
     ! typeName(i)      name of i-th atomic species                    !
     ! atomicEnergy(i)  atomic energy of i-th atomic species           !
     !                                                                 !
     ! activeType(i)    == 1, if type i is "active"                    !
     !                                                                 !
     ! netFile(i)       path to NN potential file for species i        !
     ! netArch(i)       architecture (as string) of NN for species i   !
     !                                                                 !
     ! do_forces        .true. if forces shall be calculated           !
     ! do_timing        .true. if timings shall be saved               !
     ! print_atomic_energies  if .true., atomic energies will be saved !
     !                                                                 !
     ! nStrucs          number of structures to run calculations for   !
     ! strucFile(i)     path to i-th atomic structure file             !
     !                                                                 !
     !----------------- structural fingerprint basis ------------------!
     ! setupFile(i)     path to the basis function setup for i         !
     !                                                                 !
     !--------------------------- training ----------------------------!
     ! trn_file         name of the training set file                  !
     ! trn_testset      percentage of datapoints to be used for testing!
     ! trn_maxenergy    max. formation energy per atom for training    !
     ! trn_steps        number of training iterations/epochs           !
     ! trn_sampling     sampling method ('sequential', 'random',       !
     !                  'weighted', or energy)                         !
     ! trn_method       short name of the training method/algorithm    !
     ! trn_methodName   long name of the training method/algorithm     !
     ! trn_nparams      number of parameters of the training method    !
     ! trn_param(i)     i-th parameter of the training algorithm       !
     ! do_save_energies if .true. save testing and training energies   !
     !                  of all structures in the reference data set    !
     !                                                                 !
     !--------------------- structural relaxation ---------------------!
     ! do_relax         .true. if structural relaxation requested      !
     ! relax_method     name of the optimization algorithm             !
     ! relax_steps      max number of relaxation steps                 !
     ! relax_F_conv     convergence criterion for the forces           !
     ! relax_E_conv     convergence criterion for the energy           !
     ! relax_dmax       max. change of coordinate during relax step    !
     !                  (in Angstrom)                                  !
     !                                                                 !
     !-------------------------- Monte-Carlo --------------------------!
     ! nSteps           total nuber of MC steps to be run              !
     ! nSweeps          number of trials per MC step                   !
     ! ensemble         NVT or mVT                                     !
     ! T                temperature                                    !
     ! mu(i)            chemical potential of species i                !
     ! mc_ngroups       number of groups of species                    !
     ! mc_group(i)      group of species i; all species of one group   !
     !                  are exchangeable in an MC simulation           !
     ! mc_ntypes_group(i) number of atom types in group i              !
     ! mc_group_type(i,j) j-th type in i-th group                      !
     !                                                                 !
     !-----------------------------------------------------------------!
     ! do_debug         activate debugging options                     !
     !-----------------------------------------------------------------!

     logical                                             :: init = .false.
     character(len=PATHLEN)                              :: file
     character(len=PATHLEN)                              :: outFileName
     character(len=32)                                   :: mode
     integer                                             :: verbosity

     integer                                             :: nTypes
     character(len=TYPELEN), dimension(:),   allocatable :: typeName
     double precision,       dimension(:),   allocatable :: atomicEnergy

     integer,                dimension(:),   allocatable :: activeType

     character(len=PATHLEN), dimension(:),   allocatable :: netFile
     character(len=LINELEN), dimension(:),   allocatable :: netArch

     logical                                             :: do_forces
     logical                                             :: do_timing
     logical                                             :: print_atomic_energies

     logical                                             :: do_relax
     character(len=32)                                   :: relax_method
     integer                                             :: relax_steps
     double precision                                    :: relax_F_conv
     double precision                                    :: relax_E_conv
     double precision,       dimension(3)                :: relax_dmax

     integer                                             :: nStrucs
     character(len=PATHLEN), dimension(:),   allocatable :: strucFile

     character(len=PATHLEN)                              :: trn_file
     double precision                                    :: trn_testset
     double precision                                    :: trn_maxenergy
     integer                                             :: trn_steps
     character(len=10)                                   :: trn_sampling
     character(len=32)                                   :: trn_method
     character(len=128)                                  :: trn_methodName
     integer                                             :: trn_nparams
     double precision,       dimension(:),   allocatable :: trn_param
     logical                                             :: do_save_energies

     character(len=PATHLEN), dimension(:),   allocatable :: setupFile

     integer                                             :: nSteps
     integer                                             :: nSweeps
     character(len=3)                                    :: ensemble
     double precision                                    :: T
     double precision                                    :: T_final
     logical                                             :: mc_relax_final
     double precision,       dimension(:),   allocatable :: mu
     integer                                             :: mc_ngroups
     integer,                dimension(:),   allocatable :: mc_group
     integer,                dimension(:),   allocatable :: mc_ntypes_group
     integer,                dimension(:,:), allocatable :: mc_group_type

     logical                                             :: do_debug
  end type InputData

  !--------------------------------------------------------------------!

  interface inp_read_value
     module procedure inp_read_value_c1, inp_read_value_d1, &
                      inp_read_value_i1, inp_read_value_l1
  end interface inp_read_value

  !----------------------- return status values -----------------------!

  integer, parameter, private :: S_OK    = 0
  integer, parameter, private :: S_WARN  = 1
  integer, parameter, private :: S_ERROR = 2
  integer, parameter, private :: S_NOT   = 3

contains

  !--------------------------------------------------------------------!
  !                         default parameters                         !
  !--------------------------------------------------------------------!

  subroutine inp_defaults(inp)

    implicit none

    type(InputData), intent(inout) :: inp

    inp%outFileName  = "OUT"

    inp%nTypes       = 0
    inp%nStrucs      = 0
    inp%nSteps       = 0
    inp%nSweeps      = huge(1)
    inp%ensemble     = 'mvt'

    inp%do_forces    = .false.
    inp%do_timing    = .false.
    inp%print_atomic_energies = .false.

    inp%do_relax     = .false.
    inp%relax_method = 'bfgs'
    inp%relax_steps  = 99
    inp%relax_F_conv = 0.01d0
    inp%relax_E_conv = 0.001d0
    inp%relax_dmax   = (/1.0d0, 1.0d0, 1.0d0/)

    inp%trn_testset      = 10.0d0
    inp%trn_steps        = 10
    inp%trn_sampling     = 'sequential'
    inp%trn_method       = 'bfgs'
    inp%trn_methodName   = 'Limited Memory BFGS'
    inp%trn_maxenergy    = huge(1.0d0)
    inp%trn_nparams      = 0
    inp%do_save_energies = .false.

    inp%T            = 0.0d0
    inp%T_final      = 0.0d0

    inp%mc_relax_final = .false.
    inp%mc_ngroups     = 0

    inp%do_debug     = .false.

  end subroutine inp_defaults

  !====================================================================!
  !                                                                    !
  !                 parsers for different input files                  !
  !                                                                    !
  !====================================================================!

  !---------------------------- generate.x ----------------------------!

  function read_InpGenerate(file) result(inp)

    implicit none

    character(len=*), intent(in) :: file
    type(InputData)              :: inp

    integer                :: iline, stat
    integer                :: u

    call aeio_assert_file_exists(file)
    inp%file        = trim(file)
    inp%mode        = 'generate'
    call inp_defaults(inp)
    inp%outFileName = 'refdata.train'

    u = io_unit()
    open(u, file=trim(adjustl(file)), status='old', action='read')

    call inp_read_types_and_energies(inp, u, iline, stat)
    if (stat == S_NOT) then
       write(0,*) "Error: no atom TYPES specified in file `", trim(file), "'"
       return
    end if

    call inp_read_setups(inp, u, iline, stat)
    if (stat == S_NOT) then
       write(0,*) "Error: no basis function SETUPS specified in file `", &
                  trim(file), "'"
       return
    end if

    call inp_read_value(u, 'output', inp%outFileName)
    call inp_read_value(u, 'timing', inp%do_timing)
    call inp_read_value(u, 'debug',  inp%do_debug)
    call inp_read_files(inp, u, iline, countonly=.true.)
    call inp_read_verbosity(inp, u)

    close(u)

    inp%init = .true.

  end function read_InpGenerate

  !----------------------------- train.x ------------------------------!

  function read_InpTrain(file) result(inp)

    implicit none

    character(len=*), intent(in) :: file
    type(InputData)              :: inp

    integer                :: iline
    integer                :: u

    call aeio_assert_file_exists(file)
    inp%file        = trim(file)
    inp%mode        = 'train'
    call inp_defaults(inp)
    inp%outFileName = ''

    u = io_unit()
    open(u, file=trim(adjustl(file)), status='old', action='read')

    call inp_read_method(inp, u, iline)

    call inp_read_value(u, 'trainingset',   inp%trn_file)
    call inp_read_value(u, 'testpercent',   inp%trn_testset)
    call inp_read_value(u, 'iterations',    inp%trn_steps)
    call inp_read_value(u, 'maxenergy',     inp%trn_maxenergy)
    call inp_read_value(u, 'sampling',      inp%trn_sampling)
    call inp_read_value(u, 'timing',        inp%do_timing)
    call inp_read_value(u, 'debug',         inp%do_debug)
    call inp_read_value(u, 'save_energies', inp%do_save_energies)
    call inp_read_verbosity(inp, u)

    close(u)

    inp%init = .true.

  end function read_InpTrain

  !---------------------------- predict.x -----------------------------!

  function read_InpPredict(file) result(inp)

    implicit none

    character(len=*), intent(in) :: file
    type(InputData)              :: inp

    integer                :: iline, stat
    integer                :: u

    call aeio_assert_file_exists(file)
    inp%file = trim(file)
    inp%mode = 'predict'
    call inp_defaults(inp)

    u = io_unit()
    open(u, file=trim(adjustl(file)), status='old', action='read')

    call inp_read_types(inp, u, iline, stat)
    if (stat == S_NOT) then
       write(0,*) "Error: no atom TYPES specified in file `", trim(file), "'"
       return
    end if

    call inp_read_networks(inp, u, iline, stat)
    if (stat == S_NOT) then
       write(0,*) "Error: no NETWORKS specified in file `", trim(file), "'"
       return
    end if

    call inp_read_value(u, 'forces', inp%do_forces)
    call inp_read_value(u, 'timing', inp%do_timing)
    call inp_read_value(u, 'print_atomic_energies', inp%print_atomic_energies)
    call inp_read_value(u, 'debug',  inp%do_debug)
    call inp_read_relax(inp, u, iline)
    call inp_read_files(inp, u, iline)
    call inp_read_verbosity(inp, u)

    close(u)

    inp%init = .true.

  end function read_InpPredict

  !------------------------------- mc.x -------------------------------!

  function read_InpMC(file) result(inp)

    implicit none

    character(len=*), intent(in) :: file
    type(InputData)              :: inp

    integer                :: iline, stat
    integer                :: u
    integer                :: itype

    call aeio_assert_file_exists(file)
    inp%file = trim(file)
    inp%mode = 'mc'
    call inp_defaults(inp)

    u = io_unit()
    open(u, file=trim(adjustl(file)), status='old', action='read')

    call inp_read_types(inp, u, iline, stat)
    if (stat == S_NOT) then
       write(0,*) "Error: no atom TYPES specified in file `", trim(file), "'"
       return
    end if

    call inp_read_networks(inp, u, iline, stat)
    if (stat == S_NOT) then
       write(0,*) "Error: no NETWORKS specified in file `", trim(file), "'"
       return
    end if

    allocate(inp%strucFile(1))
    inp%nStrucs = 1
    call inp_read_value(u, 'structure', inp%strucFile(1), stat)
    if (stat == S_NOT) then
       write(0,*) "Error: no input STRUCTURE specified in file `", trim(file), "'"
       return
    end if

    call inp_read_potentials(inp, u, iline, stat)
    if (stat == S_NOT) then
       ! no chemical potentials in input file --> set all to zero
       allocate(inp%mu(inp%nTypes), inp%activeType(inp%nTypes))
       do itype = 1, inp%nTypes
          inp%mu(itype) = 0.0d0
          inp%activeType(:) = 1
       end do
    end if

    ! if no type groups are specified, all types are in the same group
    call inp_read_mc_groups(inp, u, iline, stat)

    call inp_read_value(u, 'timing',      inp%do_timing)
    call inp_read_value(u, 'steps',       inp%nSteps)
    call inp_read_value(u, 'mctrials',    inp%nSweeps)
    call inp_read_value(u, 'ensemble',    inp%ensemble)
    call inp_read_value(u, 'debug',       inp%do_debug)
    call inp_read_value(u, 'relax_final', inp%mc_relax_final)
    call inp_read_temperature(inp, u, iline, stat)
    call inp_read_relax(inp, u, iline)
    call inp_read_files(inp, u, iline)
    call inp_read_verbosity(inp, u)

    inp%ensemble = io_lower(inp%ensemble)

    if (inp%mc_relax_final) inp%do_relax = .false.

    close(u)

    inp%init = .true.

  end function read_InpMC

  !--------------------------------------------------------------------!

  subroutine del_InputData(inp)

    implicit none

    type(InputData), intent(inout) :: inp

    if (.not. inp%init) return

    if (allocated(inp%typeName))        deallocate(inp%typeName)
    if (allocated(inp%atomicEnergy))    deallocate(inp%atomicEnergy)
    if (allocated(inp%netFile))         deallocate(inp%netFile)
    if (allocated(inp%netArch))         deallocate(inp%netArch)
    if (allocated(inp%strucFile))       deallocate(inp%strucFile)
    if (allocated(inp%trn_param))       deallocate(inp%trn_param)
    if (allocated(inp%setupFile))       deallocate(inp%setupFile)
    if (allocated(inp%mu))              deallocate(inp%mu)
    if (allocated(inp%activeType))      deallocate(inp%activeType)
    if (allocated(inp%mc_group))        deallocate(inp%mc_group)
    if (allocated(inp%mc_ntypes_group)) deallocate(inp%mc_ntypes_group)
    if (allocated(inp%mc_group_type))   deallocate(inp%mc_group_type)

    inp%init = .false.

  end subroutine del_InputData

  !--------------------------------------------------------------------!
  !                         input data summary                         !
  !--------------------------------------------------------------------!

  subroutine inp_print_info(inp)

    implicit none

    type(InputData) :: inp
    integer         :: itype

    if (.not. inp%init) return

    call aeio_header("Input Data Summary")
    write(*,*)

    if (inp%mode == 'mc') then
       write(*,*) 'Monte-Carlo Parameters'
       write(*,*) '----------------------'
       write(*,*)
       write(*,*) 'Temperature : ', trim(io_adjustl(inp%T)), ' K'
       write(*,*) 'MC steps    : ', trim(io_adjustl(inp%nSteps))
       write(*,*)
       write(*,*) 'Chemical Potentials:'
       write(*,*)
       do itype = 1, inp%nTypes
          write(*,'(1x,4x,A2,2x,F10.6)') inp%typeName(itype), inp%mu(itype)
       end do
       write(*,*)
    end if

  end subroutine inp_print_info

  !--------------------------------------------------------------------!
  !     simply read value(s) for specific keyword (implementation)     !
  !--------------------------------------------------------------------!

  subroutine inp_read_value_c1(unit, keyword, dest, stat)

    implicit none

    integer,           intent(in)    :: unit
    character(len=*),  intent(in)    :: keyword
    character(len=*),  intent(inout) :: dest
    integer, optional, intent(out)   :: stat

    integer                :: iline, ipos
    character(len=LINELEN) :: line
    character(len=10)      :: str

    if (present(stat)) stat = S_OK

    call inp_find_keyword(unit, keyword, iline, line=line)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       ipos = 1
       call io_readnext(line, ipos, str)
       call io_readnext(line, ipos, dest)
    end if

  end subroutine inp_read_value_c1

  !--------------------------------------------------------------------!

  subroutine inp_read_value_i1(unit, keyword, dest, stat)

    implicit none

    integer,           intent(in)    :: unit
    character(len=*),  intent(in)    :: keyword
    integer,           intent(inout) :: dest
    integer, optional, intent(out)   :: stat

    integer                :: iline, ipos
    character(len=LINELEN) :: line
    character(len=10)      :: str

    if (present(stat)) stat = S_OK

    call inp_find_keyword(unit, keyword, iline, line=line)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       ipos = 1
       call io_readnext(line, ipos, str)
       call io_readnext(line, ipos, dest)
    end if

  end subroutine inp_read_value_i1

  !--------------------------------------------------------------------!

  subroutine inp_read_value_d1(unit, keyword, dest, stat)

    implicit none

    integer,           intent(in)    :: unit
    character(len=*),  intent(in)    :: keyword
    double precision,  intent(inout) :: dest
    integer, optional, intent(out)   :: stat

    integer                :: iline, ipos
    character(len=LINELEN) :: line
    character(len=10)      :: str

    if (present(stat)) stat = S_OK

    call inp_find_keyword(unit, keyword, iline, line=line)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       ipos = 1
       call io_readnext(line, ipos, str)
       call io_readnext(line, ipos, dest)
    end if

  end subroutine inp_read_value_d1

  !--------------------------------------------------------------------!
  ! Flags can be specified in three different ways in the input file   !
  !                                                                    !
  !   (1) KEYWORD                                                      !
  !   (2) KEYWORD .true.                                               !
  !   (3) KEYWORD .false.                                              !
  !                                                                    !
  ! Both, (1) and (2) will result in the return value '.true.', while  !
  ! (3) will result in '.false.'.  If KEYWORD is not found, the input  !
  ! value of flag will be returned.                                    !
  !--------------------------------------------------------------------!

  subroutine inp_read_value_l1(unit, keyword, flag, stat)

    implicit none

    integer,           intent(in)    :: unit
    character(len=*),  intent(in)    :: keyword
    logical,           intent(inout) :: flag
    integer, optional, intent(out)   :: stat

    integer                :: iline, ipos
    character(len=LINELEN) :: line
    character(len=10)      :: str

    if (present(stat)) stat = S_OK

    call inp_find_keyword(unit, keyword, iline, line=line)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       ipos = 1
       call io_readnext(line, ipos, str)
       call io_readnext(line, ipos, str)
       if (ipos > 0) then
          read(str,*) flag
       else
          flag = .true.
       end if
    end if

  end subroutine inp_read_value_l1

  !--------------------------------------------------------------------!
  !                  procedures for specific keywords                  !
  !--------------------------------------------------------------------!

  subroutine inp_read_files(inp, unit, iline, stat, countonly)

    implicit none

    type(InputData),   intent(inout) :: inp
    integer,           intent(in)    :: unit
    integer,           intent(out)   :: iline
    integer, optional, intent(out)   :: stat
    logical, optional, intent(in)    :: countonly

    logical :: do_strucs
    integer :: i

    if (present(stat)) stat = S_OK
    if (present(countonly)) then
       do_strucs = (.not. countonly)
    else
       do_strucs = .true.
    end if

    call inp_find_keyword(unit, 'files', iline)
    if (iline > 0) then
       call aeio_readline(unit, iline, inp%nStrucs)
       if (do_strucs) then
          allocate(inp%strucFile(inp%nStrucs))
          do i = 1, inp%nStrucs
             call aeio_readline(unit, iline, inp%strucFile(i))
          end do
       end if
    else
       if (present(stat)) stat = S_NOT
    end if

  end subroutine inp_read_files

  !--------------------------------------------------------------------!

  subroutine inp_read_method(inp, unit, iline, stat)

    implicit none

    type(InputData),   intent(inout) :: inp
    integer,           intent(in)    :: unit
    integer,           intent(out)   :: iline
    integer, optional, intent(out)   :: stat

    character(LINELEN) :: line

    if (present(stat)) stat = S_OK

    call inp_find_keyword(unit, 'method', iline)
    if (iline > 0) then
       call aeio_readline(unit, iline, line)
       read(line,*) inp%trn_method
       select case(trim(io_lower(inp%trn_method)))
       case('bfgs')
          inp%trn_methodName = 'Limited Memory BFGS'
       case('ekf')
          inp%trn_methodName = 'Extended Kalman Filter'
          inp%trn_nparams = 6
          allocate(inp%trn_param(inp%trn_nparams))
          inp%trn_param(1) = 0.999d0
          call io_readval(line, 'lambda', inp%trn_param(1))
          inp%trn_param(2) = 0.97d0
          call io_readval(line, 'lambda0', inp%trn_param(2))
          inp%trn_param(3) = 10.0d0
          call io_readval(line, 'state', inp%trn_param(3))
          inp%trn_param(4) = 0.0d0
          call io_readval(line, 'pnoise', inp%trn_param(4))
          inp%trn_param(5) = 0.6d0
          call io_readval(line, 'adaptive', inp%trn_param(5))
          inp%trn_param(6) = 100.0
          call io_readval(line, 'wgroup', inp%trn_param(6))
       case('lm')
          inp%trn_methodName = 'Levenberg-Marquardt'
          inp%trn_nparams = 5
          allocate(inp%trn_param(inp%trn_nparams))
          inp%trn_param(1) = 1000.0d0
          call io_readval(line, 'batchsize', inp%trn_param(1))
          inp%trn_param(2) = 0.1d0
          call io_readval(line, 'learnrate', inp%trn_param(2))
          inp%trn_param(3) = 1.0d0
          call io_readval(line, 'iter', inp%trn_param(3))
          inp%trn_param(4) = 1.0d-3
          call io_readval(line, 'conv', inp%trn_param(4))
          inp%trn_param(5) = 10.0d0
          call io_readval(line, 'adjust', inp%trn_param(5))
       case('online_sd','online_gd')
          inp%trn_methodName = 'Gradient Descent (online)'
          inp%trn_nparams = 2
          allocate(inp%trn_param(inp%trn_nparams))
          inp%trn_param(1) = 0.01d0
          call io_readval(line, 'alpha', inp%trn_param(1))
          inp%trn_param(2) = 0.9d0
          call io_readval(line, 'gamma', inp%trn_param(2))
       case('adam')
          ! mu: learning rate
          ! b1: decay rate for gradient average
          ! b2: decay rate for squared gradient average
          ! eps: small number to prevent division by zero
          ! sample: size of stochastic sample for each epoch
          inp%trn_methodName = 'Adaptive Moment Estimation (Adam)'
          inp%trn_nparams = 6
          allocate(inp%trn_param(inp%trn_nparams))
          inp%trn_param(1) = 0.001d0
          call io_readval(line, 'mu', inp%trn_param(1))
          inp%trn_param(2) = 0.9d0
          call io_readval(line, 'b1', inp%trn_param(2))
          inp%trn_param(3) = 0.999d0
          call io_readval(line, 'b2', inp%trn_param(3))
          inp%trn_param(4) = 1.0d-8
          call io_readval(line, 'eps', inp%trn_param(4))
          inp%trn_param(5) = 500.0d0
          call io_readval(line, 'samplesize', inp%trn_param(5))
          inp%trn_param(6) = 1.0d0
          call io_readval(line, 'batchsize', inp%trn_param(6))
          ! Adam should not be used with sequential sampling
          inp%trn_sampling = 'random'
       end select
    else
       if (present(stat)) stat = S_NOT
    end if

  end subroutine inp_read_method

  !--------------------------------------------------------------------!

  subroutine inp_read_networks(inp, unit, iline, stat, readarch, file)

    implicit none

    type(InputData),            intent(inout) :: inp
    integer,          optional, intent(in)    :: unit
    integer,          optional, intent(out)   :: iline
    integer,          optional, intent(out)   :: stat
    logical,          optional, intent(in)    :: readarch
    character(len=*), optional, intent(in)    :: file

    character(len=LINELEN) :: line
    integer                :: itype
    character(len=32)      :: name
    integer                :: nnets, ipos
    integer                :: u, i, il
    logical                :: do_arch

    if (present(stat)) stat = S_OK

    if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='old', action='read')
    else if (present(unit)) then
       u = unit
    else
       write(0,*) "Error: neither unit nor file given in 'read_networks'."
       stat = S_ERROR
       return
    end if

    if (present(readarch)) then
       do_arch = readarch
    else
       do_arch = .false.
    end if

    call inp_find_keyword(u, 'networks', il)
    if (present(iline)) iline = il
    if (il == 0) then
       if (present(stat)) stat = S_NOT
    else if (inp%nTypes < 1) then
       write(0,*) "Error: atom types have to be defined before networks."
       if (present(stat)) stat = S_ERROR
    else
       if (do_arch) then
          allocate(inp%netFile(inp%nTypes), inp%netArch(inp%nTypes))
       else
          allocate(inp%netFile(inp%nTypes))
       end if
       nnets = inp%nTypes
       do i = 1, nnets
          call aeio_readline(u, il, line)
          ipos = 1
          call io_readnext(line, ipos, name)
          itype = geo_itype_of_name(name, inp%typeName)
          call io_readnext(line, ipos, inp%netFile(itype))
          if (do_arch) inp%netArch(itype) = line(ipos:len_trim(line))
       end do ! inet
    end if

    if (present(file)) close(u)

  end subroutine inp_read_networks

  !--------------------------------------------------------------------!

  subroutine inp_read_types(inp, unit, iline, stat)

    implicit none

    type(InputData),   intent(inout) :: inp
    integer,           intent(in)    :: unit
    integer,           intent(out)   :: iline
    integer, optional, intent(out)   :: stat

    integer :: i

    if (present(stat)) stat = S_OK

    call inp_find_keyword(unit, 'types', iline)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       call aeio_readline(unit, iline, inp%nTypes)
       allocate(inp%typeName(inp%nTypes))
       do i = 1, inp%nTypes
          call aeio_readline(unit, iline, inp%typeName(i))
       end do
    end if

  end subroutine inp_read_types

  !--------------------------------------------------------------------!

  subroutine inp_read_types_and_energies(inp, unit, iline, stat)

    implicit none

    type(InputData),   intent(inout) :: inp
    integer,           intent(in)    :: unit
    integer,           intent(out)   :: iline
    integer, optional, intent(out)   :: stat

    integer :: i
    character(len=LINELEN) :: line

    if (present(stat)) stat = S_OK

    call inp_find_keyword(unit, 'types', iline)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       call aeio_readline(unit, iline, inp%nTypes)
       allocate(inp%typeName(inp%nTypes), &
                inp%atomicEnergy(inp%nTypes))
       do i = 1, inp%nTypes
          call aeio_readline(unit, iline, line)
          read(line,*) inp%typeName(i), inp%atomicEnergy(i)
       end do
    end if

  end subroutine inp_read_types_and_energies

  !--------------------------------------------------------------------!

  subroutine inp_read_potentials(inp, unit, iline, stat)

    implicit none

    type(InputData),   intent(inout) :: inp
    integer,           intent(in)    :: unit
    integer,           intent(out)   :: iline
    integer, optional, intent(out)   :: stat

    integer                :: i, itype, ipos
    character(len=LINELEN) :: line, kwd
    character(len=10)      :: name

    if (present(stat)) stat = S_OK

    call inp_find_keyword(unit, 'potentials', iline)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       allocate(inp%mu(inp%nTypes), inp%activeType(inp%nTypes))
       do i = 1, inp%nTypes
          call aeio_readline(unit, iline, line)
          ipos = 1
          call io_readnext(line, ipos, name)
          itype = geo_itype_of_name(name, inp%typeName)
          call io_readnext(line, ipos, kwd)
          if (trim(io_lower(kwd)) == 'fix') then
             inp%activeType(itype) = 0
             inp%mu(itype) = 0.0d0
          else
             inp%activeType(itype) = 1
             read(kwd, *) inp%mu(itype)
          end if
       end do
    end if

  end subroutine inp_read_potentials

  !--------------------------------------------------------------------!

  subroutine inp_read_mc_groups(inp, unit, iline, stat)

    implicit none

    type(InputData),   intent(inout) :: inp
    integer,           intent(in)    :: unit
    integer,           intent(out)   :: iline
    integer, optional, intent(out)   :: stat

    integer                :: ipos, ig
    integer                :: it, itype
    character(len=LINELEN) :: line
    character(len=10)      :: name

    if (present(stat)) stat = S_OK

    allocate(inp%mc_group(inp%nTypes))

    call inp_find_keyword(unit, 'groups', iline)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
       ! all types are in the same MC group
       inp%mc_group(:) = 1
       inp%mc_ngroups  = 1
       allocate(inp%mc_ntypes_group(inp%mc_ngroups), &
                inp%mc_group_type(inp%mc_ngroups,inp%nTypes))
       inp%mc_ntypes_group(1) = inp%nTypes
       do itype = 1, inp%nTypes
          inp%mc_group_type(1,itype) = itype
       end do
    else
       call aeio_readline(unit, iline, line)
       read(line, *) inp%mc_ngroups
       allocate(inp%mc_ntypes_group(inp%mc_ngroups), &
                inp%mc_group_type(inp%mc_ngroups,inp%nTypes))
       inp%mc_group(:) = 0
       do ig = 1, inp%mc_ngroups
          call aeio_readline(unit, iline, line)
          ipos = 1
          ! first entry on line is the number of typenames to follow
          call io_readnext(line, ipos, name)
          read(name, *) inp%mc_ntypes_group(ig)
          do it = 1, inp%mc_ntypes_group(ig)
             call io_readnext(line, ipos, name)
             itype = geo_itype_of_name(name, inp%typeName)
             inp%mc_group_type(ig,it) = itype
             if (inp%mc_group(itype) /= 0) then
                write(0,*) "Error: Species ", &
                           trim(adjustl(inp%typeName(itype))), &
                           " can not be in multiple MC groups."
                stop
             end if
             inp%mc_group(itype) = ig
          end do
          ! single species groups musty be inactive
          if (inp%mc_ntypes_group(ig) == 1) then
             itype = inp%mc_group_type(ig,1)
             if (inp%activeType(itype) /= 0) then
                write(0,*) 'Warning: setting status of species ', &
                     trim(adjustl(inp%typeName(itype))), " to `inactive'."
                inp%activeType(itype) = 0
                inp%mu(itype) = 0.0d0
             end if
          end if
       end do
    end if

  end subroutine inp_read_mc_groups

  !--------------------------------------------------------------------!

  subroutine inp_read_relax(inp, unit, iline, stat)

    implicit none

    type(InputData),   intent(inout) :: inp
    integer,           intent(in)    :: unit
    integer,           intent(out)   :: iline
    integer, optional, intent(out)   :: stat

    character(len=LINELEN) :: line

    if (present(stat)) stat = S_OK

    call inp_find_keyword(unit, 'relax', iline)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       inp%do_relax  = .true.
       inp%do_forces = .true.
       call aeio_readline(unit, iline, line)
       call io_readval(line, 'method', inp%relax_method)
       call io_readval(line, 'F_conv', inp%relax_F_conv)
       call io_readval(line, 'E_conv', inp%relax_E_conv)
       call io_readval(line, 'steps',  inp%relax_steps)
       call io_readval(line, 'dmax',   inp%relax_dmax, 3)
    end if

  end subroutine inp_read_relax

  !--------------------------------------------------------------------!

  subroutine inp_read_setups(inp, unit, iline, stat)

    implicit none

    type(InputData),   intent(inout) :: inp
    integer,           intent(in)    :: unit
    integer,           intent(out)   :: iline
    integer, optional, intent(out)   :: stat

    character(len=LINELEN) :: line
    integer                :: itype1
    character(len=32)      :: name1
    integer                :: nsetups, ipos
    integer                :: i

    if (present(stat)) stat = S_OK

    call inp_find_keyword(unit, 'setups', iline)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else if (inp%nTypes < 1) then
       write(0,*) "Error: atom types have to be defined before networks."
       if (present(stat)) stat = S_ERROR
    else
       allocate(inp%setupFile(inp%nTypes))
       nsetups = inp%nTypes
       do i = 1, nsetups
          call aeio_readline(unit, iline, line)
          ipos = 1
          call io_readnext(line, ipos, name1)
          itype1 = geo_itype_of_name(name1, inp%typeName)
          call io_readnext(line, ipos, inp%setupFile(itype1))
       end do ! inet
    end if

  end subroutine inp_read_setups

  !--------------------------------------------------------------------!

  subroutine inp_read_temperature(inp, unit, iline, stat)

    implicit none

    type(InputData),   intent(inout) :: inp
    integer,           intent(in)    :: unit
    integer,           intent(out)   :: iline
    integer, optional, intent(out)   :: stat

    character(len=LINELEN) :: line, str
    integer                :: ipos

    if (present(stat)) stat = S_OK

    call inp_find_keyword(unit, 'temperature', iline, line=line)
    if (iline == 0) then
       if (present(stat)) stat = S_NOT
    else
       ipos = 1
       call io_readnext(line, ipos, str)
       call io_readnext(line, ipos, inp%T)
       call io_readnext(line, ipos, str)
       if (ipos > 0) then
          read(str,*) inp%T_final
       else
          inp%T_final = inp%T
       end if
    end if

  end subroutine inp_read_temperature

  !--------------------------------------------------------------------!

  subroutine inp_read_verbosity(inp, unit)

    implicit none

    type(InputData),   intent(inout) :: inp
    integer,           intent(in)    :: unit

    character(len=LINELEN) :: str

    str = 'normal'
    call inp_read_value(unit, 'verbosity', str)
    select case(trim(adjustl(io_lower(str))))
    case('low')
       inp%verbosity = 0
    case('normal')
       inp%verbosity = 1
    case('high')
       inp%verbosity = 2
    end select

  end subroutine inp_read_verbosity

  !--------------------------------------------------------------------!
  !                         read until keyword                         !
  !--------------------------------------------------------------------!

  subroutine inp_find_keyword(unit, keyword, iline, line)

    implicit none

    integer,                          intent(in)  :: unit
    character(len=*),                 intent(in)  :: keyword
    integer,                          intent(out) :: iline
    character(len=LINELEN), optional, intent(out) :: line

    integer                :: stat
    character(len=LINELEN) :: kwd, ln

    rewind(unit)
    kwd = ''
    iline = 0
    skip : do
       call aeio_readline(unit, iline, ln, stat)
       if (stat /= 0) then
          ! end of file and keyword not found
          iline = 0
          exit skip
       end if
       iline = iline + 1
       read(ln, *) kwd
       if (trim(io_lower(kwd))==trim(io_lower(keyword))) then
          ! keyword found
          if (present(line)) line = ln
          exit skip
       end if
    end do skip

  end subroutine inp_find_keyword

end module input
