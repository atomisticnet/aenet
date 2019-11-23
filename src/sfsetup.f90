!-----------------------------------------------------------------------
!     sfsetup.f90 - representation of the local atomic environment
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
! 2011-10-20 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

module sfsetup

  use aeio,     only: aeio_readline, &
                      TYPELEN

  use io,       only: io_lower,      &
                      io_readval,    &
                      io_adjustl,    &
                      io_unit

  use sfbasis,  only: FingerprintBasis, &
                      new_SFBasis,      &
                      del_SFBasis,      &
                      sfb_eval

  use symmfunc, only: sf_init,       &
                      sf_final,      &
                      sf_add_rad,    &
                      sf_add_ang,    &
                      sf_fingerprint

  implicit none
  private
  save

  public  :: read_Setup_parameters, &
             load_Setup,            &
             load_Setup_ASCII,      &
             save_Setup,            &
             save_Setup_ASCII,      &
             skip_Setup,            &
             new_Setup,             &
             del_Setup,             &
             stp_init,              &
             stp_final,             &
             stp_get_range,         &
             stp_eval,              &
             stp_normalize,         &
             stp_nsf_max,           &
             stp_print_info,        &
             stp_assert_init,       &
             stp_assert_moduleinit

  private :: read_basis_chebyshev,       &
             setup_basis_chebyshev,      &
             print_info_chebyshev,       &
             read_symmfunc_Behler2011,   &
             setup_symmfunc_Behler2011,  &
             print_info_Behler2011

  !--------------------------------------------------------------------!
  !                    structural fingerprint setup                    !
  !--------------------------------------------------------------------!

  type, public :: Setup

     !-----------------------------------------------------------------!
     ! init          .true., if the setup has been initialized         !
     ! neval         number of evaluations                             !
     ! description   an optional description from the setup file       !
     ! atomtype      species of the central atom                       !
     ! nenv          number of different surrounding species           !
     ! envtypes      species of the surrounding atoms                  !
     !                                                                 !
     ! The global atom type index is determined by the order of type   !
     ! names in the input file.  The *local* type index is given by    !
     ! the order the basis functions were first set up.  For the       !
     ! evaluation of the basis functions, the types from the input     !
     ! structure (global index) need to be assigned to local type IDs. !
     !                                                                 !
     ! gtype(i)      global atom type ID for local type i              !
     ! ltype(i)      local atom type ID for global type i              !
     !                                                                 !
     ! Rc_min/max    minimal interaction radius and max. cutoff        !
     ! sftype        basis function type (e.g. Behler2011)             !
     ! nsf           number of structural fingerprint basis functions  !
     ! sf(i)         function kind of the particular basis type        !
     ! nsfparam      the max. number of parameters of a basis function !
     ! sfparam(i,j)  i-th parameter of the j-th basis function         !
     !               i <= nsfparam                                     !
     ! sfenv(i,j)    i-th environment species for j-th basis function  !
     !                                                                 !
     ! sfval_min(i)  lowest so far encountered value of the i-th SF    !
     ! sfval_max(i)  largest so far encountered value of the i-th SF   !
     ! sfval_avg(i)  current average value of the i-th symm. function  !
     ! sfval_cov(i)  current covariance of the i-th symm. function     !
     ! --> min/max/avg/cov are updated whenever UNSCALED evaluation    !
     !     is requested.  This is usually during the screening of the  !
     !     training set.  During the training and prediciton a scaled  !
     !     value of the SFs is useful, that lies within [-1:1].        !
     !                                                                 !
     ! The scaling is implemented as                                   !
     !                                                                 !
     !  f(i) = max(sval_avg(i)-sfval_min(i),sfval_max(i)-sval_avg(i))  !
     !  s(i) = sfval_min(i)                                            !
     !  sfval(i) = (sfval(i)-s(i))/f(i)                                !
     !                                                                 !
     !-----------------------------------------------------------------!

     logical                                             :: init
     integer                                             :: neval

     character(len=1024)                                 :: description

     character(len=TYPELEN)                              :: atomtype
     integer                                             :: nenv
     character(len=TYPELEN), dimension(:),   allocatable :: envtypes

     integer                                             :: ntypes_global
     integer,                dimension(:),   allocatable :: gtype
     integer,                dimension(:),   allocatable :: ltype

     double precision                                    :: Rc_min
     double precision                                    :: Rc_max

     character(len=100)                                  :: sftype
     integer                                             :: nsf
     integer,                dimension(:),   allocatable :: sf
     integer                                             :: nsfparam
     double precision,       dimension(:,:), allocatable :: sfparam
     integer,                dimension(:,:), allocatable :: sfenv

     double precision,       dimension(:),   allocatable :: sfval_min
     double precision,       dimension(:),   allocatable :: sfval_max
     double precision,       dimension(:),   allocatable :: sfval_avg
     double precision,       dimension(:),   allocatable :: sfval_cov

  end type Setup

  !--------------- memory for basis function operations ---------------!
  ! sfval(i)         value of the i-th basis function                  !
  ! sfderiv_i(i,j)   i-th component of the derivative of the j-th SF   !
  !                  with respect to the central atom                  !
  !                  sfderiv_i(3,nsf_max)                              !
  ! sfderiv_j(i,j,k) i-th component of the derivative of the j-th SF   !
  !                  with respect to the coordinates of atom k         !
  !                  sfderiv_j(3,nsf_max,nnb_max)                      !
  !--------------------------------------------------------------------!

  integer,                                         public :: nsf_max
  integer,                                         public :: nnb_max

  !------------------------- Chebyshev basis --------------------------!
  ! sfb(i)      structural fingerprint basis of atom type i            !
  !--------------------------------------------------------------------!

  type(FingerprintBasis), dimension(:), allocatable, private :: sfb

  !--------------------------------------------------------------------!

  logical, private :: isInit = .false.

  !---------------------------- constants -----------------------------!
  ! NSFPARAM    maximum  number of SF parameters                       !
  ! NENV_MAX    maximum number of types involved in single function    !
  !             (e.g., 1 = distance, 2 = angle, 3 = dihedral)          !
  !--------------------------------------------------------------------!

  integer, parameter :: NSFPARAM = 4
  integer, parameter :: NENV_MAX = 2

contains

  function read_Setup_parameters(inFile, global_types) result(stp)

    implicit none

    character(len=*),               intent(in) :: inFile
    character(len=*), dimension(:), intent(in) :: global_types
    type(Setup)                                :: stp

    integer,           parameter :: nc = 1024
    character(len=nc)            :: line
    integer                      :: stat, iline
    character(len=20)            :: keyword
    integer                      :: isf

    integer                      :: ntypes_global
    integer                      :: u_stp, i

    stp%Rc_min = 0.0d0
    stp%Rc_max = 0.0d0
    stp%neval  = 0
    stp%nsf = 0

    ntypes_global = size(global_types(:))

    u_stp = io_unit()
    open(u_stp, file=trim(inFile), status='old', action='read')
    iline = 0
    do
       call aeio_readline(u_stp, iline, line, stat)
       if (stat /= 0) exit

       read(line,*) keyword
       select case(io_lower(keyword))
       case('atom')
          read(line(5:nc), *) stp%atomtype
       case('env')
          read(line(4:nc), *) stp%nenv
          allocate(stp%envtypes(stp%nenv))
          do i = 1, stp%nenv
             call aeio_readline(u_stp, iline, line, stat)
             read(line,*) stp%envtypes(i)
          end do
       case('descr')
          stp%description = ' '
          read(u_stp, '(A)') line
          iline = iline + 1
          do while(index(io_lower(line),'end descr') <= 0)
             if (len_trim(stp%description) < len(stp%description)) then
                stp%description = trim(stp%description) // trim(line) // "$$"
             end if
             read(u_stp, '(A)') line
             iline = iline + 1
          end do
       case('rmin')
          read(line(5:nc), *) stp%Rc_min
       case('symmfunc', 'functions')
          if (stp%nsf > 0) then
             write(0,*) "Error: Multiple basis definitions found in " // &
                        "structural fingerprint setup."
             stop
          end if
          call io_readval(line, 'type', stp%sftype)
          read(u_stp, *) stp%nsf
          stp%Rc_max   =    0.0d0
          stp%nsfparam = NSFPARAM
          allocate(stp%sf(stp%nsf), stp%sfparam(NSFPARAM,stp%nsf), &
                   stp%sfval_min(stp%nsf), stp%sfval_max(stp%nsf), &
                   stp%sfval_avg(stp%nsf), stp%sfval_cov(stp%nsf), &
                   stp%sfenv(NENV_MAX,stp%nsf))
          stp%sf(:)        = 0
          stp%sfparam(:,:) = 0.0d0
          stp%sfenv(:,:)   = 0
          stp%sfval_min(:) = 0.0d0
          stp%sfval_max(:) = 0.0d0
          stp%sfval_avg(:) = 0.0d0
          stp%sfval_cov(:) = 0.0d0
          do isf = 1, stp%nsf
             select case(trim(stp%sftype))
             case('Behler2011')
                call read_symmfunc_Behler2011(u_stp, isf, stp, iline)
             case default
                write(0,*) "Error: Unknown basis function type: ", &
                     trim(stp%sftype)
                deallocate(stp%sf, stp%sfparam)
                stop
             end select
          end do
       case('basis')
          if (stp%nsf > 0) then
             write(0,*) "Error: Multiple basis definitions found in " // &
                        "structural fingerprint setup."
             stop
          end if
          call io_readval(line, 'type', stp%sftype)
          select case(io_lower(stp%sftype))
          case('chebyshev')
             call read_basis_chebyshev(u_stp, stp, iline)
          case default
             write(0,*) "Error: Unknown basis type: ", trim(stp%sftype)
             deallocate(stp%sf, stp%sfparam)
             stop
          end select
       case default
          write(0,*) "Error: Unknown keyword : ", trim(keyword)
          write(0,*) "       file            : ", trim(inFile)
          write(0,*) "       line            : ", trim(io_adjustl(iline))
       end select

    end do
    close(u_stp)

    if (stp%Rc_min == 0.0d0) then
       write(0,*) "Warning: RMIN not set in : ", trim(inFile)
       write(0,*) "         --> setting Rc_min = 1.0"
       stp%Rc_min = 1.0d0
    end if

    if (stp%Rc_max == 0.0d0) then
       write(0,*) "Error: Cutoff undetermined not set in : ", trim(inFile)
       stop
    end if

    ! connect local atom type IDs with global ones
    allocate(stp%gtype(stp%nenv), stp%ltype(ntypes_global))
    call stp_set_global_types(stp, ntypes_global, global_types)

    stp%init = .true.

  end function read_Setup_parameters

  !--------------------------------------------------------------------!

  function new_Setup(nsf, nenv, ntypes_global) result(stp)

    implicit none

    integer, intent(in) :: nsf, nenv, ntypes_global

    type(Setup) :: stp

    allocate(stp%envtypes(nenv),        &
             stp%sf(nsf),               &
             stp%sfparam(NSFPARAM,nsf), &
             stp%sfenv(NENV_MAX,nsf),   &
             stp%sfval_min(nsf),        &
             stp%sfval_max(nsf),        &
             stp%sfval_avg(nsf),        &
             stp%sfval_cov(nsf),        &
             stp%gtype(nenv),          &
             stp%ltype(ntypes_global))

    stp%description = ' '
    stp%atomtype    = ' '
    stp%nenv        = nenv
    stp%envtypes    = ' '
    stp%Rc_min      = -1.0d0
    stp%Rc_max      = -1.0d0
    stp%sftype      = ' '

    stp%ntypes_global = ntypes_global
    stp%ltype         = 0
    stp%gtype         = 0

    stp%nsf         = nsf
    stp%nsfparam    = NSFPARAM

    stp%sf(:)        = 0
    stp%sfparam(:,:) = 0.0d0
    stp%sfenv(:,:)   = 0
    stp%sfval_min(:) = 0.0d0
    stp%sfval_max(:) = 0.0d0
    stp%sfval_avg(:) = 0.0d0
    stp%sfval_cov(:) = 0.0d0

    stp%neval        = 0
    stp%init         = .true.

  end function new_Setup

  !--------------------------------------------------------------------!

  subroutine del_Setup(stp)

    implicit none

    type(Setup), intent(inout) :: stp

    if (stp%init) then
       deallocate(stp%sf, stp%sfparam, stp%sfval_min, stp%sfval_max, &
                  stp%sfval_avg, stp%sfval_cov, stp%envtypes, stp%sfenv, &
                  stp%gtype, stp%ltype)
       stp%init = .false.
    end if

  end subroutine del_Setup

  !--------------------------------------------------------------------!
  !           print out information about a specific set-up            !
  !--------------------------------------------------------------------!

  subroutine stp_print_info(stp)

    implicit none

    type(Setup), intent(in) :: stp

    integer :: i

    call stp_assert_init(stp)

    write(*,*) 'Structural fingerprint (SF) set-up for ', &
         trim(adjustl(stp%atomtype))
    write(*,*)

    call stp_print_descr(stp%description)

    if (stp%nenv>0) then
       write(*,'(1x,"environment types: ")', advance='no')
       do i = 1, stp%nenv
          write(*,'(A2,1x)', advance='no') stp%envtypes(i)
       end do
    end if
    write(*,*)
    write(*,*) 'minimal distance : ', trim(io_adjustl(stp%Rc_min,2)), ' Angstrom'
    write(*,*) 'maximal cut-off  : ', trim(io_adjustl(stp%Rc_max,2)), ' Angstrom'
    write(*,*) 'size of basis    : ', trim(io_adjustl(stp%nsf))
    write(*,*) 'evaluations      : ', trim(io_adjustl(stp%neval))
    write(*,*)

    select case(trim(io_lower(stp%sftype)))
    case('chebyshev')
       call print_info_chebyshev(stp)
    case('behler2011')
       call print_info_Behler2011(stp)
    end select


  end subroutine stp_print_info

  !--------------------------------------------------------------------!

  subroutine stp_print_descr(descr)

    implicit none

    character(len=*), intent(in) :: descr

    integer :: i1, i2, l

    l = len_trim(descr)

    i1 = 1
    i2 = index(descr,'$$')
    do while(i2 >= i1)
       write(*,*) descr(i1:i2-1)
       i1 = i2 + 2
       i2 = i1 + index(descr(i1:l),'$$') - 1
    end do
    if (len_trim(descr(i1:l)) > 0) then
       write(*,*) trim(descr(i1:l))
    end if
    write(*,*)

  end subroutine stp_print_descr

  !--------------------------------------------------------------------!
  !             save setup to file / load setup from file              !
  !--------------------------------------------------------------------!

  subroutine save_Setup(stp, file, unit)

    implicit none

    type(Setup),                intent(in) :: stp
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer :: u

    call stp_assert_init(stp)

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='replace', action='write', &
            form='unformatted')
    else
       write(0,*) "Error: neither unit number nor file name given"
       write(0,*) "       in `save_Setup'."
       return
    end if

    write(u) stp%description
    write(u) stp%atomtype
    write(u) stp%nenv
    write(u) stp%envtypes(:)
    write(u) stp%Rc_min
    write(u) stp%Rc_max
    write(u) stp%sftype
    write(u) stp%nsf
    write(u) stp%nsfparam
    write(u) stp%sf(:)
    write(u) stp%sfparam(:,:)
    write(u) stp%sfenv(:,:)
    write(u) stp%neval
    write(u) stp%sfval_min(:)
    write(u) stp%sfval_max(:)
    write(u) stp%sfval_avg(:)
    write(u) stp%sfval_cov(:)

    if (.not. present(unit)) then
       close(u)
    end if

  end subroutine save_Setup

  !--------------------------------------------------------------------!

  subroutine save_Setup_ASCII(stp, file, unit)

    implicit none

    type(Setup),                intent(in) :: stp
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer :: i, j
    integer :: u

    character(len=*), parameter :: dfrmt = '(4(1x,ES24.17))'
    character(len=*), parameter :: ifrmt = '(4(1x,I17))'

    call stp_assert_init(stp)

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='replace', action='write')
    else
       write(0,*) "Error: neither unit number nor file name given"
       write(0,*) "       in `save_Setup'."
       return
    end if

    write(u,'(A)') trim(stp%description)
    write(u,'(A)') stp%atomtype
    write(u,*) stp%nenv
    write(u,'(A)') (stp%envtypes(i), i=1,stp%nenv)
    write(u,*) stp%Rc_min
    write(u,*) stp%Rc_max
    write(u,'(A)') stp%sftype
    write(u,*) stp%nsf
    write(u,*) stp%nsfparam
    write(u,ifrmt) (stp%sf(i), i=1,stp%nsf)
    write(u,dfrmt) ((stp%sfparam(i,j), i=1,stp%nsfparam), j=1,stp%nsf)
    write(u,ifrmt) ((stp%sfenv(i,j), i=1,stp%nenv), j=1,stp%nsf)
    write(u,*) stp%neval
    write(u,dfrmt) (stp%sfval_min(i), i=1,stp%nsf)
    write(u,dfrmt) (stp%sfval_max(i), i=1,stp%nsf)
    write(u,dfrmt) (stp%sfval_avg(i), i=1,stp%nsf)
    write(u,dfrmt) (stp%sfval_cov(i), i=1,stp%nsf)

    if (.not. present(unit)) close(u)

  end subroutine save_Setup_ASCII

  !--------------------------------------------------------------------!

  function load_Setup(global_types, file, unit) result(stp)

    implicit none

    character(len=*), dimension(:), intent(in) :: global_types
    character(len=*), optional,     intent(in) :: file
    integer,          optional,     intent(in) :: unit
    type(Setup)                                :: stp

    integer :: u

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='old', action='read', &
            form='unformatted')
    else
       write(0,*) "Error: neither unit number nor file name given"
       write(0,*) "       in `load_Setup'."
       return
    end if

    stp%ntypes_global = size(global_types(:))

    read(u) stp%description
    read(u) stp%atomtype
    read(u) stp%nenv
    allocate(stp%envtypes(stp%nenv))
    read(u) stp%envtypes
    read(u) stp%Rc_min
    read(u) stp%Rc_max
    read(u) stp%sftype
    read(u) stp%nsf
    read(u) stp%nsfparam
    allocate(stp%sf(stp%nsf),                    &
             stp%sfparam(stp%nsfparam, stp%nsf), &
             stp%sfenv(NENV_MAX,stp%nsf),        &
             stp%sfval_min(stp%nsf),             &
             stp%sfval_max(stp%nsf),             &
             stp%sfval_avg(stp%nsf),             &
             stp%sfval_cov(stp%nsf)              )
    read(u) stp%sf(:)
    read(u) stp%sfparam(:,:)
    read(u) stp%sfenv(:,:)
    read(u) stp%neval
    read(u) stp%sfval_min(:)
    read(u) stp%sfval_max(:)
    read(u) stp%sfval_avg(:)
    read(u) stp%sfval_cov(:)

    if (.not. present(unit)) then
       close(u)
    end if

    ! connect local atom type IDs with global ones
    allocate(stp%gtype(stp%nenv), stp%ltype(stp%ntypes_global))
    call stp_set_global_types(stp, stp%ntypes_global, global_types)

    stp%init = .true.

  end function load_Setup

  !--------------------------------------------------------------------!

  function load_Setup_ASCII(global_types, file, unit) result(stp)

    implicit none

    character(len=*), dimension(:), intent(in) :: global_types
    character(len=*), optional,     intent(in) :: file
    integer,          optional,     intent(in) :: unit
    type(Setup)                                :: stp

    integer :: i, j
    integer :: u

    character(len=*), parameter :: dfrmt = '(4(1x,ES24.17))'
    character(len=*), parameter :: ifrmt = '(4(1x,I17))'

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='old', action='read')
    else
       write(0,*) "Error: neither unit number nor file name given"
       write(0,*) "       in `load_Setup'."
       return
    end if

    stp%ntypes_global = size(global_types(:))

    read(u,'(A)') stp%description
    read(u,'(A)') stp%atomtype
    read(u,*) stp%nenv
    allocate(stp%envtypes(stp%nenv))
    read(u,'(A)') (stp%envtypes(i), i=1,stp%nenv)
    read(u,*) stp%Rc_min
    read(u,*) stp%Rc_max
    read(u,'(A)') stp%sftype
    read(u,*) stp%nsf
    read(u,*) stp%nsfparam
    allocate(stp%sf(stp%nsf),                    &
             stp%sfparam(stp%nsfparam, stp%nsf), &
             stp%sfenv(NENV_MAX,stp%nsf),        &
             stp%sfval_min(stp%nsf),             &
             stp%sfval_max(stp%nsf),             &
             stp%sfval_avg(stp%nsf),             &
             stp%sfval_cov(stp%nsf)              )
    read(u,ifrmt) (stp%sf(i), i=1,stp%nsf)
    read(u,dfrmt) ((stp%sfparam(i,j), i=1,stp%nsfparam), j=1,stp%nsf)
    read(u,ifrmt) ((stp%sfenv(i,j), i=1,stp%nenv), j=1,stp%nsf)
    read(u,*) stp%neval
    read(u,dfrmt) (stp%sfval_min(i), i=1,stp%nsf)
    read(u,dfrmt) (stp%sfval_max(i), i=1,stp%nsf)
    read(u,dfrmt) (stp%sfval_avg(i), i=1,stp%nsf)
    read(u,dfrmt) (stp%sfval_cov(i), i=1,stp%nsf)

    if (.not. present(unit)) then
       close(u)
    end if

    ! connect local atom type IDs with global ones
    allocate(stp%gtype(stp%nenv), stp%ltype(stp%ntypes_global))
    call stp_set_global_types(stp, stp%ntypes_global, global_types)

    stp%init = .true.

  end function load_Setup_ASCII

  !--------------------------------------------------------------------!

  subroutine skip_Setup(u)

    implicit none

    integer, intent(in) :: u

    read(u) ! stp%description
    read(u) ! stp%atomtype
    read(u) ! stp%nenv
    read(u) ! stp%envtypes
    read(u) ! stp%Rc_min
    read(u) ! stp%Rc_max
    read(u) ! stp%sftype
    read(u) ! stp%nsf
    read(u) ! stp%nsfparam
    read(u) ! stp%sf(:)
    read(u) ! stp%sfparam(:,:)
    read(u) ! stp%sfenv(:,:)
    read(u) ! stp%neval
    read(u) ! stp%sfval_min(:)
    read(u) ! stp%sfval_max(:)
    read(u) ! stp%sfval_avg(:)
    read(u) ! stp%sfval_cov(:)

  end subroutine skip_Setup

  !--------------------------------------------------------------------!

  subroutine stp_set_global_types(stp, ntypes_global, global_types)

    implicit none

    type(Setup),                                intent(inout) :: stp
    integer,                                    intent(in)    :: ntypes_global
    character(len=*), dimension(ntypes_global), intent(in)    :: global_types

    integer :: i, j

    ! atom type indices that have no corresponding entry are set to 0

    do i = 1, ntypes_global
       stp%ltype(i) = 0
       env : do j = 1, stp%nenv
          if (trim(global_types(i)) == trim(stp%envtypes(j))) then
             stp%ltype(i) = j
             exit env
          end if
       end do env
    end do

    ! reverse direction, because we do not know, if the sets of types
    ! are identical
    do i = 1, stp%nenv
       stp%gtype(i) = 0
       global : do j = 1, ntypes_global
          if (trim(global_types(j)) == trim(stp%envtypes(i))) then
             stp%gtype(i) = j
             exit global
          end if
       end do global
    end do

  end subroutine stp_set_global_types



  !====================================================================!
  !                                                                    !
  !       structural fingerprint basis function set-up module          !
  !                                                                    !
  !====================================================================!



  !--------------------------------------------------------------------!
  !                    initialization/finalization                     !
  !--------------------------------------------------------------------!

  subroutine stp_init(ntypes, stp, N_nb_max)

    implicit none

    integer,                        intent(in)  :: ntypes
    type(Setup), dimension(ntypes), intent(in)  :: stp
    integer,                        intent(in)  :: N_nb_max

    integer            :: itype
    character(len=100) :: sftype
    integer            :: nG_max

    if (isInit) then
       write(0,*) "Error: module already initialized in `stp_init'."
       stop
    end if

    sftype = trim(stp(1)%sftype)
    do itype = 1, ntypes
       if (trim(stp(itype)%sftype) /= trim(sftype)) then
          write(0,*) "Error: Mixing of basis functions of different " &
                         // "types not yet implemented."
          write(0,*) trim(sftype), trim(stp(itype)%sftype)
          stop
       end if
    end do

    select case(trim(io_lower(sftype)))
    case('chebyshev')
       allocate(sfb(ntypes))
       do itype = 1, ntypes
          call setup_basis_chebyshev(stp(itype), sfb(itype))
       end do
    case('behler2011')
       nG_max = 0
       do itype = 1, ntypes
          nG_max = max(nG_max, stp(itype)%nsf)
       end do
       call sf_init(ntypes, nG_max)
       do itype = 1, ntypes
          call setup_symmfunc_Behler2011(itype, stp(itype))
       end do
    case default
       write(0,*) "Error: Unknown basis function type : ", trim(sftype)
       stop
    end select

    ! store max number of SFs in module
    nsf_max = stp_nsf_max(stp)

    ! max number of atoms in interaction range
    nnb_max = N_nb_max

    isInit = .true.

  end subroutine stp_init

  !--------------------------------------------------------------------!

  subroutine stp_final(nTypes, stp)

    implicit none

    integer,                        intent(in) :: nTypes
    type(Setup), dimension(nTypes), intent(in) :: stp

    integer :: itype

    if (.not. isInit) return

    do itype = 1, nTypes
       select case(trim(io_lower(stp(itype)%sftype)))
       case('chebyshev')
          if (allocated(sfb)) deallocate(sfb)
       case('behler2011')
          ! multiple calls to sf_final() do not cause harm
          call sf_final()
       end select
    end do

    isInit = .false.

  end subroutine stp_final

  !--------------------------------------------------------------------!
  !              determine interaction range from setups               !
  !--------------------------------------------------------------------!

  subroutine stp_get_range(nTypes, stp, Rc_min, Rc_max)

    implicit none

    integer,                        intent(in)  :: nTypes
    type(Setup), dimension(nTypes), intent(in)  :: stp
    double precision,               intent(out) :: Rc_min
    double precision,               intent(out) :: Rc_max

    integer :: itype

    call stp_assert_init(stp(1))

    Rc_min = stp(1)%Rc_min
    Rc_max = stp(1)%Rc_max

    do itype = 1, nTypes
       call stp_assert_init(stp(itype))
       Rc_min = min(Rc_min, stp(itype)%Rc_min)
       Rc_max = max(Rc_max, stp(itype)%Rc_max)
    end do

  end subroutine stp_get_range

  !--------------------------------------------------------------------!
  !                      basis function evaluation                     !
  !--------------------------------------------------------------------!

  subroutine stp_eval(itype0, coo0, n, coo1, type1, stp, sfval, &
                      sfderiv_i, sfderiv_j, scaled)

    implicit none

    integer,                                      intent(in)    :: itype0
    double precision, dimension(3),               intent(in)    :: coo0
    integer,                                      intent(in)    :: n
    double precision, dimension(3,n),             intent(in)    :: coo1
    integer,          dimension(n),               intent(in)    :: type1
    type(Setup),                                  intent(inout) :: stp
    double precision, dimension(:),               intent(out)   :: sfval
    double precision, dimension(:,:),   optional, intent(out)   :: sfderiv_i
    double precision, dimension(:,:,:), optional, intent(out)   :: sfderiv_j
    logical,                            optional, intent(in)    :: scaled

    integer               :: type0_loc
    integer, dimension(n) :: type1_loc

    integer          :: i, nsf
    logical          :: do_deriv, do_scale

    call stp_assert_moduleinit()
    call stp_assert_init(stp)

    if (present(sfderiv_i) .and. present(sfderiv_j)) then
       do_deriv = .true.
    else
       do_deriv = .false.
    end if

    if (present(scaled)) then
       do_scale = scaled
    else
       do_scale = .false.
    end if

    ! convert global atom type IDs to setup local IDs
    type0_loc = stp%ltype(itype0)
    do i = 1, n
       type1_loc(i) = stp%ltype(type1(i))
    end do

    select case(trim(io_lower(stp%sftype)))
    case('chebyshev')
       nsf = stp%nsf
       if (do_deriv) then
          call sfb_eval(sfb(itype0), type0_loc, coo0, n, type1_loc, coo1, &
                        nsf, sfval(1:nsf), sfderiv_i(1:3,1:nsf), &
                        sfderiv_j(1:3,1:nsf,1:n))
       else
          call sfb_eval(sfb(itype0), type0_loc, coo0, n, type1_loc, coo1, &
                        nsf, sfval(1:nsf))
       end if
    case('behler2011')
       nsf = stp%nsf
       if (do_deriv) then
          call sf_fingerprint(type0_loc, coo0, n, coo1, type1_loc, nsf, &
                              sfval(1:nsf), sfderiv_i(1:3,1:nsf), &
                              sfderiv_j(1:3,1:nsf,1:n))
       else
          call sf_fingerprint(type0_loc, coo0, n, coo1, type1_loc, nsf, &
                              sfval(1:nsf))
       end if
    end select

    if (do_scale) then
       if (do_deriv) then
          call stp_normalize(stp, sfval, sfderiv_i, sfderiv_j, n)
       else
          call stp_normalize(stp, sfval)
       end if
    else
       if (stp%neval == 0) then
          stp%sfval_min(1:nsf) = sfval(1:nsf)
          stp%sfval_max(1:nsf) = sfval(1:nsf)
          stp%sfval_avg(1:nsf) = sfval(1:nsf)
          stp%sfval_cov(1:nsf) = sfval(1:nsf)*sfval(1:nsf)
       else
          do i = 1, stp%nsf
             stp%sfval_min(i) = min(stp%sfval_min(i), sfval(i))
             stp%sfval_max(i) = max(stp%sfval_max(i), sfval(i))
             stp%sfval_avg(i) = (dble(stp%neval)*stp%sfval_avg(i) &
                                + sfval(i))/(dble(stp%neval+1))
             stp%sfval_cov(i) = (dble(stp%neval)*stp%sfval_cov(i) &
                                + sfval(i)*sfval(i))/(dble(stp%neval+1))
          end do
       end if
    end if

    stp%neval = stp%neval + 1

  end subroutine stp_eval

  !--------------------------------------------------------------------!
  !          normalization of basis function values to [-1,1]          !
  !--------------------------------------------------------------------!

  subroutine stp_normalize(stp, sfval, sfderiv_i, sfderiv_j, n)

    implicit none

    type(Setup),                                  intent(inout) :: stp
    double precision, dimension(:),               intent(inout) :: sfval
    double precision, dimension(:,:),   optional, intent(inout) :: sfderiv_i
    double precision, dimension(:,:,:), optional, intent(inout) :: sfderiv_j
    integer,                            optional, intent(in)    :: n

    double precision :: scale, shift, s
    integer          :: isf
    logical          :: do_deriv

    if (present(sfderiv_i) .and. present(sfderiv_j) .and. present(n)) then
       do_deriv = .true.
    else
       do_deriv = .false.
    end if

    do isf = 1, stp%nsf
       shift = stp%sfval_avg(isf)
       ! scale covariance to 1
       ! s = sqrt(stp%sfval_cov(isf) + shift*shift - 2.0d0*shift*stp%sfval_avg(isf))
       s = sqrt(stp%sfval_cov(isf) - shift*shift)
       if (s <= 1.0d-10) then
          write(0,*) "Warning: Invalid scaling encountered in ", &
                               "'stp_normalize()'."
          write(0,*) "         This means at least one fingerprint ", &
                               "function for ", trim(adjustl(stp%atomtype)), &
                               " is always equal to zero!"
          write(0,*) "         Maybe an atomic species is not present ", &
                               "in the reference set?"
          write(0,*) "         type       = ", trim(adjustl(stp%sftype)), &
                               " ", trim(io_adjustl(stp%sf(isf)))
          write(0,*) "         covariance = ", stp%sfval_cov(isf)
          write(0,*) "         average    = ", stp%sfval_avg(isf)
          write(0,*) "         min, max   = ", stp%sfval_min(isf), &
                                               stp%sfval_max(isf)
          scale = 0.0d0
       else
          scale = 1.0d0/s
       end if
       sfval(isf) = scale*(sfval(isf)-shift)
       if (do_deriv) then
          sfderiv_i(1:3,isf) = scale*sfderiv_i(1:3,isf)
          sfderiv_j(1:3,isf,1:n) = scale*sfderiv_j(1:3,isf,1:n)
       end if
    end do

  end subroutine stp_normalize

  !--------------------------------------------------------------------!
  !            maximum number of basis functions per setup             !
  !--------------------------------------------------------------------!

  function stp_nsf_max(stp) result(N_sf_max)

    implicit none

    type(Setup), dimension(:), optional, intent(in) :: stp
    integer                                         :: N_sf_max

    integer :: itype, nTypes

    if (present(stp)) then
       nTypes = size(stp(:))
       N_sf_max = 0
       do itype = 1, nTypes
          call stp_assert_init(stp(itype))
          N_sf_max = max(N_sf_max, stp(itype)%nsf)
       end do
    else if (isInit) then
       N_sf_max = nsf_max
    else
       write(0,*) "Error: module not initialized in `stp_nsf_max()'."
       stop
    end if

  end function stp_nsf_max

  !--------------------------------------------------------------------!
  !                        auxiliary procedures                        !
  !--------------------------------------------------------------------!

  subroutine stp_assert_init(stp)
    implicit none
    type(Setup), intent(in) :: stp
    if (.not. stp%init) then
       write(*,*) 'Error: The basis function set-up is not initialized.'
       write(*,*)
       stop
    end if
  end subroutine stp_assert_init

  !--------------------------------------------------------------------!

  subroutine stp_assert_moduleinit()
    implicit none
    if (.not. isInit) then
       write(*,*) 'Error: Structural fingerprint setup module NOT initialized.'
       write(*,*)
       stop
    end if
  end subroutine stp_assert_moduleinit





  !====================================================================!
  !                                                                    !
  !            procedures for different basis function types           !
  !                                                                    !
  !====================================================================!


  subroutine read_basis_chebyshev(u_stp, stp, iline)

    implicit none

    integer,     intent(in)    :: u_stp
    type(Setup), intent(inout) :: stp
    integer,     intent(inout) :: iline

    character(len=1024) :: line
    integer             :: r_N, a_N
    double precision    :: r_Rc, a_Rc

    read(u_stp, '(A)') line
    iline = iline + 1

    r_Rc = 0.0d0
    a_Rc = 0.0d0
    r_N = 0
    a_N = 0

    call io_readval(line, 'radial_Rc', r_Rc)
    call io_readval(line, 'radial_N', r_N)
    call io_readval(line, 'angular_Rc', a_Rc)
    call io_readval(line, 'angular_N', a_N)

    stp%nsf = r_N + a_N + 2
    if (stp%nenv > 1) then
       stp%nsf = 2*stp%nsf
    end if

    ! FIXME: most of these allocations are not required for the
    ! Chebyshev basis, but we need to allocate the memory to stay
    ! compatible with the current 'save' and 'load' routines
    allocate(stp%sf(stp%nsf), stp%sfparam(NSFPARAM,stp%nsf), &
         stp%sfval_min(stp%nsf), stp%sfval_max(stp%nsf), &
         stp%sfval_avg(stp%nsf), stp%sfval_cov(stp%nsf), &
         stp%sfenv(NENV_MAX,stp%nsf))

    stp%nsfparam     = NSFPARAM
    stp%Rc_max       = 0.0d0
    stp%sf(:)        = 0
    stp%sfparam(:,:) = 0.0d0
    stp%sfenv(:,:)   = 0
    stp%sfval_min(:) = 0.0d0
    stp%sfval_max(:) = 0.0d0
    stp%sfval_avg(:) = 0.0d0
    stp%sfval_cov(:) = 0.0d0

    ! use only the first row of sfparam to store the actual parameters
    stp%sfparam(1,1) = r_Rc
    stp%sfparam(2,1) = dble(r_N)
    stp%sfparam(3,1) = a_Rc
    stp%sfparam(4,1) = dble(a_N)

    stp%Rc_max = max(r_Rc, a_Rc)

  end subroutine read_basis_chebyshev

  !--------------------------------------------------------------------!

  subroutine setup_basis_chebyshev(stp, sfb)

    implicit none

    type(Setup),            intent(in)  :: stp
    type(FingerprintBasis), intent(out) :: sfb

    double precision :: r_Rc, a_Rc
    integer          :: r_N, a_N

    r_Rc = stp%sfparam(1,1)
    r_N = nint(stp%sfparam(2,1))
    a_Rc = stp%sfparam(3,1)
    a_N = nint(stp%sfparam(4,1))

    sfb = new_SFBasis(stp%nenv, stp%envtypes, r_N, a_N, r_Rc, a_Rc)

  end subroutine setup_basis_chebyshev

  !--------------------------------------------------------------------!

  subroutine print_info_chebyshev(stp)

    implicit none

    type(Setup), intent(in) :: stp

    double precision :: r_Rc, a_Rc
    integer          :: r_N, a_N

    r_Rc = stp%sfparam(1,1)
    r_N = nint(stp%sfparam(2,1))
    a_Rc = stp%sfparam(3,1)
    a_N = nint(stp%sfparam(4,1))

    write(*,*) 'Basis function type Chebyshev'
    write(*,*) '[N. Artrith and A. Urban (2016)]'
    write(*,*)
    write(*,*) 'Radial Rc     : ' // trim(io_adjustl(r_Rc))
    write(*,*) 'Angular Rc    : ' // trim(io_adjustl(a_Rc))
    write(*,*) 'Radial order  : ' // trim(io_adjustl(r_N))
    write(*,*) 'Angular order : ' // trim(io_adjustl(a_N))
    write(*,*)

  end subroutine print_info_chebyshev



  !--------------------------------------------------------------------!
  ! J. Behler, J. Chem. Phys. 134 (2011) 074106                        !
  !                                                                    !
  !      parameters                                                    !
  !                                                                    !
  ! G1:  Rc                                                            !
  ! G2:  Rc, Rs, eta                                                   !
  ! G3:  Rc, kappa                                                     !
  ! G4:  Rc, lambda, zeta, eta                                         !
  ! G5:  Rc, lambda, zeta, eta                                         !
  !--------------------------------------------------------------------!

  subroutine read_symmfunc_Behler2011(u_stp, isf, stp, iline)

    implicit none

    integer,     intent(in)    :: u_stp
    integer,     intent(in)    :: isf
    type(Setup), intent(inout) :: stp
    integer,     intent(inout) :: iline

    character(len=1024)        :: line, typename

    read(u_stp, '(A)') line
    iline = iline + 1

    call io_readval(line, 'G', stp%sf(isf))
    select case(stp%sf(isf))
    case(1)
       call io_readval(line, 'type2',  typename)
       stp%sfenv(1,isf) = stp_typeindex(stp, typename)
       call io_readval(line, 'Rc',     stp%sfparam(1,isf))
    case(2)
       call io_readval(line, 'type2',  typename)
       stp%sfenv(1,isf) = stp_typeindex(stp, typename)
       call io_readval(line, 'Rc',     stp%sfparam(1,isf))
       call io_readval(line, 'Rs',     stp%sfparam(2,isf))
       call io_readval(line, 'eta',    stp%sfparam(3,isf))
    case(3)
       call io_readval(line, 'type2',  typename)
       stp%sfenv(1,isf) = stp_typeindex(stp, typename)
       call io_readval(line, 'Rc',     stp%sfparam(1,isf))
       call io_readval(line, 'kappa',  stp%sfparam(2,isf))
    case(4)
       call io_readval(line, 'type2',  typename)
       stp%sfenv(1,isf) = stp_typeindex(stp, typename)
       call io_readval(line, 'type3',  typename)
       stp%sfenv(2,isf) = stp_typeindex(stp, typename)
       call io_readval(line, 'Rc',     stp%sfparam(1,isf))
       call io_readval(line, 'lambda', stp%sfparam(2,isf))
       call io_readval(line, 'zeta',   stp%sfparam(3,isf))
       call io_readval(line, 'eta',    stp%sfparam(4,isf))
    case(5)
       call io_readval(line, 'type2',  typename)
       stp%sfenv(1,isf) = stp_typeindex(stp, typename)
       call io_readval(line, 'type3',  typename)
       stp%sfenv(2,isf) = stp_typeindex(stp, typename)
       call io_readval(line, 'Rc',     stp%sfparam(1,isf))
       call io_readval(line, 'lambda', stp%sfparam(2,isf))
       call io_readval(line, 'zeta',   stp%sfparam(3,isf))
       call io_readval(line, 'eta',    stp%sfparam(4,isf))
    case default
       write(0,*) "Error: Symmetry function type not implemented : ", &
                  trim(io_adjustl(stp%sf(isf)))
       write(0,*) "in   : `read_symmfunc_Behler2011'"
       stop
    end select

    stp%Rc_max = max(stp%Rc_max, stp%sfparam(1,isf))

  end subroutine read_symmfunc_Behler2011

  !--------------------------------------------------------------------!
  ! set up symmetry functions of the Behler2011 type as specified by   !
  ! the provided set-up `stp'                                          !
  !--------------------------------------------------------------------!

  subroutine setup_symmfunc_Behler2011(itype, stp)

    implicit none

    integer,     intent(in) :: itype
    type(Setup), intent(in) :: stp

    integer          :: type1, type2, type3
    double precision :: Rc, Rs, eta, kappa, lambda, zeta
    integer          :: isf, funct

    type1 = itype

    SFs : do isf = 1, stp%nsf

       funct = stp%sf(isf)

       select case(funct)
       case(1)
          type2  = stp%sfenv(1,isf)
          Rc     = stp%sfparam(1,isf)
          call sf_add_rad(funct, type1, type2, Rc=Rc)
       case(2)
          type2  = stp%sfenv(1,isf)
          Rc     = stp%sfparam(1,isf)
          Rs     = stp%sfparam(2,isf)
          eta    = stp%sfparam(3,isf)
          call sf_add_rad(funct, type1, type2, Rc=Rc, Rs=Rs, eta=eta)
       case(3)
          type2  = stp%sfenv(1,isf)
          Rc     = stp%sfparam(1,isf)
          kappa  = stp%sfparam(2,isf)
          call sf_add_rad(funct, type1, type2, Rc=Rc, kappa=kappa)
       case(4)
          type2  = stp%sfenv(1,isf)
          type3  = stp%sfenv(2,isf)
          Rc     = stp%sfparam(1,isf)
          lambda = stp%sfparam(2,isf)
          zeta   = stp%sfparam(3,isf)
          eta    = stp%sfparam(4,isf)
          call sf_add_ang(funct, type1, type2, type3, Rc=Rc, lambda=lambda, &
                      zeta=zeta, eta=eta)
       case(5)
          type2  = stp%sfenv(1,isf)
          type3  = stp%sfenv(2,isf)
          Rc     = stp%sfparam(1,isf)
          lambda = stp%sfparam(2,isf)
          zeta   = stp%sfparam(3,isf)
          eta    = stp%sfparam(4,isf)
          call sf_add_ang(funct, type1, type2, type3, Rc=Rc, lambda=lambda, &
                      zeta=zeta, eta=eta)
       case default
          write(0,*) "Error: Symmetry function type not implemented : ", &
                     trim(io_adjustl(stp%sf(isf)))
          write(0,*) "in   : `setup_symmfunc_Behler2011'"
          stop
       end select

   end do SFs

  end subroutine setup_symmfunc_Behler2011

  !--------------------------------------------------------------------!
  !                        print info to stdout                        !
  !--------------------------------------------------------------------!

  subroutine print_info_Behler2011(stp)

    implicit none

    type(Setup), intent(in) :: stp

    integer            :: isf, iG
    double precision   :: Rc, Rs, eta, kappa, lambda, zeta
    character(len=128) :: frmt
    integer            :: type2, type3

    write(*,*) 'Basis function type Behler2011'
    write(*,*)
    write(*,*) '[see also: J. Behler, J. Chem. Phys. 134 (2011) 074106]'
    write(*,*)
    write(*,*)

    write(*,'(6x,"G",2x,"parameters")')
    write(*,*)

    do isf = 1, stp%nsf
       iG = stp%sf(isf)
       select case(iG)
       case(1)
          type2 = stp%sfenv(1,isf)
          Rc = stp%sfparam(1,isf)
          write(*,'(1x,I3,2x,I1,2x,"type2=",A2,"  Rc = ",F7.3)') &
               isf, iG, stp%envtypes(type2), Rc
       case(2)
          type2 = stp%sfenv(1,isf)
          Rc  = stp%sfparam(1,isf)
          Rs  = stp%sfparam(2,isf)
          eta = stp%sfparam(3,isf)
          write(*,'(1x,I3,2x,I1,2x,"type2=",A2,"  Rc = ",F7.3,"  Rs = ",F7.3,"  eta = ",F7.3)') &
               isf, iG, stp%envtypes(type2), Rc, Rs, eta
       case(3)
          type2 = stp%sfenv(1,isf)
          Rc    = stp%sfparam(1,isf)
          kappa = stp%sfparam(2,isf)
          write(*,'(1x,I3,2x,I1,2x,"type2=",A2,"  Rc = ",F7.3,"  kappa = ",F7.3)') &
               isf, iG, stp%envtypes(type2), Rc, kappa
       case(4,5)
          type2  = stp%sfenv(1,isf)
          type3  = stp%sfenv(2,isf)
          Rc     = stp%sfparam(1,isf)
          lambda = stp%sfparam(2,isf)
          zeta   = stp%sfparam(3,isf)
          eta    = stp%sfparam(4,isf)
          frmt = '(1x,I3,2x,I1,2x,"type2=",A2,"  type3=",A2,' &
               // '"  Rc = ",F7.3,"  lambda = ",F7.3,' &
               // '"  zeta = ",F7.3,"  eta = ",F7.3)'
          write(*,frmt) isf, iG, stp%envtypes(type2), &
               stp%envtypes(type3), Rc, lambda, zeta, eta
       end select
    end do
    write(*,*)

  end subroutine print_info_Behler2011

  !--------------------------------------------------------------------!

  function stp_typeindex(stp, typename) result(itype)

    implicit none

    type(Setup),      intent(in) :: stp
    character(len=*), intent(in) :: typename
    integer                      :: itype

    integer :: i

    itype = 0
    do i = 1, stp%nenv
       if (trim(typename) == trim(stp%envtypes(i))) then
          itype = i
          exit
       end if
    end do

  end function stp_typeindex

end module sfsetup
