!-----------------------------------------------------------------------
!           trainset.f90 - handling of the training set file
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
! 2013-05-09 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

module trainset

  use aeio,    only: aeio_header,           &
                     TYPELEN, PATHLEN

  use io,      only: io_adjustl,            &
                     io_unit

  use sfsetup, only: Setup,                 &
                     save_Setup,            &
                     load_Setup,            &
                     del_Setup,             &
                     stp_init,              &
                     stp_final,             &
                     stp_normalize,         &
                     stp_assert_moduleinit, &
                     stp_nsf_max
  implicit none
  private
  save

  public  :: new_TrnSet,              &
             open_TrnSet,             &
             close_TrnSet,            &
             rewind_TrnSet,           &
             new_TrnSet_info,         &
             save_TrnSet_info,        &
             save_TrnSet_info_ASCII,  &
             load_TrnSet_info,        &
             load_TrnSet_info_ASCII,  &
             save_TrnSet_ASCII,       &
             load_TrnSet_ASCII,       &
             ts_count_atoms,          &
             ts_load_Setups,          &
             ts_normalize,            &
             ts_print_info,           &
             ts_read_atom_info,       &
             ts_read_header,          &
             ts_read_sf_info,         &
             ts_read_sf_values,       &
             ts_read_structure_info,  &
             ts_skip_header,          &
             ts_skip_structure,       &
             ts_skip_atoms,           &
             ts_unload_Setups,        &
             ts_write_atom_info,      &
             ts_write_header,         &
             ts_write_sf_info,        &
             ts_write_structure_info, &
             ts_read_footer,          &
             ts_write_footer


  private :: ts_assert_init,          &
             ts_assert_notinit,       &
             ts_assert_writemode,     &
             ts_assert_readmode,      &
             ts_energy_norm_simple

  type, public :: TrnSet

     !-----------------------------------------------------------------!
     ! init        .true., if the training set has been initialized    !
     ! file        name of the corresponding training set file         !
     ! unit        unit number of that file                            !
     ! mode        current access mode; 'read', 'write', 'info'        !
     ! normalized  .true., if the input and output values have been    !
     !             normalized to the interval [-1,1] ('read' mode only)!
     !                                                                 !
     ! if (normalized == .true.)                                       !
     ! scale       energy scaling factor used for the normalization    !
     ! shift       atomic energy shift used for energy normalization   !
     !                                                                 !
     ! nTypes      number of atomic species in the training set        !
     ! nAtomsTot   total number of atoms in the training set           !
     ! typeName(i) name of i-th atomic species                         !
     ! nStrucs     total number of structures in the training set      !
     !             --> when open in 'write' mode, not necessarily all  !
     !                 files have yet been parsed                      !
     ! iStruc      current file record position (0=before first file)  !
     !-----------------------------------------------------------------!

     logical                                           :: init = .false.

     character(len=PATHLEN)                            :: file
     integer                                           :: unit
     character(len=5)                                  :: mode

     logical                                           :: normalized
     double precision                                  :: scale, shift

     integer                                           :: nTypes
     integer                                           :: nAtomsTot
     character(len=TYPELEN), dimension(:), allocatable :: typeName
     double precision,       dimension(:), allocatable :: E_atom
     integer                                           :: nStrucs
     integer                                           :: iStruc

     double precision                                  :: E_min, E_max, E_av

  end type TrnSet

  !--------------------------------------------------------------------!
  ! Basis function values and derivatives may be read and written      !
  ! either using a basis function setup [type(Setup)] or directly into !
  ! double precision arrays of the correct dimensions.                 !
  !--------------------------------------------------------------------!

contains

  function new_TrnSet(nTypes, typeName, E_atom, nStrucs, file, scale, &
                      shift) result(ts)

    implicit none

    integer,                             intent(in) :: nTypes
    character(len=*), dimension(nTypes), intent(in) :: typeName
    double precision, dimension(nTypes), intent(in) :: E_atom
    integer,                             intent(in) :: nStrucs
    character(len=*),                    intent(in) :: file
    double precision, optional,          intent(in) :: scale
    double precision, optional,          intent(in) :: shift
    type(TrnSet)                                    :: ts

    logical :: fexists

    inquire(file=trim(adjustl(file)), exist=fexists)
    if (fexists) then
       write(0,*) 'Error: file already exists: ', trim(adjustl(file))
       stop
    end if

    allocate(ts%typeName(nTypes), ts%E_atom(nTypes))

    ts%file               = trim(adjustl(file))
    ts%nTypes             = nTypes
    ts%typeName(1:nTypes) = typeName(1:nTypes)
    ts%E_atom(1:nTypes)   = E_atom(1:nTypes)
    ts%nStrucs            = nStrucs
    ts%iStruc             = 0
    if (present(scale) .and. present(shift)) then
       ts%normalized = .true.
       ts%scale = scale
       ts%shift = shift
    else
       ts%normalized = .false.
       ts%scale = 1.0d0
       ts%shift = 0.0d0
    end if

    ts%unit   = io_unit()
    open(ts%unit, file=trim(ts%file), status='new', action='write', &
         form='unformatted')

    ts%nAtomsTot = 0

    ts%mode = 'write'
    ts%init = .true.

    call ts_write_header(ts)

  end function new_TrnSet

  !--------------------------------------------------------------------!

  function open_TrnSet(file, maxenergy, readfooter, raw) result(ts)

    implicit none

    !------------------------------------------------------------------!
    ! file         path to the training set file                       !
    ! maxenergy    maximum atomic energy (do not consider structures   !
    !              with higher energy)                                 !
    ! readfooter   skip over all structures to read footer first       !
    ! raw          if true, do not normalize the trainset              !
    !------------------------------------------------------------------!

    character(len=*),           intent(in) :: file
    double precision, optional, intent(in) :: maxenergy
    logical,          optional, intent(in) :: readfooter
    logical,          optional, intent(in) :: raw
    type(TrnSet)                           :: ts

    logical :: do_footer, do_raw

    if (present(raw)) then
       do_raw = raw
    else
       do_raw = .false.
    end if    
    
    call ts_read_header(ts, file, do_raw)

    ! To read the footer we have to skip over all structures.
    ! The flag `readfooter=.false.' can thus save some time, if
    ! the stats from the footer are not needed.
    if (present(readfooter)) then
       do_footer = readfooter
    else
       do_footer = .true.
    end if

    ts%iStruc = 0
    ts%mode   = 'read'
    ts%init   = .true.

    if (.not. do_raw) then
       if (.not. ts%normalized) then
          if (present(maxenergy)) then
             call ts_normalize(ts, maxenergy)
          else
             call ts_normalize(ts, huge(1.0d0))
          end if
       end if
    end if

    if (do_footer) then
       call ts_read_footer(ts)
       call rewind_TrnSet(ts)
    else if (ts%normalized) then
       ts%E_av  =  0.0d0
       ts%E_min = -1.0d0
       ts%E_max =  1.0d0
    end if

  end function open_TrnSet

  !--------------------------------------------------------------------!

  subroutine close_TrnSet(ts, stp, status)

    implicit none

    type(TrnSet),                        intent(inout) :: ts
    type(Setup), dimension(:), optional, intent(in)    :: stp
    character(len=*),          optional, intent(in)    :: status

    if (.not. ts%init) return

    if (trim(ts%mode) == 'write') then
       if (present(stp)) then
          call ts_write_footer(ts, stp=stp)
       else
          call ts_write_footer(ts)
       end if
    end if

    if ((trim(ts%mode)=='read') .or. (trim(ts%mode)=='write')) then
       if (present(status)) then
          close(ts%unit, status=trim(status))
       else
          close(ts%unit)
       end if
    end if

    deallocate(ts%typeName, ts%E_atom)
    ts%init   = .false.

  end subroutine close_TrnSet

  !--------------------------------------------------------------------!

  subroutine rewind_TrnSet(ts, rec)

    implicit none

    type(TrnSet),      intent(inout) :: ts
    integer, optional, intent(in)    :: rec

    integer :: irec

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    if (present(rec) .and. (rec > 0)) then
       irec = rec - 1
    else
       irec = 0
    end if

    if (ts%iStruc > irec) then
       rewind(ts%unit)
       ! skip header
       call ts_skip_header(ts)
       ! reset structure record pointer
       ts%iStruc = 0
    end if

    do while(ts%iStruc < irec)
       ! fast-forward to desired structure
       call ts_skip_structure(ts)
    end do

  end subroutine rewind_TrnSet

  !--------------------------------------------------------------------!
  !              only training set info - no actual data               !
  !--------------------------------------------------------------------!

  subroutine save_TrnSet_info(ts, file, unit)

    implicit none

    type(TrnSet),               intent(in) :: ts
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer :: u

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='replace', &
            form='unformatted', action='write')
    else
       write(0,*) "Error: neither unit nor file specified in `save_TrnSet_info()'."
       stop
    end if

    write(u) ts%file
    write(u) ts%normalized
    write(u) ts%scale
    write(u) ts%shift
    write(u) ts%nTypes
    write(u) ts%typeName(1:ts%nTypes)
    write(u) ts%E_atom(1:ts%nTypes)
    write(u) ts%nAtomsTot
    write(u) ts%nStrucs
    write(u) ts%E_min, ts%E_max, ts%E_av

    if (.not. present(unit)) close(u)

  end subroutine save_TrnSet_info

  !--------------------------------------------------------------------!

  subroutine save_TrnSet_info_ASCII(ts, file, unit)

    implicit none

    type(TrnSet),               intent(in) :: ts
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer :: i
    integer :: u

    character(len=*), parameter :: dfrmt = '(4(1x,ES24.17))'

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='replace', action='write')
    else
       write(0,*) "Error: neither unit nor file specified in ", &
                  "`save_TrnSet_info()'."
       stop
    end if

    write(u,'(A)') trim(ts%file)
    write(u,*) ts%normalized
    write(u,*) ts%scale
    write(u,*) ts%shift
    write(u,*) ts%nTypes
    write(u,'(A)') (ts%typeName(i), i=1,ts%nTypes)
    write(u,dfrmt) (ts%E_atom(i), i=1,ts%nTypes)
    write(u,*) ts%nAtomsTot
    write(u,*) ts%nStrucs
    write(u,*) ts%E_min, ts%E_max, ts%E_av

    if (.not. present(unit)) close(u)

  end subroutine save_TrnSet_info_ASCII

  !--------------------------------------------------------------------!

  function load_TrnSet_info(file, unit) result(ts)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(TrnSet)                           :: ts

    integer :: u

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='old', &
            form='unformatted', action='read')
    else
       write(0,*) "Error: neither unit nor file specified in ", &
                  "`load_TrnSet_info()'."
       stop
    end if

    read(u) ts%file
    read(u) ts%normalized
    read(u) ts%scale
    read(u) ts%shift
    read(u) ts%nTypes
    allocate(ts%typeName(ts%nTypes), ts%E_atom(ts%nTypes))
    read(u) ts%typeName(1:ts%nTypes)
    read(u) ts%E_atom(1:ts%nTypes)
    read(u) ts%nAtomsTot
    read(u) ts%nStrucs
    read(u) ts%E_min, ts%E_max, ts%E_av

    ts%iStruc = ts%nStrucs
    ts%init = .true.
    ts%mode = 'info'

    if (.not. present(unit)) close(u)

  end function load_TrnSet_info

  !--------------------------------------------------------------------!

  function load_TrnSet_info_ASCII(file, unit) result(ts)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(TrnSet)                           :: ts

    integer :: i
    integer :: u

    character(len=*), parameter :: dfrmt = '(4(1x,ES24.17))'

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='old', action='read')
    else
       write(0,*) "Error: neither unit nor file specified in `load_TrnSet_info()'."
       stop
    end if

    read(u,'(A)') ts%file
    read(u,*) ts%normalized
    read(u,*) ts%scale
    read(u,*) ts%shift
    read(u,*) ts%nTypes
    allocate(ts%typeName(ts%nTypes), ts%E_atom(ts%nTypes))
    read(u,'(A)') (ts%typeName(i), i=1,ts%nTypes)
    read(u,dfrmt) (ts%E_atom(i), i=1,ts%nTypes)
    read(u,*) ts%nAtomsTot
    read(u,*) ts%nStrucs
    read(u,*) ts%E_min, ts%E_max, ts%E_av

    ts%iStruc = ts%nStrucs
    ts%init = .true.
    ts%mode = 'info'

    if (.not. present(unit)) close(u)

  end function load_TrnSet_info_ASCII

  !--------------------------------------------------------------------!

  function new_TrnSet_info(nTypes) result(ts)

    implicit none

    integer,                             intent(in) :: nTypes
    type(TrnSet)                                    :: ts

    allocate(ts%typeName(nTypes), ts%E_atom(nTypes))

    ts%nTypes     = nTypes
    ts%normalized = .false.
    ts%scale      = 1.0d0
    ts%shift      = 1.0d0
    ts%file       = ''
    ts%unit       = -1
    ts%nAtomsTot  = 0
    ts%nStrucs    = 0
    ts%iStruc     = 0

    ts%init = .true.
    ts%mode = 'info'

  end function new_TrnSet_info

  !--------------------------------------------------------------------!
  !               save entire training set to ASCII file               !
  !--------------------------------------------------------------------!

  subroutine save_TrnSet_ASCII(ts, file, unit)

    implicit none

    type(TrnSet),               intent(inout) :: ts
    character(len=*), optional, intent(in)    :: file
    integer,          optional, intent(in)    :: unit

    integer                                     :: istruc
    integer                                     :: iatom, nAtoms
    integer                                     :: itype, nTypes
    character(len=PATHLEN)                      :: filename
    double precision                            :: E_coh
    double precision, dimension(3)              :: forCart, cooCart
    integer                                     :: nsf, nsf_max
    double precision, dimension(:), allocatable :: sfval

    integer :: i
    integer :: u

    character(len=*), parameter :: dfrmt = '(1000(1x,ES24.17))'

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='replace', action='write')
    else
       write(0,*) "Error: neither unit nor file specified in `save_TrnSet()'."
       stop
    end if

    call save_TrnSet_info_ASCII(ts, unit=u)
    call rewind_TrnSet(ts)

    nsf_max = 100
    allocate(sfval(nsf_max))

    do istruc = 1, ts%nStrucs
       call ts_read_structure_info(ts, filename, nAtoms, nTypes, E_coh)
       write(u,'(A)') trim(filename)
       write(u,*) nAtoms, nTypes
       write(u,*) E_coh
       do iatom = 1, nAtoms
          call ts_read_atom_info(ts, itype, cooCart, forCart)
          write(u,*) itype
          write(u,*) cooCart
          write(u,*) forCart
          call ts_read_sf_info(ts, nsf)
          write(u,*) nsf
          if (nsf > nsf_max) then
             nsf_max = nsf
             deallocate(sfval)
             allocate(sfval(nsf_max))
          end if
          call ts_read_sf_values(ts, nsf, sfval(1:nsf))
          write(u,dfrmt,advance='no') (sfval(i), i=1,nsf)
          write(u,*)
       end do       
    end do

    deallocate(sfval)
    if (.not. present(unit)) close(u)

  end subroutine save_TrnSet_ASCII

  !--------------------------------------------------------------------!

  function load_TrnSet_ASCII(trnset_file, file, unit) result(ts)

    implicit none

    !------------------------------------------------------------------!
    ! This will not load the training set into memory, but will        !
    ! transcribe the ASCII format to the binary format.  A new trnset  !
    ! file will be created at path 'trnset_file'.                      !
    !------------------------------------------------------------------!

    character(len=*),           intent(in) :: trnset_file
    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    type(TrnSet)                           :: ts1, ts

    integer                                     :: istruc
    integer                                     :: iatom, nAtoms
    integer                                     :: itype, nTypes
    character(len=PATHLEN)                      :: filename
    double precision                            :: energy
    double precision, dimension(3)              :: forCart, cooCart
    integer                                     :: nsf, nsf_max
    double precision, dimension(:), allocatable :: sfval

    integer :: i
    integer :: u
    logical :: fexists

    character(len=*), parameter :: dfrmt = '(1000(1x,ES24.17))'

    inquire(file=trim(trnset_file), exist=fexists)
    if (fexists) then
       write(0,*) 'Error: Will not overwrite file: ', trim(trnset_file)
       stop
    end if

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='old', action='read')
    else
       write(0,*) "Error: neither unit nor file specified in `load_TrnSet()'."
       stop
    end if

    ts1 = load_TrnSet_info_ASCII(unit=u)
    ts = new_TrnSet(ts1%nTypes, ts1%typeName, ts1%E_atom, ts1%nStrucs, &
                    trnset_file, ts1%scale, ts1%shift)
    call close_TrnSet(ts1)

    nsf_max = 100
    allocate(sfval(nsf_max))

    do istruc = 1, ts%nStrucs
       read(u,'(A)') filename
       read(u,*) nAtoms, nTypes
       read(u,*) energy
       call ts_write_structure_info(ts, filename, nAtoms, nTypes, energy)
       do iatom = 1, nAtoms
          read(u,*) itype
          read(u,*) cooCart
          read(u,*) forCart
          call ts_write_atom_info(ts, itype, cooCart, forCart)
          read(u,*) nsf
          if (nsf > nsf_max) then
             nsf_max = nsf
             deallocate(sfval)
             allocate(sfval(nsf_max))
          end if
          read(u,dfrmt) (sfval(i), i=1,nsf)
          call ts_write_sf_info(ts, nsf, sfval(1:nsf))
       end do
    end do

    deallocate(sfval)
    if (.not. present(unit)) close(u)

  end function load_TrnSet_ASCII


  !--------------------------------------------------------------------!
  !                    info about the training set                     !
  !--------------------------------------------------------------------!

  subroutine ts_print_info(ts)

    implicit none

    type(TrnSet), intent(inout) :: ts

    integer :: itype

    call ts_assert_init(ts)

    call aeio_header("Training set info.")
    write(*,*)

    write(*,*) 'Training set file                   : ', trim(adjustl(ts%file))
    write(*,*) 'Number of structures in the data set: ', trim(io_adjustl(ts%nStrucs))
    if (ts%iStruc /= ts%nStrucs) then
       if (trim(ts%mode) == 'write') then
          write(*,*) '  Structures included so far        : ', trim(io_adjustl(ts%iStruc))
       else
          write(*,*) '  Structures read so far            : ', trim(io_adjustl(ts%iStruc))
       end if
    end if
    write(*,*)

    write(*,*) 'Atomic species in training set      : ', trim(io_adjustl(ts%nTypes))
    write(*,'(1x,"  Species :")', advance='no')
    do itype = 1, ts%nTypes
       write(*,'(1x,A)', advance='no') trim(ts%typeName(itype))
    end do
    write(*,*)
    write(*,*)

    if (ts%normalized .or. (ts%iStruc == ts%nStrucs)) then
       write(*,*) 'Average energy (eV/atom) : ', trim(io_adjustl(ts%E_av,6))
       write(*,*) 'Minimum energy (eV/atom) : ', trim(io_adjustl(ts%E_min,6))
       write(*,*) 'Maximum energy (eV/atom) : ', trim(io_adjustl(ts%E_max,6))
       write(*,*)
    end if

    if (ts%normalized) then
       write(*,*) 'The input and output values have been normalized to [-1.0, 1.0].'
       write(*,*) 'Structures outside of this interval will not be used for training.'
       write(*,*) '  Energy scaling factor: ', trim(io_adjustl(ts%scale,6))
       write(*,*) '  Atomic energy shift  : ', trim(io_adjustl(ts%shift,6))
    else
       write(*,*) 'The input and output values have not yet been normalized.'
    end if
    write(*,*)

  end subroutine ts_print_info

  !--------------------------------------------------------------------!
  !                      training set file header                      !
  !--------------------------------------------------------------------!

  subroutine ts_write_header(ts)

    implicit none

    type(TrnSet), intent(inout) :: ts

    call ts_assert_init(ts)
    call ts_assert_writemode(ts)

    write(ts%unit) ts%nTypes
    write(ts%unit) ts%nStrucs
    write(ts%unit) ts%typeName(:)
    write(ts%unit) ts%E_atom(:)
    write(ts%unit) ts%normalized
    write(ts%unit) ts%scale
    write(ts%unit) ts%shift

  end subroutine ts_write_header

  !--------------------------------------------------------------------!

  subroutine ts_read_header(ts, file, raw)

    implicit none

    type(TrnSet),     intent(inout) :: ts
    character(len=*), intent(in)    :: file
    logical,          intent(in)    :: raw
    
    character(len=1024) :: file2
    logical             :: fexists

    call ts_assert_notinit(ts)

    ! check, if normalized training set file exists:
    file2 = trim(adjustl(file)) // '.scaled'
    inquire(file=trim(file2), exist=fexists)
    if ((.not. raw) .and. fexists) then
       write(*,*) 'Loading scaled training set file: ' // trim(file2)
       write(*,*)
    else
       file2 = trim(adjustl(file))
       inquire(file=trim(file2), exist=fexists)
       if (.not. fexists) then
          write(0,*) 'Error: file not found: ', trim(adjustl(file))
          stop
       end if
    end if

    ts%file = trim(adjustl(file2))
    ts%unit = io_unit()
    open(ts%unit, file=trim(ts%file), status='old', action='read', &
         form='unformatted')

    read(ts%unit) ts%nTypes
    read(ts%unit) ts%nStrucs
    allocate(ts%typeName(ts%nTypes), ts%E_atom(ts%nTypes))
    read(ts%unit) ts%typeName(:)
    read(ts%unit) ts%E_atom(:)
    read(ts%unit) ts%normalized
    read(ts%unit) ts%scale
    read(ts%unit) ts%shift

  end subroutine ts_read_header

  !--------------------------------------------------------------------!

  subroutine ts_skip_header(ts)

    implicit none

    type(TrnSet), intent(inout) :: ts

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    read(ts%unit) ! ts%nTypes
    read(ts%unit) ! ts%nStrucs
    read(ts%unit) ! ts%typeName(:)
    read(ts%unit) ! ts%E_atom(:)
    read(ts%unit) ! ts%normalized
    read(ts%unit) ! ts%scale
    read(ts%unit) ! ts%shift

  end subroutine ts_skip_header

  !--------------------------------------------------------------------!
  !           training set file footer containing statistics           !
  !--------------------------------------------------------------------!

  subroutine ts_write_footer(ts, stp)

    implicit none

    type(TrnSet),                        intent(inout) :: ts
    type(Setup), dimension(:), optional, intent(in)    :: stp

    integer :: itype, nTypes
    logical :: has_setups

    call ts_assert_init(ts)
    call ts_assert_writemode(ts)

    if (ts%iStruc < ts%nStrucs) then
       write(0,*) "Warning: writing footer to incomplete training set file."
    end if

    write(ts%unit) ts%nAtomsTot
    write(ts%unit) ts%E_av, ts%E_min, ts%E_max
    if (present(stp)) then
       nTypes = size(stp(:))
       if (nTypes /= ts%nTypes) then
          write(0,*) "Error: wrong size of array stp in `ts_read_footer()'."
          stop
       end if
       has_setups = .true.
       write(ts%unit) has_setups
       do itype = 1, ts%nTypes
          write(ts%unit) itype
          call save_Setup(stp(itype), unit=ts%unit)
       end do
    else
       has_setups = .false.
       write(ts%unit) has_setups
    end if

  end subroutine ts_write_footer

  !--------------------------------------------------------------------!

  subroutine ts_read_footer(ts, stp)

    implicit none

    type(TrnSet),                        intent(inout) :: ts
    type(Setup), dimension(:), optional, intent(out)   :: stp

    integer :: i, itype, nTypes
    logical :: has_setups

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    do while(ts%iStruc < ts%nStrucs)
       call ts_skip_structure(ts)
    end do

    read(ts%unit) ts%nAtomsTot
    read(ts%unit) ts%E_av, ts%E_min, ts%E_max
    read(ts%unit) has_setups

    if (present(stp) .and. has_setups) then
       nTypes = size(stp(:))
       if (nTypes /= ts%nTypes) then
          write(0,*) "Error: wrong size of array stp in `ts_read_footer()'."
          stop
       end if
       do i = 1, nTypes
          read(ts%unit) itype
          stp(itype) = load_Setup(ts%typeName, unit=ts%unit)
       end do
    else if (present(stp) .and. (.not. has_setups)) then
       write(0,*) "Error: no structural fingerprint basis setups in training set file!"
       stop
    end if

  end subroutine ts_read_footer

  !--------------------------------------------------------------------!
  !              data from structures in the training set              !
  !--------------------------------------------------------------------!

  subroutine ts_write_structure_info(ts, filename, nAtoms, nTypes, energy)

    implicit none

    type(TrnSet),     intent(inout) :: ts
    character(len=*), intent(in)    :: filename
    integer,          intent(in)    :: nAtoms
    integer,          intent(in)    :: nTypes
    double precision, intent(in)    :: energy

    double precision :: E_atom

    call ts_assert_init(ts)
    call ts_assert_writemode(ts)

    if (ts%iStruc >= ts%nStrucs) then
       write(0,*) "Error: too many files for training set."
       stop
    else
       ts%iStruc = ts%iStruc + 1
    end if

    write(ts%unit) len_trim(filename)
    write(ts%unit) trim(filename)
    write(ts%unit) nAtoms, nTypes
    write(ts%unit) energy

    ! energy stats
    E_atom = energy/dble(nAtoms)
    if (ts%iStruc > 1) then
       ts%E_min = min(ts%E_min, E_atom)
       ts%E_max = max(ts%E_max, E_atom)
       ts%E_av  = ts%E_av + E_atom/dble(ts%nStrucs)
    else
       ts%E_min = E_atom
       ts%E_max = E_atom
       ts%E_av  = E_atom/dble(ts%nStrucs)
    end if

    ! keep track of the atoms in the training set
    ts%nAtomsTot = ts%nAtomsTot + nAtoms

  end subroutine ts_write_structure_info

  !--------------------------------------------------------------------!

  subroutine ts_read_structure_info(ts, filename, nAtoms, nTypes, energy)

    implicit none

    type(TrnSet),     intent(inout) :: ts
    character(len=*), intent(out)   :: filename
    integer,          intent(out)   :: nAtoms
    integer,          intent(out)   :: nTypes
    double precision, intent(out)   :: energy

    integer :: l

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    if (ts%iStruc >= ts%nStrucs) then
       write(0,*) "Error: no more file record to read."
       stop
    else
       ts%iStruc = ts%iStruc + 1
    end if

    read(ts%unit) l
    filename = ' '
    l = min(l,len(filename))
    read(ts%unit) filename(1:l)
    read(ts%unit) nAtoms, nTypes
    read(ts%unit) energy

  end subroutine ts_read_structure_info

  !--------------------------------------------------------------------!

  subroutine ts_write_atom_info(ts, itype, cooCart, forCart)

    implicit none

    type(TrnSet),                   intent(inout) :: ts
    integer,                        intent(in)    :: itype
    double precision, dimension(3), intent(in)    :: cooCart
    double precision, dimension(3), intent(in)    :: forCart

    call ts_assert_init(ts)
    call ts_assert_writemode(ts)

    write(ts%unit) itype
    write(ts%unit) cooCart(1:3)
    write(ts%unit) forCart(1:3)

  end subroutine ts_write_atom_info

  !--------------------------------------------------------------------!

  subroutine ts_read_atom_info(ts, itype, cooCart, forCart)

    implicit none

    type(TrnSet),                   intent(inout) :: ts
    integer,                        intent(out)   :: itype
    double precision, dimension(3), intent(out)   :: cooCart
    double precision, dimension(3), intent(out)   :: forCart

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    read(ts%unit) itype
    read(ts%unit) cooCart(1:3)
    read(ts%unit) forCart(1:3)

  end subroutine ts_read_atom_info

  !--------------------------------------------------------------------!

  subroutine ts_write_sf_info(ts, nsf, sfval)

    implicit none

    type(TrnSet),                       intent(inout) :: ts
    integer,                            intent(in)    :: nsf
    double precision, dimension(nsf),   intent(in)    :: sfval

    call ts_assert_init(ts)
    call ts_assert_writemode(ts)

    write(ts%unit) nsf
    write(ts%unit) sfval(1:nsf)

  end subroutine ts_write_sf_info

  !--------------------------------------------------------------------!

  subroutine ts_read_sf_info(ts, nsf)

    implicit none

    type(TrnSet), intent(inout) :: ts
    integer,      intent(out)   :: nsf

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    read(ts%unit) nsf

  end subroutine ts_read_sf_info

  !--------------------------------------------------------------------!

  subroutine ts_read_sf_values(ts, nsf, sfval)

    implicit none

    type(TrnSet),                       intent(inout) :: ts
    integer,                            intent(in)    :: nsf
    double precision, dimension(nsf),   intent(out)   :: sfval

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    read(ts%unit) sfval(1:nsf)

  end subroutine ts_read_sf_values

  !--------------------------------------------------------------------!

  subroutine ts_skip_structure(ts)

    implicit none

    type(TrnSet), intent(inout) :: ts

    character(len=1024) :: filename
    integer             :: nAtoms
    integer             :: nTypes
    double precision    :: energy

    ! sanity checks & increment of iFile done by read_structure_info
    call ts_read_structure_info(ts, filename, nAtoms, nTypes, energy)

    ! skip atomic data
    call ts_skip_atoms(ts, nAtoms)

  end subroutine ts_skip_structure

  !--------------------------------------------------------------------!

  subroutine ts_skip_atoms(ts, nAtoms)

    implicit none

    type(TrnSet), intent(inout) :: ts
    integer,      intent(in)    :: nAtoms

    integer :: iatom

    do iatom = 1, nAtoms
       read(ts%unit) ! itype
       read(ts%unit) ! cooCart(1:3)
       read(ts%unit) ! forCart(1:3)
       read(ts%unit) ! stp%nsf
       read(ts%unit) ! stp%sfval(1:stp%nsf)
       !$$ read(ts%unit) ! stp%sfderiv_i(1:3,1:stp%nsf)
    end do

  end subroutine ts_skip_atoms

  !--------------------------------------------------------------------!
  !                 save/restore basis function set-ups                !
  !--------------------------------------------------------------------!

  subroutine ts_load_Setups(ts, stp)

    implicit none

    type(TrnSet),              intent(inout) :: ts
    type(Setup), dimension(:), intent(out)   :: stp

    call ts_read_footer(ts, stp=stp)

  end subroutine ts_load_Setups

  !--------------------------------------------------------------------!

  subroutine ts_unload_Setups(ts, stp)

    implicit none

    type(TrnSet),              intent(inout) :: ts
    type(Setup), dimension(:), intent(inout) :: stp

    integer :: itype

    call ts_assert_init(ts)

    do itype = 1, ts%nTypes
       if (stp(itype)%init) call del_Setup(stp(itype))
    end do

  end subroutine ts_unload_Setups

  !--------------------------------------------------------------------!
  !                       data set normalization                       !
  !--------------------------------------------------------------------!

  subroutine ts_normalize(ts, maxenergy)

    implicit none

    type(TrnSet),     intent(inout) :: ts
    double precision, intent(in)    :: maxenergy

    type(TrnSet)                                :: ts2
    character(len=1024)                         :: file2
    type(Setup),      dimension(:), allocatable :: stp
    double precision, dimension(:), allocatable :: sfval

    double precision               :: scale, shift
    character(len=1024)            :: filename
    integer                        :: iatom, nAtoms
    integer                        :: itype, nTypes
    integer                        :: nStrucs2
    integer                        :: nsf, nsf_max
    double precision               :: energy
    double precision, dimension(3) :: cooCart, forCart

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    ! file name for normalized training set
    file2 = trim(ts%file) // ".scaled"

    call aeio_header('Training set normalization')
    write(*,*)
    write(*,*) 'The training set will be normalized now.  Depending on its size this'
    write(*,*) 'process can take a while.  The normalized data set will be written to'
    write(*,*) 'another file. Load that file in future to avoid this step.'
    write(*,*)
    write(*,*) 'Name of the new training set file: ' // trim(file2)
    write(*,*)

    allocate(stp(ts%nTypes))

    call ts_read_footer(ts, stp=stp)
    call rewind_TrnSet(ts)

    !$ ! initialize the setup module to allocate needed memory
    !$ ! however, we set nnb_max=1, as more is not necessary
    !$ call stp_init(ts%nTypes, stp, 1)
    nsf_max = stp_nsf_max(stp=stp)
    allocate(sfval(nsf_max))
    sfval(:) = 0.0d0

    ! Further normalization methods could be implemented here
    call ts_energy_norm_simple(ts, maxenergy, scale, shift)
    write(*,*) 'The network output energy will be normalized to the interval [-1,1].'
    write(*,*) '  Energy scaling factor: f = ' // trim(io_adjustl(scale,6))
    write(*,*) '  Atomic energy shift  : s = ' // trim(io_adjustl(shift,6))
    write(*,*)

    ! If maxenergy is lower than some structure in the initial training
    ! set, these structures will not be included.  Thus, we need to count
    ! all included structures.
    if (maxenergy < ts%E_max) then
       nStrucs2 = 0
       do while(ts%iStruc < ts%nStrucs)
          call ts_read_structure_info(ts, filename, nAtoms, nTypes, energy)
          energy = scale*(energy - dble(nAtoms)*shift)
          if (energy/dble(nAtoms) <= 1.001d0) nStrucs2 = nStrucs2 + 1
          call ts_skip_atoms(ts, nAtoms)
       end do
       call rewind_TrnSet(ts)
       if (nStrucs2 <= ts%nStrucs) then
          write(*,*) trim(io_adjustl(ts%nStrucs-nStrucs2)) // &
               ' high-energy structures will be removed from the scaled training set.'
          write(*,*)
       end if
    else
       nStrucs2 = ts%nStrucs
    end if

    ! new normalized (and possibly smaller) training set
    ts2 = new_TrnSet(ts%nTypes, ts%typeName, ts%E_atom, nStrucs2, &
                     trim(file2), scale=scale, shift=shift)

    do while(ts%iStruc < ts%nStrucs)
       call ts_read_structure_info(ts, filename, nAtoms, nTypes, energy)
       energy = scale*(energy - dble(nAtoms)*shift)
       if (energy/dble(nAtoms) <= 1.001d0) then
          call ts_write_structure_info(ts2, filename, nAtoms, nTypes, energy)
          do iatom = 1, nAtoms
             call ts_read_atom_info(ts, itype, cooCart, forCart)
             ! scale force to match energy
             forCart(1:3) = scale*forCart(1:3)
             call ts_write_atom_info(ts2, itype, cooCart, forCart)
             call ts_read_sf_info(ts, nsf)
             call ts_read_sf_values(ts, nsf, sfval(1:nsf))
             call stp_normalize(stp(itype), sfval)
             call ts_write_sf_info(ts2, nsf, sfval(1:nsf))
          end do
       else
          ! energy too high, do not include this structure
          call ts_skip_atoms(ts, nAtoms)
       end if
    end do

    call stp_final(ts%nTypes, stp)
    call close_TrnSet(ts2, stp=stp)
    call ts_unload_Setups(ts,stp)
    deallocate(stp, sfval)
    call close_TrnSet(ts)
    ts = open_TrnSet(trim(file2))

    write(*,*) 'Training set normalization done.'
    write(*,*)

  end subroutine ts_normalize

  !--------------------------------------------------------------------!

  subroutine ts_energy_norm_simple(ts, maxenergy, scale, shift)

    implicit none

    type(TrnSet),     intent(inout) :: ts
    double precision, intent(in)    :: maxenergy
    double precision, intent(out)   :: scale
    double precision, intent(out)   :: shift

    double precision :: E_max

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    E_max = min(ts%E_max, maxenergy)

    ! scale = 1.0d0/(max(abs(ts%E_max), abs(ts%E_min)))
    scale =  2.0d0/(E_max - ts%E_min)
    shift =  0.5d0*(E_max + ts%E_min)

  end subroutine ts_energy_norm_simple

  !--------------------------------------------------------------------!
  !                 count atoms in train and test set                  !
  !--------------------------------------------------------------------!

  subroutine ts_count_atoms(ts, isTest, nAtomsTrain, nAtomsTest)

    implicit none

    type(TrnSet),          intent(inout) :: ts
    logical, dimension(:), intent(in)    :: isTest
    integer,               intent(out)   :: nAtomsTrain
    integer,               intent(out)   :: nAtomsTest

    character(len=1024) :: filename
    integer             :: iatom, nAtoms
    integer             :: nTypes
    double precision    :: energy

    call ts_assert_init(ts)
    call ts_assert_readmode(ts)

    if (size(isTest) < ts%nStrucs) then
       write(0,*) "Error: vector `isTest' too short in `ts_count_atoms()'."
       stop
    end if

    if (ts%iStruc > 0) call rewind_TrnSet(ts)

    nAtomsTrain = 0
    nAtomsTest  = 0
    do while(ts%iStruc < ts%nStrucs)
       call ts_read_structure_info(ts, filename, nAtoms, nTypes, energy)
       do iatom = 1, nAtoms
          read(ts%unit) ! itype
          read(ts%unit) ! cooCart(1:3)
          read(ts%unit) ! forCart(1:3)
          read(ts%unit) ! stp%nsf
          read(ts%unit) ! stp%sfval(1:stp%nsf)
          !$$ read(ts%unit) ! stp%sfderiv_i(1:3,1:stp%nsf)
       end do
       if (isTest(ts%iStruc)) then
          nAtomsTest = nAtomsTest + nAtoms
       else
          nAtomsTrain = nAtomsTrain + nAtoms
       end if
    end do

  end subroutine ts_count_atoms

  !--------------------------------------------------------------------!
  !                            state checks                            !
  !--------------------------------------------------------------------!

  subroutine ts_assert_init(ts)
    implicit none
    type(TrnSet), intent(in) :: ts
    if (.not. ts%init) then
       write(0,*) "Error: training set not initialized."
       stop
    end if
  end subroutine ts_assert_init

  subroutine ts_assert_notinit(ts)
    implicit none
    type(TrnSet), intent(in) :: ts
    if (ts%init) then
       write(0,*) "Error: training set already initialized."
       stop
    end if
  end subroutine ts_assert_notinit

  subroutine ts_assert_writemode(ts)
    implicit none
    type(TrnSet), intent(in) :: ts
    if (trim(ts%mode) /= 'write') then
       write(0,*) "Error: training set not in 'write' mode."
       stop
    end if
  end subroutine ts_assert_writemode

  subroutine ts_assert_readmode(ts)
    implicit none
    type(TrnSet), intent(in) :: ts
    if (trim(ts%mode) /= 'read') then
       write(0,*) "Error: training set not in 'read' mode."
       stop
    end if
  end subroutine ts_assert_readmode

end module trainset
