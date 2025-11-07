!-----------------------------------------------------------------------
!       potential.f90 - handling of ANN potential parameter files
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
! 2013-05-27 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

module potential

  use aeio,        only: PATHLEN, TYPELEN, &
                         aeio_assert_file_exists

  use geometry,    only: geo_itype_of_name

  use io,          only: io_unit

  use feedforward, only: Network,          &
                         load_Network,     &
                         del_Network,      &
                         ff_print_info

  use sfsetup,     only: Setup,            &
                         load_Setup,       &
                         del_Setup,        &
                         stp_init,         &
                         stp_final,        &
                         stp_print_info

  use trainset,    only: TrnSet,           &
                         load_TrnSet_info, &
                         close_TrnSet,     &
                         ts_print_info

  implicit none
  save

  public :: load_NNPot,       &
            load_NNPot_ASCII, &
            del_NNPot,        &
            pot_init,         &
            pot_final,        &
            pot_get_range,    &
            pot_print_info,   &
            pot_assert_init

  !--------------------------------------------------------------------!
  !                      Neural Network potential                      !
  !--------------------------------------------------------------------!

  type, public :: NNPot

     !-----------------------------------------------------------------!
     ! init        .true. if potential has been initialized            !
     ! file        path to the potential file, if known                !
     ! unit        unit number, if not directly loaded from file 'unit'!
     !             will be set to -1, if the potential was loaded      !
     !             directly from a file                                !
     ! typeName    species of the central atom                         !
     !                                                                 !
     ! E_scale     energy scaling factor                               !
     ! E_shift     shift of the atomic energy                          !
     ! E_atom      atomic reference energy of the central atom         !
     !             This is the shifted atomic energy!  If you need the !
     !             reference atmic energy, use E_atom - E_shift .      !
     !                                                                 !
     ! stp         structural fingerprint basis setup (type: Setup)    !
     ! net         trained neural network (type: Network)              !
     ! ts          training set info (type: TrnSet)                    !
     !-----------------------------------------------------------------!

     logical                :: init = .false.

     character(len=PATHLEN) :: file
     integer                :: unit

     character(len=TYPELEN) :: typeName
     double precision       :: E_scale
     double precision       :: E_shift
     double precision       :: E_atom

     type(Setup)            :: stp
     type(Network)          :: net
     type(TrnSet)           :: ts

  end type NNPot

  !--------------------------------------------------------------------!

  logical :: isInit = .false.

contains

  function load_NNPot(global_types, file, unit) result(pot)

    implicit none

    character(len=*), dimension(:), intent(in) :: global_types
    character(len=*), optional,     intent(in) :: file
    integer,          optional,     intent(in) :: unit
    type(NNPot)                                :: pot

    integer :: itype
    integer :: u

    if (.not. (present(file) .or. present(unit))) then
       write(0,*) "Error: neither file nor file unit specified in `load_NNPot()'."
       stop
    end if

    if (present(unit)) then
       u = unit
       pot%unit = u
       pot%file = ''
    else
       u = io_unit()
       call aeio_assert_file_exists(file)
       open(u, file=trim(adjustl(file)), status='old', action='read', &
            form='unformatted')
       pot%file = trim(adjustl(file))
       pot%unit = -1
    end if

    pot%net = load_Network(unit=u)
    pot%stp = load_Setup(global_types, unit=u)
    pot%ts  = load_TrnSet_info(unit=u)

    pot%E_scale     = 1.0d0/pot%ts%scale

    ! central atom
    pot%typeName = pot%stp%atomtype
    pot%E_shift  = pot%ts%shift
    itype = geo_itype_of_name(pot%typeName, pot%ts%typeName)
    pot%E_atom   = pot%ts%E_atom(itype)

    if (.not. present(unit)) close(u)

    pot%init = .true.

  end function load_NNPot

  !--------------------------------------------------------------------!

  subroutine del_NNPot(pot)

    implicit none

    type(NNPot), intent(inout) :: pot

    if (.not. pot%init) return

    call close_TrnSet(pot%ts)
    call del_Setup(pot%stp)
    call del_Network(pot%net)

    pot%init = .false.

  end subroutine del_NNPot

  !--------------------------------------------------------------------!
  !             load ANN potential from ASCII format                   !
  !                                                                    !
  ! This subroutine uses the ASCII format introduced by                !
  ! Jon Lopez-Zorilla's aenet-PyTorch training implementation.  The    !
  ! code is partially based on aenet-PyTorch's `nnASCII2bin.f90' tool, !
  ! which is released under the terms of the MIT license with          !
  ! copyright held by 2023 J. Lopez-Zorrilla (jon.lopezz@ehu.eus) and  !
  ! N. Artrith (n.artrith@uu.nl).  The aenet-PyTorch code can be       !
  ! obtained from https://github.com/atomisticnet/aenet-pytorch        !
  ! When using the ASCII format, please cite:                          !
  ! J. Lopez-Zorrilla et al., J. Chem. Phys. 158, 164105 (2023).       !
  ! DOI: https://doi.org/10.1063/5.0146803                             !
  !--------------------------------------------------------------------!

  function load_NNPot_ASCII(global_types, file, unit, output_file) result(pot)

    implicit none

    character(len=*), dimension(:), intent(in) :: global_types
    character(len=*), optional,     intent(in) :: file
    integer,          optional,     intent(in) :: unit
    character(len=*), optional,     intent(in) :: output_file
    type(NNPot)                                :: pot

    character(len=PATHLEN) :: outfile

    integer :: itype
    integer :: u_in, u_out

    integer              :: nlayers, nnodesmax, Wsize, nvalues
    integer, allocatable :: nnodes(:), fun(:), iw(:), iv(:)
    double precision, allocatable  :: W(:)

    character(len=1024)           :: description
    character(len=100)            :: sftype
    character(len=2)              :: atomtype
    character(len=2), allocatable :: envtypes(:)
    double precision              :: rc_min, rc_max
    integer                       :: nsf, nsfparam, neval, nenv
    integer, allocatable          :: sf(:), sfenv(:,:)
    double precision, allocatable :: sfparam(:,:), sfval_min(:)
    double precision, allocatable :: sfval_max(:), sfval_avg(:), sfval_cov(:)
    character(len=PATHLEN)        :: filename
    logical                       :: normalized
    double precision              :: scale, shift, E_min, E_max, E_avg
    integer                       :: ntypes, natomtot, nstrucs
    character(len=2), allocatable :: type_names(:)
    double precision, allocatable :: E_atom(:)

    if (.not. (present(file) .or. present(unit))) then
       write(0,*) "Error: neither file nor file unit " &
            // "specified in `load_NNPot_ascii()'."
       stop
    end if

    if (present(unit)) then
       u_in = unit
    else
       call aeio_assert_file_exists(file)
       u_in = io_unit()
       open(u_in, action="read", status="old", file=file)
    end if

    ! Network information
    read(u_in,*) nlayers
    read(u_in,*) nnodesmax
    read(u_in,*) Wsize
    read(u_in,*) nvalues

    allocate(nnodes(nlayers), fun(nlayers-1), &
         iw(nlayers), iv(nlayers), W(Wsize))

    read(u_in,*) nnodes(:)
    read(u_in,*) fun(:)
    read(u_in,*) iw(:)
    read(u_in,*) iv(:)
    read(u_in,*) W(:)

    ! Structural Fingerprint setup information
    read(u_in,*) description
    read(u_in,*) atomtype
    read(u_in,*) nenv

    allocate(envtypes(nenv))

    read(u_in,*) envtypes(:)
    read(u_in,*) rc_min
    read(u_in,*) rc_max
    read(u_in,*) sftype
    read(u_in,*) nsf
    read(u_in,*) nsfparam

    allocate(sf(nsf), sfparam(nsfparam,nsf), sfenv(2,nsf), sfval_min(nsf), &
             sfval_max(nsf), sfval_avg(nsf), sfval_cov(nsf))

    read(u_in,*) sf(:)
    read(u_in,*) sfparam(:,:)
    read(u_in,*) sfenv(:,:)
    read(u_in,*) neval
    read(u_in,*) sfval_min
    read(u_in,*) sfval_max
    read(u_in,*) sfval_avg
    read(u_in,*) sfval_cov

    ! Trainset information
    read(u_in,*) filename
    read(u_in,*) normalized
    read(u_in,*) scale
    read(u_in,*) shift
    read(u_in,*) ntypes

    allocate(type_names(ntypes), E_atom(ntypes))

    read(u_in,*) type_names(:)
    read(u_in,*) E_atom(:)
    read(u_in,*) natomtot
    read(u_in,*) nstrucs
    read(u_in,*) E_min, E_max, E_avg

    if (.not. present(unit)) close(u_in)

    if (present(output_file)) then
       outfile = trim(adjustl(output_file))
    else
       outfile = trim(adjustl(atomtype)) // ".nn"
    end if

    ! open binary/unformatted output file
    u_out = io_unit()
    open(u_out, action="write", status="replace", file=outfile, &
         form="unformatted")

    ! network
    write(u_out) nlayers
    write(u_out) nnodesmax
    write(u_out) Wsize
    write(u_out) nvalues
    write(u_out) nnodes(:)
    write(u_out) fun(:)
    write(u_out) iw(:)
    write(u_out) iv(:)
    write(u_out) W(:)

    ! sf setup
    write(u_out) description
    write(u_out) atomtype
    write(u_out) nenv
    write(u_out) envtypes(:)
    write(u_out) rc_min
    write(u_out) rc_max
    write(u_out) sftype
    write(u_out) nsf
    write(u_out) nsfparam
    write(u_out) sf(:)
    write(u_out) sfparam(:,:)
    write(u_out) sfenv(:,:)
    write(u_out) neval
    write(u_out) sfval_min
    write(u_out) sfval_max
    write(u_out) sfval_avg
    write(u_out) sfval_cov

    ! trainset info
    write(u_out) filename
    write(u_out) normalized
    write(u_out) scale
    write(u_out) shift
    write(u_out) ntypes
    write(u_out) type_names(:)
    write(u_out) E_atom(:)
    write(u_out) natomtot
    write(u_out) nstrucs
    write(u_out) E_min, E_max, E_avg

    close(u_out)

    deallocate(nnodes, fun, iw, iv, W)
    deallocate(sf, sfparam, sfenv, sfval_min, sfval_max, &
         sfval_avg, sfval_cov)
    deallocate(type_names, E_atom)

    pot = load_NNPot(global_types, file=outfile)

  end function load_NNPot_ASCII


  !--------------------------------------------------------------------!
  !          initialize memory for a collection of potentials          !
  !--------------------------------------------------------------------!

  subroutine pot_init(nTypes, pot, nnb_max)

    implicit none

    integer,                                intent(in) :: nTypes
    type(NNPot), dimension(nTypes), target, intent(in) :: pot
    integer,                                intent(in) :: nnb_max

    if (isInit) then
       write(0,*) "Warning: module already initialized in `pot_init()'"
       return
    end if

    call stp_init(nTypes, pot(1:nTypes)%stp, nnb_max)

    isInit = .true.

  end subroutine pot_init

  !--------------------------------------------------------------------!

  subroutine pot_final(nTypes, pot)

    implicit none

    integer,                                intent(in) :: nTypes
    type(NNPot), dimension(nTypes), target, intent(in) :: pot

    if (.not. isInit) return

    call stp_final(nTypes, pot(1:nTypes)%stp)

    isInit = .false.

  end subroutine pot_final

  !--------------------------------------------------------------------!
  !                   print info about NN potential                    !
  !--------------------------------------------------------------------!

  subroutine pot_print_info(pot)

    implicit  none

    type(NNPot), intent(inout) :: pot

    call pot_assert_init(pot)

    write(*,*) 'Atomic species : ', trim(pot%typeName)
    write(*,*) 'File name      : ', trim(pot%file)
    write(*,*)

    call ts_print_info(pot%ts)
    call ff_print_info(pot%net)
    call stp_print_info(pot%stp)

  end subroutine pot_print_info

  !--------------------------------------------------------------------!
  !      determine interaction range for over several potentials       !
  !                (see also stp_get_range in sfsetup)                 !
  !--------------------------------------------------------------------!

  subroutine pot_get_range(nTypes, pot, Rc_min, Rc_max)

    implicit none

    integer,                        intent(in)  :: nTypes
    type(NNPot), dimension(nTypes), intent(in)  :: pot
    double precision,               intent(out) :: Rc_min
    double precision,               intent(out) :: Rc_max

    integer :: itype

    call pot_assert_init(pot(1))

    Rc_min = pot(1)%stp%Rc_min
    Rc_max = pot(1)%stp%Rc_max

    do itype = 1, nTypes
       call pot_assert_init(pot(itype))
       Rc_min = min(Rc_min, pot(itype)%stp%Rc_min)
       Rc_max = max(Rc_max, pot(itype)%stp%Rc_max)
    end do

  end subroutine pot_get_range

  !--------------------------------------------------------------------!
  !                        auxiliary procedures                        !
  !--------------------------------------------------------------------!

  subroutine pot_assert_init(pot)
    implicit none
    type(NNPot), intent(in) :: pot
    if (.not. pot%init) then
       write(*,*) 'Error: NN potential is not initialized.'
       write(*,*)
       stop
    end if
  end subroutine pot_assert_init

end module potential
