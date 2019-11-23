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

  public :: load_NNPot,     &
            del_NNPot,      &
            pot_init,       &
            pot_final,      &
            pot_get_range,  &
            pot_print_info, &
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
