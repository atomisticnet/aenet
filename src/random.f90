!-----------------------------------------------------------------------
!       random.f90 - everything related to pseudo random numbers
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
! This module is a simple interface to the default Fortran pseudo-RNG.
! There is no guarantee that the (compiler provided) RNG is 'good', so
! do not use this module for anything critical (such as MC simulations)
! unless you are certain that your compiler provides a reliable RNG.
!-----------------------------------------------------------------------
! 2013-07-03 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

module random

  implicit none
  private
  save

  public :: random_init,       &
            random_reinit,     &
            random_final,      &
            random_integer,    &
            random_gauss,      &
            random_save_state, &
            random_load_state


  double precision,        parameter, private :: PI = 3.141592653589793d0
  integer,                 parameter, private :: DEFAULT_UNIT = 67

  !--------------------------------------------------------------------!

  integer, dimension(:), allocatable, private :: iSeed

  double precision,                   private :: gauss_sigma = 1.0d0
  double precision,                   private :: gauss_x0    = 0.0d0
  integer,                            private :: gauss_N     = 20

  logical,                            private :: isInit = .false.

contains

  subroutine random_init()

    implicit none

    integer :: i, n
    integer :: clock

    if (isInit) return

    i = 0
    call random_seed(size=n)
    allocate(iSeed(n))
    call system_clock(count=clock)
    iSeed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put=iSeed)

    isInit = .true.

  end subroutine random_init

  !--------------------------------------------------------------------!

  subroutine random_reinit()

    implicit none

    call random_final()
    call random_init()

  end subroutine random_reinit

  !--------------------------------------------------------------------!

  subroutine random_final()

    implicit none

    if (.not. isInit) return

    deallocate(iSeed)
    isInit = .false.

  end subroutine random_final

  !--------------------------------------------------------------------!
  !                     save/load state of the RNG                     !
  !--------------------------------------------------------------------!

  subroutine random_save_state(file, unit)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer                            :: u
    integer                            :: n
    integer, dimension(:), allocatable :: curr_seed

    if (.not. (present(file) .or. present(unit))) then
       write(0,*) "Warning: neither file, nor unit given in 'save_state()'."
       return
    end if

    if (present(unit)) then
       u = unit
    else
       u = DEFAULT_UNIT
    end if

    if (present(file)) then
       open(u, file=trim(file), status='replace', action='write', &
               form='unformatted')
    end if

    call random_seed(size=n)
    allocate(curr_seed(n))
    call random_seed(get=curr_seed)
    write(u) n
    write(u) curr_seed
    deallocate(curr_seed)

    if (present(file)) close(u)

  end subroutine random_save_state

  !--------------------------------------------------------------------!

  subroutine random_load_state(file, unit)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit

    integer                            :: u
    integer                            :: n1, n2
    integer, dimension(:), allocatable :: curr_seed

    if (.not. (present(file) .or. present(unit))) then
       write(0,*) "Warning: neither file, nor unit given in 'load_state()'."
       return
    end if

    if (present(unit)) then
       u = unit
    else
       u = DEFAULT_UNIT
    end if

    if (present(file)) then
       open(u, file=trim(file), status='old', action='read', &
               form='unformatted')
    end if

    read(u) n1
    allocate(curr_seed(n1))
    read(u) curr_seed
    call random_seed(size=n2)
    if (n1 /= n2) then
       write(0,*) "Warning: unable to restart RNG state (incompatible seed)."
    else
       call random_seed(put=curr_seed)
    end if
    deallocate(curr_seed)

    if (present(file)) close(u)

  end subroutine random_load_state

  !--------------------------------------------------------------------!
  !            random integer irand with 1 <= irand <= imax            !
  !--------------------------------------------------------------------!

  function random_integer(imax) result(irand)

    implicit none

    integer, intent(in) :: imax
    integer             :: irand

    double precision :: r

    call random_number(r)
    irand = ceiling(r*dble(imax))
    ! in case (r == 0)
    if (irand == 0) irand = imax

  end function random_integer

  !--------------------------------------------------------------------!
  !                        normal distribution                         !
  !    (simple implementation based on the central limit theoreme)     !
  !--------------------------------------------------------------------!

  function random_gauss(sigma, x0, N) result(r)

    implicit none

    double precision, optional, intent(in) :: sigma
    double precision, optional, intent(in) :: x0
    integer,          optional, intent(in) :: N
    double precision                       :: r

    double precision :: x
    integer          :: i

    if (present(sigma)) gauss_sigma = sigma
    if (present(x0))    gauss_x0    = x0
    if (present(N))     gauss_N     = N

    r = 0.0d0
    do i = 1, gauss_N
       call random_number(x)
       r = r + x
    end do

    r = r - 0.5d0*dble(gauss_N)
    r = r*sqrt(12.0d0/dble(gauss_N))

    ! apply requested mean and variance
    r = gauss_x0 + gauss_sigma*r

  end function random_gauss



end module random
