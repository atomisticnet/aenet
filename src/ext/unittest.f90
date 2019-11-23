!-----------------------------------------------------------------------
!               testing.f90 -- routines for unit testing
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
! 2014-09-24 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

module unittest

  use io, only: io_trim

  implicit none
  private

  public :: tst_new, tst_check_passed, tst_equal, tst_assert

  integer, public :: tst_msg_len = 61
  logical, public :: has_passed

  logical, private :: has_EPS = .false.

  interface tst_equal
     module procedure tst_equal_d, tst_equal_dn,    &
                      tst_equal_dn2, tst_equal_dn3, &
                      tst_equal_i, tst_equal_in,    &
                      tst_equal_in2, tst_equal_c,   &
                      tst_equal_cn, tst_equal_cn2
  end interface tst_equal

contains

  !--------------------------------------------------------------------!
  !           message describing the test without line break           !
  !--------------------------------------------------------------------!

  subroutine tst_new(msg)

    implicit none

    character(len=*), intent(in) :: msg

    write(*,'(1x,A,": ")', advance='no') io_trim(msg,tst_msg_len)
    has_passed = .true.

  end subroutine tst_new

  !--------------------------------------------------------------------!
  !           print 'passed.' or 'FAILED!' depending on test           !
  !--------------------------------------------------------------------!

  subroutine tst_check_passed(check)

    implicit none

    logical, optional, intent(in) :: check

    if (present(check)) has_passed = (has_passed .and. check)

    if (has_passed) then
       write(*,*) 'passed.'
    else
       write(*,*) 'FAILED!'
    end if

  end subroutine tst_check_passed

  !--------------------------------------------------------------------!
  !                           assert or fail                           !
  !--------------------------------------------------------------------!

  subroutine tst_assert(condition, message, status)

    implicit none

    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    logical, optional, intent(out) :: status

    if (.not. condition) then
       write(0, *) "FAILED assertion: " // trim(message)
       if (present(status)) then
          status = condition
       else
          stop
       end if
    end if

  end subroutine tst_assert

  !--------------------------------------------------------------------!
  !                          assert equality                           !
  !--------------------------------------------------------------------!

  function tst_equal_i(a, b) result(is_equal)

    implicit none

    integer, intent(in) :: a, b
    logical             :: is_equal

    is_equal = (a == b)
    has_passed = (has_passed .and. is_equal)

  end function tst_equal_i

  !--------------------------------------------------------------------!

  function tst_equal_in(a, b) result(is_equal)

    implicit none

    integer, dimension(:), intent(in) :: a, b
    logical :: is_equal
    integer :: i, n

    is_equal = .true.
    if (size(a) /= size(b)) then
       is_equal = .false.
       has_passed = .false.
       return
    end if

    n = size(a)
    do i = 1, n
       is_equal = (is_equal .and. tst_equal_i(a(i), b(i)))
    end do

    has_passed = (has_passed .and. is_equal)

  end function tst_equal_in

  !--------------------------------------------------------------------!

  function tst_equal_in2(a, b) result(is_equal)

    implicit none

    integer, dimension(:,:), intent(in) :: a, b
    logical :: is_equal
    integer :: i, n

    is_equal = .true.
    if (any(shape(a) /= shape(b))) then
       is_equal = .false.
       has_passed = .false.
       return
    end if

    n = product(shape(a))
    is_equal = (is_equal .and. &
                tst_equal_in(reshape(a,[n]), reshape(b,[n])))

    has_passed = (has_passed .and. is_equal)

  end function tst_equal_in2

  !--------------------------------------------------------------------!

  function tst_equal_d(a, b, prec, tol) result(is_equal)

    implicit none

    double precision,           intent(in) :: a, b
    double precision, optional, intent(in) :: prec
    double precision, optional, intent(in) :: tol

    double precision :: abs_diff
    double precision :: largest
    double precision :: eps

    logical :: is_equal

    if (present(tol)) then
       eps = tol
    else
       ! eps = epsilon(1.0d0)*1000.0d0
       eps = 1.0d-10
    end if

    is_equal = (abs(a - b) <= abs(eps))

    if ((.not. is_equal) .and. present(prec)) then
       abs_diff = abs(a - b)
       largest = max(abs(a), abs(b))
       is_equal = (abs_diff <= largest*prec)
    end if

    has_passed = (has_passed .and. is_equal)

  end function tst_equal_d

  !--------------------------------------------------------------------!

  function tst_equal_dn(a, b, prec, tol) result(is_equal)

    implicit none

    double precision, dimension(:), intent(in) :: a, b
    double precision, optional,     intent(in) :: prec
    double precision, optional,     intent(in) :: tol

    integer :: i, n

    logical :: is_equal

    is_equal = .true.

    if (size(a) /= size(b)) then
       is_equal = .false.
       has_passed = .false.
       return
    end if

    n = size(a)
    do i = 1, n
       if (present(prec) .and. present(tol)) then
          is_equal = (is_equal .and. tst_equal_d(a(i), b(i), prec=prec, tol=tol))
       else if (present(prec)) then
          is_equal = (is_equal .and. tst_equal_d(a(i), b(i), prec=prec))
       else if (present(tol)) then
          is_equal = (is_equal .and. tst_equal_d(a(i), b(i), tol=tol))
       else
          is_equal = (is_equal .and. tst_equal_d(a(i), b(i)))
       end if
    end do

    has_passed = (has_passed .and. is_equal)

  end function tst_equal_dn

  !--------------------------------------------------------------------!

  function tst_equal_dn2(a, b, prec, tol) result(is_equal)

    implicit none

    double precision, dimension(:,:), intent(in) :: a, b
    double precision, optional,       intent(in) :: prec
    double precision, optional,       intent(in) :: tol

    integer :: n

    logical :: is_equal

    is_equal = .true.

    if (any(shape(a) /= shape(b))) then
       is_equal = .false.
       has_passed = .false.
       return
    end if

    n = product(shape(a))
    if (present(prec) .and. present(tol)) then
       is_equal = (is_equal .and. tst_equal_dn(&
            reshape(a,[n]), reshape(b,[n]), prec=prec, tol=tol))
    else if (present(prec)) then
       is_equal = (is_equal .and. tst_equal_dn(&
            reshape(a,[n]), reshape(b,[n]), prec=prec))
    else if (present(tol)) then
       is_equal = (is_equal .and. tst_equal_dn(&
            reshape(a,[n]), reshape(b,[n]), tol=tol))
    else
       is_equal = (is_equal .and. tst_equal_dn(&
            reshape(a,[n]), reshape(b,[n])))
    end if

    has_passed = (has_passed .and. is_equal)

  end function tst_equal_dn2

  !--------------------------------------------------------------------!

  function tst_equal_dn3(a, b, prec, tol) result(is_equal)

    implicit none

    double precision, dimension(:,:,:), intent(in) :: a, b
    double precision, optional,         intent(in) :: prec
    double precision, optional,         intent(in) :: tol

    integer :: n

    logical :: is_equal

    is_equal = .true.

    if (any(shape(a) /= shape(b))) then
       is_equal = .false.
       has_passed = .false.
       return
    end if

    n = product(shape(a))
    if (present(prec) .and. present(tol)) then
       is_equal = (is_equal .and. tst_equal_dn(&
            reshape(a,[n]), reshape(b,[n]), prec=prec, tol=tol))
    else if (present(prec)) then
       is_equal = (is_equal .and. tst_equal_dn(&
            reshape(a,[n]), reshape(b,[n]), prec=prec))
    else if (present(tol)) then
       is_equal = (is_equal .and. tst_equal_dn(&
            reshape(a,[n]), reshape(b,[n]), tol=tol))
    else
       is_equal = (is_equal .and. tst_equal_dn(&
            reshape(a,[n]), reshape(b,[n])))
    end if

    has_passed = (has_passed .and. is_equal)

  end function tst_equal_dn3

  !--------------------------------------------------------------------!

  function tst_equal_c(a, b) result(is_equal)

    implicit none

    character(len=*), intent(in) :: a, b
    logical             :: is_equal

    is_equal = (trim(a) == trim(b))
    has_passed = (has_passed .and. is_equal)

  end function tst_equal_c

  !--------------------------------------------------------------------!

  function tst_equal_cn(a, b) result(is_equal)

    implicit none

    character(len=*), dimension(:), intent(in) :: a, b
    logical :: is_equal
    integer :: i, n

    is_equal = .true.
    if (size(a) /= size(b)) then
       is_equal = .false.
       has_passed = .false.
       return
    end if

    n = size(a)
    do i = 1, n
       is_equal = (is_equal .and. tst_equal_c(a(i), b(i)))
    end do

    has_passed = (has_passed .and. is_equal)

  end function tst_equal_cn

  !--------------------------------------------------------------------!

  function tst_equal_cn2(a, b) result(is_equal)

    implicit none

    character(len=*), dimension(:,:), intent(in) :: a, b
    logical :: is_equal
    integer :: i, n

    is_equal = .true.
    if (any(shape(a) /= shape(b))) then
       is_equal = .false.
       has_passed = .false.
       return
    end if

    n = product(shape(a))
    is_equal = (is_equal .and. &
                tst_equal_cn(reshape(a,[n]), reshape(b,[n])))

    has_passed = (has_passed .and. is_equal)

  end function tst_equal_cn2

end module unittest
