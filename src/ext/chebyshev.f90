!-----------------------------------------------------------------------
!     chebyshev.f90 - Chebyshev polynomials and their derivatives
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
! 2015-12-27 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

module chebyshev

  implicit none
  private

  public :: chebyshev_polynomial, &
            chebyshev_polynomial_d1

contains

  !--------------------------------------------------------------------!
  !                   Evaluate Chebyshev polynomials                   !
  !--------------------------------------------------------------------!

  function chebyshev_polynomial(r, r0, r1, n) result(T)
    ! Arguments:
    !    r        function argument
    !    r0, r1   the Chebyshev polynomials will be rescaled from [-1,1]
    !             to the interval [r0,r1]
    !    n        maximum polynomial order
    ! Returns:
    !    T(i)  with i=1,n+1  where T(i) is the Chebyshev polynomial of
    !    order (i-1)
    !
    ! The Chebyshev polynomials obey the following recurrence relation:
    !    T[0](x) = 1
    !    T[1](x) = x
    !    T[n+1](x) = 2x T[n](x) - T[n-1](x)

    implicit none

    double precision, intent(in)     :: r, r0, r1
    integer,          intent(in)     :: n
    double precision, dimension(n+1) :: T

    integer          :: i
    double precision :: x

    x = (2.0d0*r - r0 - r1)/(r1 - r0)

    T(1) = 1.0d0
    if (n > 0) then
       T(2) = x
       do i = 3, n+1
          T(i) = 2.0d0*x*T(i-1) - T(i-2)
       end do
    end if

  end function chebyshev_polynomial

  !--------------------------------------------------------------------!
  !             First derivative of Chebyshev polynomials              !
  !--------------------------------------------------------------------!

  function chebyshev_polynomial_d1(r, r0, r1, n) result(dT)
    ! Arguments:
    !    r        function argument
    !    r0, r1   the Chebyshev polynomials will be rescaled from [-1,1]
    !             to the interval [r0,r1]
    !    n        maximum polynomial order
    ! Returns:
    !    dT(i)  with i=1,n+1  where dT(i) is the first derivative of the
    !    Chebyshev polynomial of order (i-1)
    !
    ! The derivatives of the Chebyshev polynomials obey the following
    ! recurrence relation:
    !    dT[n](x)/dx = n U[n-1](x) with n = 1,...
    ! where
    !    U[0](x) = 1
    !    U[1](x) = 2x
    !    U[n+1](x) = 2x U[n](x) - U[n-1](x)
    ! are the Chebyshev polynomials of the second kind.

    implicit none

    double precision, intent(in)     :: r, r0, r1
    integer,          intent(in)     :: n
    double precision, dimension(n+1) :: dT

    integer          :: i
    double precision :: x
    double precision :: U1, U2, U3

    x = (2.0d0*r - r0 - r1)/(r1 - r0)

    dT(1) = 0.0d0
    if (n > 0) then
       U1 = 1.0d0
       dT(2) = U1
       U2 = 2.0d0*x
       do i = 3, n+1
          dT(i) = U2*dble(i-1)
          U3 = 2.0d0*x*U2 - U1
          U1 = U2
          U2 = U3
       end do
    end if

    ! inner derivative (from rescaling)
    dT = dT*2.0d0/(r1 - r0)

  end function chebyshev_polynomial_d1

end module chebyshev

!----------------------------------------------------------------------!
!                              Unit test                               !
!----------------------------------------------------------------------!

!!$ program chebyshev_test
!!$
!!$   use chebyshev, only: chebyshev_polynomial, chebyshev_polynomial_d1
!!$
!!$   implicit none
!!$
!!$   integer, parameter          :: N = 100
!!$   integer, parameter          :: O = 10
!!$   double precision, parameter :: R0 = -2.0d0
!!$   double precision, parameter :: R1 =  5.0d0
!!$
!!$   double precision, dimension(O+1) :: T, dT, T_prev
!!$
!!$   double precision    :: r, dr
!!$   integer             :: i
!!$   character(len=1024) :: frmt
!!$
!!$   write(frmt, *) O + 2
!!$   frmt = "(1x," // trim(adjustl(frmt)) // "(ES18.8,1x))"
!!$
!!$   open(20, file="chebytest-values.dat", status="replace", action="write")
!!$   open(21, file="chebytest-derivs.dat", status="replace", action="write")
!!$   open(22, file="chebytest-nderivs.dat", status="replace", action="write")
!!$
!!$   dr = (R1 - R0)/dble(N - 1)
!!$   r = r0
!!$   T_prev = 0.0d0
!!$   do i = 1, N
!!$      T = chebyshev_polynomial(r, R0, R1, O)
!!$      dT = chebyshev_polynomial_d1(r, R0, R1, O)
!!$      write(20, frmt) r, T(1:O+1)
!!$      write(21, frmt) r, dT(1:O+1)
!!$      if (i > 1) then
!!$         write(22, frmt) r - 0.5d0*dr, (T - T_prev)/dr
!!$      end if
!!$      T_prev = T
!!$      r = r + dr
!!$   end do
!!$
!!$   close(20)
!!$   close(21)
!!$   close(22)
!!$
!!$ end program chebyshev_test
