!-----------------------------------------------------------------------
!         Unit tests for the BP symmetry functions module
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
! 2014-09-24 Alexander Urban (AU) and Nongnuch Artrith (NA)
!-----------------------------------------------------------------------
program test_symmfunc

  use io,       only: io_unlink
  use random,   only: random_init, random_final, random_integer
  use unittest, only: tst_new, tst_check_passed, tst_equal

  use symmfunc, only: sf_init,    &
                      sf_final,   &
                      sf_add_rad, &
                      sf_add_ang, &
                      sf_nG_type, &
                      sf_fingerprint

  implicit none

  call test_setup()
  call test_fcc_rotation()
  call test_derivatives()

contains

  subroutine test_setup()

    implicit none

    integer :: ntypes, nG
    logical :: has_passed

    call tst_new("Symmfunc Test 1: set-up")

    ntypes = 2
    nG = 1
    call sf_init(ntypes, nG)

    call sf_add_rad(1, 1, 1, Rc = 3.0d0)
    call sf_add_rad(2, 1, 1, Rc = 6.0d0, Rs = 3.0d0, eta = 1.0d0)
    call sf_add_rad(3, 1, 1, Rc = 6.0d0, kappa = 1.0d0)
    call sf_add_ang(4, 1, 1, 1, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 1, 1, 1, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)

    nG = sf_nG_type(1)

    call sf_final()

    has_passed = tst_equal(nG, 5)
    call tst_check_passed(has_passed)

  end subroutine test_setup

  !--------------------------------------------------------------------!

  subroutine test_fcc_rotation()

    implicit none

    integer, parameter :: NMAX = 1000
    integer, parameter :: GMAX = 100

    integer                             :: ntypes, nG
    double precision                    :: Rc
    double precision                    :: alat
    double precision, dimension(3,3)    :: avec
    double precision, dimension(3,NMAX) :: X
    double precision, dimension(3)      :: Xi
    integer,          dimension(NMAX)   :: types
    double precision, dimension(GMAX)   :: G1, G2
    double precision, dimension(3,GMAX)      :: dGi1, dGi2
    double precision, dimension(3,GMAX,NMAX) :: dGj1, dGj2
    integer                             :: nx

    double precision :: ax
    double precision, dimension(3,3) :: R


    logical :: has_passed

    call tst_new("Symmfunc Test 2: FCC rotation")
    has_passed = .true.

    ntypes = 2
    nG = 1
    call sf_init(ntypes, nG)

    call sf_add_rad(1, 1, 1, Rc = 3.0d0)
    call sf_add_rad(2, 1, 1, Rc = 6.0d0, Rs = 3.0d0, eta = 1.0d0)
    call sf_add_rad(3, 1, 1, Rc = 6.0d0, kappa = 1.0d0)
    call sf_add_rad(1, 1, 2, Rc = 3.0d0)
    call sf_add_rad(2, 1, 2, Rc = 6.0d0, Rs = 3.0d0, eta = 1.0d0)
    call sf_add_rad(3, 1, 2, Rc = 6.0d0, kappa = 1.0d0)
    call sf_add_ang(4, 1, 1, 1, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 1, 1, 1, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(4, 1, 1, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 1, 1, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(4, 1, 2, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 1, 2, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)

    call sf_add_rad(1, 2, 1, Rc = 3.0d0)
    call sf_add_rad(2, 2, 1, Rc = 6.0d0, Rs = 3.0d0, eta = 1.0d0)
    call sf_add_rad(3, 2, 1, Rc = 6.0d0, kappa = 1.0d0)
    call sf_add_rad(1, 2, 2, Rc = 3.0d0)
    call sf_add_rad(2, 2, 2, Rc = 6.0d0, Rs = 3.0d0, eta = 1.0d0)
    call sf_add_rad(3, 2, 2, Rc = 6.0d0, kappa = 1.0d0)
    call sf_add_ang(4, 2, 1, 1, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 2, 1, 1, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(4, 2, 1, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 2, 1, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(4, 2, 2, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 2, 2, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)

    nG = sf_nG_type(1)

    ! (1) FCC normal orientation

    Rc = 6.0d0

    alat = 2.5d0
    avec(1:3,1) = (/ 0.0d0, 0.5d0, 0.5d0 /)*alat
    avec(1:3,2) = (/ 0.5d0, 0.0d0, 0.5d0 /)*alat
    avec(1:3,3) = (/ 0.5d0, 0.5d0, 0.0d0 /)*alat

    Xi(1:3) = (/ 0.0d0, 0.0d0, 0.0d0 /)
    nx = NMAX
    call get_coo(Rc, avec, Xi, nx, X)
    types(1:nx/2)    = 1
    types(nx/2+1:nx) = 2

    Xi = Xi + (/ 0.1d0, 0.0d0, 0.0d0 /)
    call sf_fingerprint(1, Xi, nx, X, types, GMAX, G1, dGi1, dGj1)

    ! (2) FCC rotated by 45 degrees around x axis and translated

    ax = sqrt(0.5d0)
    R(1:3,1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
    R(1:3,2) = (/ 0.0d0,    ax,    ax /)
    R(1:3,3) = (/ 0.0d0,   -ax,    ax /)

    avec = matmul(R, avec)

    ! translation:
    Xi(1:3) = (/ 0.2d0, 0.4d0, 0.6d0 /)
    nx = NMAX
    call get_coo(Rc, avec, Xi, nx, X)

    ! we can only look at the x-axis:
    Xi = Xi + (/ 0.1d0, 0.0d0, 0.0d0 /)
    call sf_fingerprint(1, Xi, nx, X, types, GMAX, G2, dGi2, dGj2)

    ! (3) check whether both structures gave the same results

    call sf_final()

    has_passed = (has_passed .and. tst_equal(G1, G2, prec=1.0d-6))
    ! derivatives: only x direction has to be equal, since that is the
    !              axis we rotated about
    has_passed = (has_passed .and. tst_equal(&
                  dGi1(1,1:nG), dGi2(1,1:nG), prec=1.0d-6))
    has_passed = (has_passed .and. tst_equal(&
                  dGj1(1,1:nG,1:nx), dGj2(1,1:nG,1:nx), prec=1.0d-6))

    call tst_check_passed(has_passed)

  end subroutine test_fcc_rotation

  !--------------- check symmetry function derivatives ----------------!

  subroutine test_derivatives()

    implicit none

    integer, parameter :: NMAX = 1000
    integer, parameter :: GMAX = 100

    integer                             :: ntypes, nG
    double precision                    :: Rc
    double precision                    :: alat
    double precision, dimension(3,3)    :: avec
    double precision, dimension(3,NMAX) :: X
    double precision, dimension(3)      :: Xi
    integer,          dimension(NMAX)   :: types
    double precision, dimension(GMAX)   :: G0, G1, G2
    double precision, dimension(3,GMAX)      :: dGi0, dGi1
    double precision, dimension(3,GMAX,NMAX) :: dGj0, dGj1
    integer                             :: nx

    double precision :: ax
    double precision, dimension(3,3) :: R

    double precision :: d
    double precision, dimension(3,3) :: dd
    integer :: i, j, iat

    logical :: has_passed

    call tst_new("Symmfunc Test 3: derivatives")
    has_passed = .true.

    ntypes = 2
    nG = 1
    call sf_init(ntypes, nG)

    call sf_add_rad(1, 1, 1, Rc = 3.0d0)
    call sf_add_rad(2, 1, 1, Rc = 6.0d0, Rs = 3.0d0, eta = 1.0d0)
    call sf_add_rad(3, 1, 1, Rc = 6.0d0, kappa = 1.0d0)
    call sf_add_rad(1, 1, 2, Rc = 3.0d0)
    call sf_add_rad(2, 1, 2, Rc = 6.0d0, Rs = 3.0d0, eta = 1.0d0)
    call sf_add_rad(3, 1, 2, Rc = 6.0d0, kappa = 1.0d0)
    call sf_add_ang(4, 1, 1, 1, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 1, 1, 1, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(4, 1, 1, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 1, 1, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(4, 1, 2, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 1, 2, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)

    call sf_add_rad(1, 2, 1, Rc = 3.0d0)
    call sf_add_rad(2, 2, 1, Rc = 6.0d0, Rs = 3.0d0, eta = 1.0d0)
    call sf_add_rad(3, 2, 1, Rc = 6.0d0, kappa = 1.0d0)
    call sf_add_rad(1, 2, 2, Rc = 3.0d0)
    call sf_add_rad(2, 2, 2, Rc = 6.0d0, Rs = 3.0d0, eta = 1.0d0)
    call sf_add_rad(3, 2, 2, Rc = 6.0d0, kappa = 1.0d0)
    call sf_add_ang(4, 2, 1, 1, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 2, 1, 1, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(4, 2, 1, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 2, 1, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(4, 2, 2, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)
    call sf_add_ang(5, 2, 2, 2, Rc = 6.0d0, eta = 1.0d0, lambda = 1.0d0, zeta = 1.0d0)

    nG = sf_nG_type(1)

    ! rotated FCC basis
    alat = 2.5d0
    avec(1:3,1) = (/ 0.0d0, 0.5d0, 0.5d0 /)*alat
    avec(1:3,2) = (/ 0.5d0, 0.0d0, 0.5d0 /)*alat
    avec(1:3,3) = (/ 0.5d0, 0.5d0, 0.0d0 /)*alat
    ax = sqrt(0.3d0)
    R(1:3,1) = (/ 1.0d0, 0.0d0, 0.0d0 /)
    R(1:3,2) = (/ 0.0d0,    ax,    ax /)
    R(1:3,3) = (/ 0.0d0,   -ax,    ax /)
    avec = matmul(R, avec)
    ax = sqrt(0.7d0)
    R(1:3,1) = (/    ax, 0.0d0,    ax /)
    R(1:3,2) = (/ 0.0d0, 1.0d0, 0.0d0 /)
    R(1:3,3) = (/   -ax, 0.0d0,    ax /)
    avec = matmul(R, avec)

    ! translation
    Xi(1:3) = (/ 0.2d0, 0.4d0, 0.6d0 /)
    nx = NMAX
    Rc = 6.0d0
    call get_coo(Rc, avec, Xi, nx, X)

    types(1:nx/2)    = 1
    types(nx/2+1:nx) = 2
    Xi = Xi + (/ 0.1d0, -0.3d0, 0.2d0 /)
    call sf_fingerprint(1, Xi, nx, X, types, GMAX, G0, dGi0, dGj0)

    ! numerical derivatives wrt. central atom
    d = 0.01d0
    dd(1,:) = [d, 0.0d0, 0.0d0]
    dd(2,:) = [0.0d0, d, 0.0d0]
    dd(3,:) = [0.0d0, 0.0d0, d]
    do i = 1, 3
       Xi = Xi - dd(i,:)
       call sf_fingerprint(1, Xi, nx, X, types, GMAX, G1)
       Xi = Xi + 2.0d0*dd(i,:)
       call sf_fingerprint(1, Xi, nx, X, types, GMAX, G2)
       Xi = Xi - dd(i,:)
       dGi1(i,1:nG) = (G2(1:nG) - G1(1:nG))/(2.0d0*d)
    end do

    has_passed = tst_equal(dGi0(:,1:nG), dGi1(:,1:nG), prec=0.01d0)
    if (.not. has_passed) then
       open(99, file='TEST_DERIVATIVES_I', status='replace', action='write')
       do i = 1, nG
          write(99, '(9(1x,ES15.8))') dGi0(:,i), dGi1(:,i), dGi0(:,i) - dGi1(:,i)
       end do
       close(99)
    end if

    ! numerical derivatives wrt. other atoms
    call random_init()
    check_j : do j = 1, 10
       iat = random_integer(nx)
       do i = 1, 3
          X(:,iat) = X(:,iat) - dd(i,:)
          call sf_fingerprint(1, Xi, nx, X, types, GMAX, G1)
          X(:,iat) = X(:,iat) + 2.0d0*dd(i,:)
          call sf_fingerprint(1, Xi, nx, X, types, GMAX, G2)
          X(:,iat) = X(:,iat) - dd(i,:)
          dGj1(i,1:nG,iat) = (G2(1:nG) - G1(1:nG))/(2.0d0*d)
       end do
       has_passed = tst_equal(&
            dGj0(:,1:nG,iat), dGj1(:,1:nG,iat), tol=1.0d-5, prec=0.01d0)
       if (.not. has_passed) then
          open(99, file='TEST_DERIVATIVES_J', status='replace', action='write')
          do i = 1, nG
             write(99, '(9(1x,ES15.8))') dGj0(:,i,iat), dGj1(:,i,iat), &
                                         dGj0(:,i,iat) - dGj1(:,i,iat)
          end do
          close(99)
          exit check_j
       end if
    end do check_j
    call random_final()

    call tst_check_passed(has_passed)

    call sf_final()

  end subroutine test_derivatives

  !------------ auxiliary routine to generate coordinates -------------!

  subroutine get_coo(Rc, avec, X0, nx, X)

    implicit none

    double precision,                  intent(in)    :: Rc
    double precision, dimension(3,3),  intent(in)    :: avec
    double precision, dimension(3),    intent(inout) :: X0
    integer,                           intent(inout) :: nx
    double precision, dimension(3,nx), intent(inout) :: X

    integer                        :: ix, iy, iz
    double precision               :: d2, Rc2
    double precision, dimension(3) :: t
    integer                        :: ipt
    logical                        :: nextx

    Rc2 = Rc*Rc
    ipt = 0

    ! convert X0 to cartesian coordinates:
    X0(:) = X0(1)*avec(:,1) + X0(2)*avec(:,2) + X0(3)*avec(:,3)

    iz = 0
    zloop : do
       iy = 0
       yloop : do
          ix = 0
          xloop : do

             nextx = .false.
             if (all((/ix,iy,iz/) == 0)) then
                ix = ix + 1
                cycle xloop
             end if

             t(:) = dble((/ix, iy, iz/))
             t(:) = t(1)*avec(:,1) + t(2)*avec(:,2) + t(3)*avec(:,3)
             d2 = sum(t*t)
             if (d2 < Rc2) then
                ipt = ipt + 1
                X(:,ipt) = X0(:) + t(:)
                if (iz /= 0) then
                   ipt = ipt + 1
                   X(:,ipt) = X0(:) - t(:)
                end if
                nextx = .true.
             end if

             if (ix /= 0)  then

                t(:) = dble((/-ix, iy, iz/))
                t(:) = t(1)*avec(:,1) + t(2)*avec(:,2) + t(3)*avec(:,3)
                d2 = sum(t*t)
                if (d2 < Rc2) then
                   ipt = ipt + 1
                   X(:,ipt) = X0(:) + t(:)
                   if (iz /= 0) then
                      ipt = ipt + 1
                      X(:,ipt) = X0(:) - t(:)
                   end if
                   nextx = .true.
                end if

                if (iy /= 0) then
                   t(:) = dble((/ix, -iy, iz/))
                   t(:) = t(1)*avec(:,1) + t(2)*avec(:,2) + t(3)*avec(:,3)
                   d2 = sum(t*t)
                   if (d2 < Rc2) then
                      ipt = ipt + 1
                      X(:,ipt) = X0(:) + t(:)
                      if (iz /= 0) then
                         ipt = ipt + 1
                         X(:,ipt) = X0(:) - t(:)
                      end if
                      nextx = .true.
                   end if
                   t(:) = dble((/-ix, -iy, iz/))
                   t(:) = t(1)*avec(:,1) + t(2)*avec(:,2) + t(3)*avec(:,3)
                   d2 = sum(t*t)
                   if (d2 < Rc2) then
                      ipt = ipt + 1
                      X(:,ipt) = X0(:) + t(:)
                      if (iz /= 0) then
                         ipt = ipt + 1
                         X(:,ipt) = X0(:) - t(:)
                      end if
                      nextx = .true.
                   end if
                end if ! iy /= 0

             else ! ix == 0

                if (iy /= 0) then
                   t(:) = dble((/ix, -iy, iz/))
                   t(:) = t(1)*avec(:,1) + t(2)*avec(:,2) + t(3)*avec(:,3)
                   d2 = sum(t*t)
                   if (d2 < Rc2) then
                      ipt = ipt + 1
                      X(:,ipt) = X0(:) + t(:)
                      if (iz /= 0) then
                         ipt = ipt + 1
                         X(:,ipt) = X0(:) - t(:)
                      end if
                      nextx = .true.
                   end if
                end if ! iy /= 0

             end if ! ix /= 0

             if (.not. nextx) exit xloop
             ix = ix + 1
          end do xloop
          if (ix == 0) exit yloop
          iy = iy + 1
       end do yloop
       if (iy == 0) exit zloop
       iz = iz + 1
    end do zloop

    nx = ipt

  end subroutine get_coo

end program test_symmfunc
