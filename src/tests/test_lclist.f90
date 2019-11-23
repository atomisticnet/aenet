!-----------------------------------------------------------------------
!            untest_lclist.f90 - Unit tests for lclist.f90
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
! 2012-06-04 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

program test_lclist

  use lclist, only: lcl_init,          &
                    lcl_final,         &
                    lcl_print_info,    &
                    lcl_nmax_cell,     &
                    lcl_nmax_nblist,   &
                    lcl_nmax_nbdist,   &
                    lcl_nblist,        &
                    lcl_nbdist,        &
                    lcl_nbdist_cart

  use unittest, only: tst_new, tst_check_passed

  implicit none

  call test_fcc_primitive()
  call test_fcc_conventional()
  call test_fcc_isolated()

contains

  !--------------------------------------------------------------------!
  !             Test 1 - FCC lattice (primitive unit cell)             !
  !--------------------------------------------------------------------!

  subroutine test_fcc_primitive()

    implicit none

    double precision                              :: a, d_NN
    double precision, dimension(3,3)              :: avec
    integer                                       :: natoms
    double precision, dimension(:,:), allocatable :: coords
    integer,          dimension(:),   allocatable :: types
    double precision                              :: r_min, r_max
    integer                                       :: nmax
    double precision, dimension(:,:), allocatable :: nbcoords
    double precision, dimension(:),   allocatable :: nbdist
    integer                                       :: stat

    logical :: has_passed

    call tst_new("LClist Test 1 - FCC primitive unit cell")
    has_passed = .true.

    a = 5.0d0
    d_NN = 0.5d0*sqrt(2.0d0)*a

    avec(:,1) = (/ 0.0d0, 0.5d0, 0.5d0  /)*a
    avec(:,2) = (/ 0.5d0, 0.0d0, 0.5d0  /)*a
    avec(:,3) = (/ 0.5d0, 0.5d0, 0.0d0  /)*a

    natoms = 1
    allocate(coords(3,natoms), types(natoms))
    coords(:,1) = (/ 0.0d0, 0.0d0, 0.0d0 /)
    types(1)    = 1

    r_min = 1.0d0
    r_max = d_NN + 0.1d0

    call lcl_init(r_min, r_max, avec, natoms, types(1:natoms), coords(:,1:natoms))
    nmax = lcl_nmax_nbdist(r_min, r_max)
    allocate(nbcoords(3,nmax), nbdist(nmax))
    call lcl_nbdist(1, nmax, nbcoords, nbdist, stat=stat)

    has_passed = (has_passed .and. (stat == 0))
    has_passed = (has_passed .and. (nmax == 12))

    call tst_check_passed(has_passed)
    if (.not. has_passed) then
       write(*,*) "stat = ", stat
       call lcl_print_info()
    end if

    deallocate(nbcoords, nbdist, coords, types)
    call lcl_final()

  end subroutine test_fcc_primitive

  !--------------------------------------------------------------------!
  !            Test 2 - FCC conventional cell 4x4x4 supercell          !
  !--------------------------------------------------------------------!

  subroutine test_fcc_conventional()

    implicit none

    double precision                              :: a, d_NN
    double precision, dimension(3,3)              :: avec
    integer                                       :: natoms
    double precision, dimension(:,:), allocatable :: coords, coords2
    integer,          dimension(:),   allocatable :: types
    integer                                       :: ix, iy, iz, i, iatom
    double precision                              :: r_min, r_max
    integer                                       :: nmax
    double precision, dimension(:,:), allocatable :: nbcoords
    double precision, dimension(:),   allocatable :: nbdist
    integer                                       :: stat

    logical :: has_passed

    call tst_new("LClist Test 2 - FCC conventional unit cell")
    has_passed = .true.

    a = 5.0d0
    d_NN = 0.5d0*sqrt(2.0d0)*a

    avec(:,1) = (/ 1.0d0, 0.0d0, 0.0d0  /)*a
    avec(:,2) = (/ 0.0d0, 1.0d0, 0.0d0  /)*a
    avec(:,3) = (/ 0.0d0, 0.0d0, 1.0d0  /)*a
    avec      = 4.0d0*avec

    natoms = 4 * 4*4*4
    allocate(coords2(3,4), coords(3,natoms), types(natoms))
    coords2(:,1) = (/ 0.0d0, 0.0d0, 0.0d0 /)
    coords2(:,2) = (/ 0.0d0, 0.5d0, 0.5d0 /)
    coords2(:,3) = (/ 0.5d0, 0.0d0, 0.5d0 /)
    coords2(:,4) = (/ 0.5d0, 0.5d0, 0.0d0 /)
    types(1:natoms) = 1

    iatom = 0
    do ix = 0, 3
       do iy = 0, 3
          do iz = 0, 3
             do i = 1, 4
                iatom = iatom + 1
                coords(:,iatom) = 0.25d0*(coords2(:,i) &
                     + (/ dble(ix), dble(iy), dble(iz) /))
             end do
          end do
       end do
    end do
    deallocate(coords2)

    r_min = 1.0d0
    r_max = d_NN + 0.1d0

    call lcl_init(r_min, r_max, avec, natoms, types(1:natoms), coords(:,1:natoms))
    nmax = lcl_nmax_nbdist(r_min, r_max)
    allocate(nbcoords(3,nmax), nbdist(nmax))
    atoms : do iatom = 1, natoms
       call lcl_nbdist(iatom, nmax, nbcoords, nbdist, stat=stat)
       has_passed = (has_passed .and. (stat == 0))
       has_passed = (has_passed .and. (nmax == 12))
    end do atoms

    deallocate(nbcoords, nbdist, coords, types)
    call lcl_final()

    call tst_check_passed(has_passed)

  end subroutine test_fcc_conventional

  !--------------------------------------------------------------------!
  !            Test 3 - isolated structure (2x2x2 FCC cube)            !
  !--------------------------------------------------------------------!

  subroutine test_fcc_isolated()

    implicit none

    double precision                              :: a, d_NN
    double precision, dimension(3,3)              :: avec
    integer                                       :: natoms
    double precision, dimension(:,:), allocatable :: coords, coords2
    integer,          dimension(:),   allocatable :: types
    integer                                       :: ix, iy, iz, i, iatom
    double precision                              :: r_min, r_max
    integer                                       :: nmax, nmax_in
    double precision, dimension(:,:), allocatable :: nbcoords
    double precision, dimension(:),   allocatable :: nbdist
    integer                                       :: stat

    logical :: has_passed

    call tst_new("LClist Test 3 - isolated structure")
    has_passed = .true.

    a = 5.0d0
    d_NN = 0.5d0*sqrt(2.0d0)*a

    avec(:,1) = (/ 1.0d0, 0.0d0, 0.0d0  /)*a
    avec(:,2) = (/ 0.0d0, 1.0d0, 0.0d0  /)*a
    avec(:,3) = (/ 0.0d0, 0.0d0, 1.0d0  /)*a
    avec      = 2.0d0*avec

    natoms = 4 * 2*2*2
    allocate(coords2(3,4), coords(3,natoms), types(natoms))
    coords2(:,1) = (/ 0.0d0, 0.0d0, 0.0d0 /)
    coords2(:,2) = (/ 0.0d0, 0.5d0, 0.5d0 /)
    coords2(:,3) = (/ 0.5d0, 0.0d0, 0.5d0 /)
    coords2(:,4) = (/ 0.5d0, 0.5d0, 0.0d0 /)
    types(1:natoms) = 1

    iatom = 0
    do ix = 0, 1
       do iy = 0, 1
          do iz = 0, 1
             do i = 1, 4
                iatom = iatom + 1
                coords(:,iatom) = 0.5d0*(coords2(:,i) &
                     + (/ dble(ix), dble(iy), dble(iz) /))
             end do
          end do
       end do
    end do
    deallocate(coords2)

    r_min = 1.0d0
    r_max = d_NN + 0.1d0

    call lcl_init(r_min, r_max, avec, natoms, types(1:natoms), &
         coords(:,1:natoms), pbc_in=.false.)
    nmax_in = lcl_nmax_nbdist(r_min, r_max)
    allocate(nbcoords(3,nmax_in), nbdist(nmax_in))
    atoms2 : do iatom = 1, natoms
       nmax = nmax_in
       call lcl_nbdist(iatom, nmax, nbcoords, nbdist, stat=stat)
       has_passed = (has_passed .and. (stat == 0))
       has_passed = (has_passed .and. (nmax <= 12))
    end do atoms2
    call lcl_final()
    deallocate(nbcoords, nbdist, coords, types)

    call tst_check_passed(has_passed)

  end subroutine test_fcc_isolated

end program test_lclist
