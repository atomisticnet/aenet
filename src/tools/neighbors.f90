!-----------------------------------------------------------------------
!      Get all atoms within a specified distance of another atom
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
! 2016-05-15 Alexander Urban (AU) and Nongnuch Artrith (NA)
!-----------------------------------------------------------------------
program neighbors

  use geometry, only: geo_init,     &
                      geo_final,    &
                      geo_write,    &
                      pbc,          &
                      latticeVec,   &
                      nAtoms,       &
                      atomType,     &
                      atomTypeName, &
                      cooLatt

  use io,       only: io_adjustl

  use lclist,   only: lcl_init,        &
                      lcl_final,       &
                      lcl_nmax_nbdist, &
                      lcl_nbdist_cart

  implicit none

  double precision, parameter :: RMIN = 0.1d0

  character(len=1024) :: filename, frmt
  integer             :: atomid
  double precision    :: rcut
  integer             :: nnb, nnb_max
  integer             :: i

  double precision, dimension(:,:), allocatable :: nbcoo
  double precision, dimension(:),   allocatable :: nbdist
  integer,          dimension(:),   allocatable :: nblist
  integer,          dimension(:),   allocatable :: nbtype

  call initialize(filename, atomid, rcut)
  call geo_init(filename, 'xsf')
  call lcl_init(RMIN, rcut, latticeVec, nAtoms, atomType, cooLatt, pbc)

  nnb_max = lcl_nmax_nbdist(RMIN, rcut)
  allocate(nbcoo(3,nnb_max), nbdist(nnb_max), nblist(nnb_max), nbtype(nnb_max))

  nnb = nnb_max
  call lcl_nbdist_cart(atomid, nnb, nbcoo, nbdist, nblist=nblist, nbtype=nbtype)

  write(*,*) nnb
  write(*,*) "Neighbors of atom " // trim(io_adjustl(atomid)) // &
             " within " // "a cutoff of " // trim(io_adjustl(rcut)) // " A"
  frmt = "(1x,A2,3(2x,F15.8))"
  do i = 1, nnb
     write(*, frmt) atomTypeName(nbtype(i)), nbcoo(1,i), nbcoo(2,i), nbcoo(3,i)
  end do

  deallocate(nbcoo, nbdist, nblist, nbtype)

  call lcl_final()
  call geo_final()
  call finalize()

contains

  subroutine initialize(filename, atomid, rcut)

    implicit none

    character(len=*), intent(out) :: filename
    integer,          intent(out) :: atomid
    double precision, intent(out) :: rcut

    character(len=1024) :: str
    integer             :: nargs

    nargs = command_argument_count()
    if (nargs < 3) then
       write(0,*) "Error: Three arguments expected."
       call print_usage()
       call finalize()
       stop
    end if

    call get_command_argument(1, value=filename)
    call get_command_argument(2, value=str)
    read(str, *) atomid
    call get_command_argument(3, value=str)
    read(str, *) rcut

  end subroutine initialize

  !--------------------------------------------------------------------!

  subroutine finalize()

    implicit none

  end subroutine finalize

  !--------------------------------------------------------------------!

  subroutine print_usage()

    implicit none

    write(*,*) "Usage: neighbors.x structure_file atom_id cutoff"
    write(*,*)
    write(*,*) "  structure_file   Path to an atomic structure file in XSF format."
    write(*,*) "  atom_id          Number of the centeral atom starting with 1."
    write(*,*) "  cutoff           Cutoff radius in Angstrom."
    write(*,*)

  end subroutine print_usage

end program neighbors
