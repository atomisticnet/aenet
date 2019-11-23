!-----------------------------------------------------------------------
!                     Compute structural fingerprint
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
! 2016-06-02 Alexander Urban (AU) and Nongnuch Artrith (NA)
!-----------------------------------------------------------------------
program fingerprint

  use geometry, only: geo_init,     &
                      geo_final,    &
                      geo_write,    &
                      pbc,          &
                      latticeVec,   &
                      nAtoms,       &
                      atomType,     &
                      atomTypeName, &
                      cooLatt, cooCart

  use io,       only: io_adjustl

  use lclist,   only: lcl_init,        &
                      lcl_final,       &
                      lcl_nmax_nbdist, &
                      lcl_nbdist_cart

  use sfbasis,  only: FingerprintBasis, &
                      new_SFBasis,      &
                      del_SFBasis,      &
                      sfb_eval

  implicit none

  double precision, parameter :: RMIN = 0.1d0

  character(len=1024) :: filename, frmt
  integer             :: order1, order2
  double precision    :: rcut1, rcut2, rcut
  integer             :: ntypes
  integer             :: nnb, nnb_max
  integer             :: i, iatom
  type(FingerprintBasis) :: sfb

  character(len=2), dimension(:),   allocatable :: atom_types
  double precision, dimension(:,:), allocatable :: nbcoo
  double precision, dimension(:),   allocatable :: nbdist
  integer,          dimension(:),   allocatable :: nblist
  integer,          dimension(:),   allocatable :: nbtype
  double precision, dimension(:),   allocatable :: values, avg_values

  call initialize()
  sfb = new_SFBasis(ntypes, atom_types, order1, order2, rcut1, rcut2)
  call geo_init(filename, 'xsf')
  rcut = max(rcut1, rcut2)
  call lcl_init(RMIN, rcut, latticeVec, nAtoms, atomType, cooLatt, pbc)

  nnb_max = lcl_nmax_nbdist(RMIN, rcut)
  allocate(nbcoo(3,nnb_max), &
           nbdist(nnb_max),  &
           nblist(nnb_max),  &
           nbtype(nnb_max),  &
           values(sfb%N),    &
           avg_values(sfb%N))

  avg_values(:) = 0.0d0
  do iatom = 1, nAtoms
     nnb = nnb_max
     call lcl_nbdist_cart(iatom, nnb, nbcoo, nbdist, nblist=nblist, &
                          nbtype=nbtype)
     write(0,*) trim(filename), nnb
     call sfb_eval(sfb, atomType(iatom), cooCart(iatom), nnb, &
                   nbtype(1:nnb), nbcoo(:, 1:nnb), sfb%N, values)
     avg_values(1:sfb%N) = avg_values(1:sfb%N) + values(1:sfb%N)
  end do

  frmt = '(A,1x,' // trim(io_adjustl(sfb%N)) // '(ES15.8,1x))'
  write(*, frmt) trim(filename), avg_values(1:sfb%N)/dble(nAtoms)

  deallocate(nbcoo, nbdist, nblist, nbtype, values, avg_values)

  call lcl_final()
  call geo_final()
  call del_SFBasis(sfb)
  call finalize()

contains

  subroutine initialize()

    implicit none

    character(len=1024) :: str
    integer             :: i, nargs

    nargs = command_argument_count()
    if (nargs < 6) then
       write(0,*) "Error: At least 6 arguments expected."
       call print_usage()
       call finalize()
       stop
    end if

    call get_command_argument(1, value=filename)
    call get_command_argument(2, value=str)
    read(str, *) rcut1
    call get_command_argument(3, value=str)
    read(str, *) order1
    call get_command_argument(4, value=str)
    read(str, *) rcut2
    call get_command_argument(5, value=str)
    read(str, *) order2

    ntypes = nargs - 5
    allocate(atom_types(ntypes))
    do i = 6, nargs
       call get_command_argument(i, atom_types(i-5))
    end do

  end subroutine initialize

  !--------------------------------------------------------------------!

  subroutine finalize()

    implicit none

    if (allocated(atom_types)) deallocate(atom_types)

  end subroutine finalize

  !--------------------------------------------------------------------!

  subroutine print_usage()

    implicit none

    write(*,*) "Usage: fingerprint.x structure_file rcut1 N1 rcut2 N2 <species>"
    write(*,*)
    write(*,*) "  structure_file   Path to an atomic structure file in XSF format."
    write(*,*) "  rcut1 N1         Radial cutoff and expansion order of the radial basis."
    write(*,*) "  rcut2 N2         Radial cutoff and expansion order of the angular basis."
    write(*,*) "  <species>        List of chemical symbols."
    write(*,*)

  end subroutine print_usage

end program fingerprint
