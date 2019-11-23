!-----------------------------------------------------------------------
! xsflib.f90 -  A library with I/O routines to parse the XCrysDen
!               Structure File format (XSF).
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
! Currently, the XSF format is supported by (at least) the following
! free molecular graphics tools: XCrysDen, JMol, Vesta, VMD
!-----------------------------------------------------------------------
! 2011-02-11 Alexander Urban (AU), Nongnuch Artrith (NA)
! 2011-11-07 -- changed names and syntax to conform with other
!               packages, e.g. cooCart instead of coorat
! 2013-05-23 -- (1) allow unsorted atom types; atoms of the same
!               do no longer need to be consecutive in the file
!               (2) use io_unit() instead of fix file unit number
! 2013-08-24 -- support for isolated structures
!-----------------------------------------------------------------------

module xsflib

  use io, only: io_adjustl,  &
                io_isfigure, &
                io_isletter, &
                io_readnext, &
                io_readval,  &
                io_split,    &
                io_lower,    &
                io_unit,     &
                io_upper

  implicit none
  private

  public :: xsf_init,  &
            xsf_final


  !--------------------------------------------------------------------!
  ! LENLINE : max. length of lines (number of char's) in input file    !
  !--------------------------------------------------------------------!

  integer, parameter, private :: LENLINE = 1024

  !---------------------------- structure -----------------------------!
  ! nAtoms             : total number of atoms                         !
  ! nTypes             : number of atom types / atomic species         !
  ! unitCellVol        : unit cell volume                              !
  ! latticeVec(i,j)    : i-th component of the j-th lattice vector     !
  ! atomTypeName(ityp) : name of atomic species ityp; ityp = 1, nTypes !
  ! nAtomsOfType(ityp) : number of atoms of atomic species ityp        !
  ! cooCart(i,iat)     : i-th component of the cartesian coordinates   !
  !                      of atom number iat                            !
  ! forCart(i,iat)     : same for the cartesian forces                 !
  !                                                                    !
  !------------------------------ output ------------------------------!
  ! hasOutput          : .true., if file contains output information   !
  ! hasForces          : .true., if forces are present                 !
  ! E_tot              : total energy                                  !
  ! E_coh              : cohesive energy                               !
  !--------------------------------------------------------------------!

  logical,                                       public :: pbc
  integer,                                       public :: nAtoms
  integer,                                       public :: nTypes
  double precision,                              public :: unitCellVol
  double precision, dimension(3,3),              public :: latticeVec
  character(len=2), dimension(:),   allocatable, public :: atomTypeName
  integer,          dimension(:),   allocatable, public :: nAtomsOfType
  integer,          dimension(:),   allocatable, public :: atomType
  double precision, dimension(:,:), allocatable, public :: cooCart
  double precision, dimension(:,:), allocatable, public :: forCart

  logical,                                       public :: hasOutput
  logical,                                       public :: hasForces
  double precision,                              public :: E_tot
  double precision,                              public :: E_coh

  !--------------------------------------------------------------------!
  !                              PRIVATE                               !
  !                                                                    !
  ! iline   : number of the current line of the input file             !
  ! line    : input buffer for single line from the input file         !
  !--------------------------------------------------------------------!

  integer,                private :: iline
  character(len=LENLINE), private :: line
  integer,                private :: u_in

contains

  subroutine xsf_init(file)

    implicit none

    character(len=*), intent(in) :: file

    logical                      :: fexists
    integer                      :: status

    character(len=LENLINE)       :: species

    inquire(file=trim(file), exist=fexists)
    if (.not. fexists) then
       write(0,*) 'Error: File not found: ', trim(file)
       stop
    end if

    ! initial values:
    hasOutput  = .false.
    hasForces  = .false.
    pbc        = .false.
    E_coh      = 0.0d0
    E_tot      = 0.0d0

    u_in = io_unit()
    open(u_in, file=trim(file), status='old', action='read')

    iline = 1
    read_input : do

       read(u_in, '(A)', iostat=status) line
       if (status /= 0) exit read_input
       iline = iline + 1

       ! skip blank lines
       if (len_trim(line) == 0) cycle read_input

       ! cohesive energy
       if (index(line, 'cohesive energy') > 0) then
          call io_readval(line, name='cohesive energy', val=E_coh)
          hasOutput = .true.
          cycle read_input
       end if

       ! total energy
       if (index(line, 'total energy') > 0) then
          call io_readval(line, name='total energy', val=E_tot)
          hasOutput = .true.
          cycle read_input
       end if

       ! skip all remaining comments
       line = trim(adjustl(line))
       if (line(1:1) == '#') cycle read_input

       ! lattice vectors
       if (index(line, 'PRIMVEC') > 0) then
          pbc = .true.
          read(u_in, *) latticeVec(1:3, 1)
          read(u_in, *) latticeVec(1:3, 2)
          read(u_in, *) latticeVec(1:3, 3)
          iline = iline + 3
          cycle read_input
       end if

       ! coordinates of periodic structures
       primcoord : if (index(line, 'PRIMCOORD') > 0) then
          read(u_in, *) nAtoms
          iline = iline + 1
          call xsf_count_atoms_and_types(nAtoms, nTypes, species, hasForces)
          call xsf_malloc()
          call io_split(species, atomTypeName)
          call xsf_read_coordinates()
          cycle read_input
       end if primcoord

       ! coordinates of isolated structures
       atoms : if (index(line, 'ATOMS') > 0) then
          call xsf_count_atoms_and_types(nAtoms, nTypes, species, hasForces)
          call xsf_malloc()
          call io_split(species, atomTypeName)
          call xsf_read_coordinates()
          cycle read_input
       end if atoms

    end do read_input

    close(u_in)

    if (pbc) then
       ! unit cell volume:
       unitCellVol = latticeVec(1,1)*latticeVec(2,2)*latticeVec(3,3) &
                   + latticeVec(2,1)*latticeVec(3,2)*latticeVec(1,3) &
                   + latticeVec(3,1)*latticeVec(1,2)*latticeVec(2,3) &
                   - latticeVec(3,1)*latticeVec(2,2)*latticeVec(1,3) &
                   - latticeVec(1,1)*latticeVec(3,2)*latticeVec(2,3) &
                   - latticeVec(2,1)*latticeVec(1,2)*latticeVec(3,3)
    end if

  end subroutine xsf_init

  !--------------------------------------------------------------------!

  subroutine xsf_final()

    implicit none

    if (allocated(cooCart)) then
       deallocate(cooCart, atomType, atomTypeName, nAtomsOfType)
    end if
    if (allocated(forCart)) deallocate(forCart)

  end subroutine xsf_final



  !====================================================================!
  !                                                                    !
  !                        auxiliary procedures                        !
  !                                                                    !
  !====================================================================!



  !--------------------------------------------------------------------!
  !                          allocate memory                           !
  !--------------------------------------------------------------------!

  subroutine xsf_malloc()

    implicit none

    !-----------------------  module variables ------------------------!
    ! nAtoms, nTypes, hasForces                                        !
    !------------------------------------------------------------------!

    if (allocated(cooCart)) then
       write(0,*) 'Warning: attempted double memory allocation in XSF module.'
       return
    end if

    if (hasForces) then
       allocate(cooCart(3,nAtoms), nAtomsOfType(nTypes), &
                atomTypeName(nTypes), atomType(nAtoms),  &
                forCart(3,nAtoms))
    else
       allocate(cooCart(3,nAtoms), nAtomsOfType(nTypes), &
                atomTypeName(nTypes), atomType(nAtoms))
    end if

  end subroutine xsf_malloc

  !--------------------------------------------------------------------!
  !              atom type index for given atom type name              !
  !--------------------------------------------------------------------!

  function xsf_itype(el) result(itype)

    implicit none

    !------------------------ module variables ------------------------!
    ! nTypes, atomTypeName                                             !
    !------------------------------------------------------------------!

    character(len=*), intent(in) :: el
    integer                      :: itype

    itype = 1
    do while(itype <= nTypes)
       if (trim(adjustl(el)) == trim(adjustl(atomTypeName(itype)))) then
          exit
       end if
       itype = itype + 1
    end do
    if (itype > nTypes) itype = -1

  end function xsf_itype

  !--------------------------------------------------------------------!
  !                   count atoms and atomic species                   !
  !                                                                    !
  ! Assumes that the next line contains the species, coordinates, and  !
  ! possibly forces of the FIRST atom.                                 !
  !                                                                    !
  ! The string 'species' is used to store all chemical symbols in one  !
  ! pass.  512 characters are more than enough to store all chemical   !
  ! symbols of the PSE separated by a blank.                           !
  !--------------------------------------------------------------------!

  subroutine xsf_count_atoms_and_types(natoms, ntypes, species, hasforces)

    implicit none

    !------------------------------------------------------------------!
    ! natoms     number of atoms                                       !
    ! ntypes     number of atomic species                              !
    ! species    character string continaing all atomic species        !
    ! hasforces  .true., iff atomic forces are present                 !
    !------------------------------------------------------------------!

    integer,          intent(out) :: natoms
    integer,          intent(out) :: ntypes
    character(len=*), intent(out) :: species
    logical,          intent(out) :: hasforces

    !------------------------ module variables ------------------------!
    ! u_in, iline, line                                                !
    !------------------------------------------------------------------!

    character(len=2)       :: el
    character(len=LENLINE) :: str
    integer                :: stat
    integer                :: i, k

    ntypes    = 0
    natoms    = 0
    species   = ' '
    hasforces = .false.

    ! count atoms and types
    count : do
       read(u_in, '(A)', iostat=stat) str
       if (stat /= 0) exit count
       if (len_trim(str) == 0) exit count
       line = str
       read(line, *) el
       ! blanks are used to distinguish, e.g., 'S' from 'Sn'
       if (index(species,' ' // trim(el) // ' ' ) == 0) then
          species = trim(species) // ' ' // trim(el) // ' '
          ntypes   = ntypes + 1
       end if
       natoms = natoms + 1
    end do count

    ! are forces present?
    ! --> count records in final atom line
    k = 0
    i = 1
    do while (i > 0)
       k = k + 1
       call io_readnext(line, i, str)
    end do
    if (k >= 7) hasForces = .true.

    ! rewind to original position in file
    rewind(u_in)
    do i = 1, iline-1
       read(u_in, *)
    end do

  end subroutine xsf_count_atoms_and_types

  !--------------------------------------------------------------------!
  !                    read coordinates and forces                     !
  !--------------------------------------------------------------------!

  subroutine xsf_read_coordinates()

    implicit none

    !------------------------ module variables ------------------------!
    ! line, iline, nAtoms, nTypes, atomType, nAtomsOfType, cooCart,    !
    ! forCart                                                          !
    !------------------------------------------------------------------!

    integer          :: iat, itype
    character(len=2) :: el

    ! read Cartesian coordinates and forces
    do iat = 1, nAtoms
       read(u_in, '(A)') line
       iline = iline + 1
       if (hasForces) then
          read(line, *) el, cooCart(1:3,iat), forCart(1:3,iat)
       else
          read(line, *) el, cooCart(1:3,iat)
       end if
       itype = xsf_itype(el)
       atomType(iat) = itype
       nAtomsOfType(itype) = nAtomsOfType(itype) + 1
    end do

  end subroutine xsf_read_coordinates

end module xsflib
