!-----------------------------------------------------------------------
!         geometry.f90 - interface to atomic structure formats
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
! 2011-10-20 Alexander Urban (AU), Nongnuch Artrith (NA)
! 2013-08-25 AU -- support for isolated structures
!-----------------------------------------------------------------------

module geometry

  use aeio,      only: TYPELEN, PATHLEN, LINELEN, STDOUT
  use constants, only: PI
  use io,        only: io_adjustl, &
                       io_lower,   &
                       io_unit

  implicit none
  save

  public  :: geo_init,          &
             geo_final,         &
             geo_cell_volume,   &
             geo_recip_lattice, &
             geo_itype_of_name, &
             geo_type_conv,     &
             geo_update_bounds, &
             cooCart,           &
             forCart

  !-------------------------- lattice basis ---------------------------!
  ! pbc                .true. for periodic structures, .false. else    !
  ! latticeVec(i,j)    i-th component of the j-th lattice vector       !
  ! recLattVec(i,j)    i-th component of the j-th recip. latt. vector  !
  ! origin(i)          i-th Cartesian component of the origin of the   !
  !                    basis; only relevant for isolated structures    !
  !---------------------------- structure -----------------------------!
  ! nAtoms             number of atoms in the structure                !
  ! nTypes             number of different atomic species in structure !
  ! atomTypeName(i)    name (symbol) of atom type ID i                 !
  ! atomType(i)        atom type ID of atom i                          !
  ! cooLatt(i,j)       i-th component of the j-th atomic coordinates   !
  !                    in the lattice basis (latticeVec); in case of   !
  !                    isolated structures the coordinates are scaled  !
  !                    to [0,1) and latticeVec contains the scaling.   !
  !------------------------------ output ------------------------------!
  ! hasForces          .true. iff the atomic forces are available and  !
  !                    memory for the forces has been allocated        !
  ! forLatt(i,j)       i-th component of the atomic forces acting on   !
  !                    the j-th atom in the lattice basis              !
  ! hasEnergy          .true. iff the energy of the structure is known !
  ! cohesiveEnergy     cohesive energy (eV)                            !
  ! total energy       total energy (eV) relative to reference method  !
  !---------------------------- input file ----------------------------!
  ! structureFile      path to the input structure file                !
  ! structureFormat    name of the structure format                    !
  !--------------------------------------------------------------------!

  logical,                                             public :: pbc
  double precision,       dimension(3,3),              public :: latticeVec
  double precision,       dimension(3,3),              public :: recLattVec
  double precision,       dimension(3),                public :: origin

  integer,                                             public :: nAtoms
  integer,                                             public :: nTypes
  character(len=TYPELEN), dimension(:),   allocatable, public :: atomTypeName
  integer,                dimension(:),   allocatable, public :: atomType
  double precision,       dimension(:,:), allocatable, public :: cooLatt

  logical,                                             public :: hasForces
  double precision,       dimension(:,:), allocatable, public :: forLatt
  logical,                                             public :: hasEnergy
  double precision,                                    public :: cohesiveEnergy
  double precision,                                    public :: totalEnergy

  character(len=PATHLEN),                              public :: structureFile
  character(len=LINELEN),                              public :: structureFormat

  !--------------------------------------------------------------------!

  logical, private :: isInit

contains

  subroutine geo_init(file, form)

    implicit none

    character(len=*), intent(in) :: file
    character(len=*), intent(in) :: form

    logical                      :: fexists

    if (isInit) then
       write(0,*) "Error: module already initialized in `geo_init'."
    end if

    inquire(file=file, exist=fexists)
    if (.not. fexists) then
       write(0,*) "Error: file not found in `geo_init': ", trim(file)
       stop
    end if

    ! defaults for output values
    hasForces = .false.
    hasEnergy = .false.

    select case(trim(io_lower(form)))
    case('xsf')
       call geo_init_xsf(file)
    case default
       write(0,*) "Error: unknown file format in `geo_init': ", trim(form)
       stop
    end select

    structureFile   = trim(file)
    structureFormat = trim(form)

    isInit = .true.

  end subroutine geo_init

  !--------------------------------------------------------------------!

  subroutine geo_final()

    implicit none

    if (.not. isInit) return

    if (allocated(cooLatt)) deallocate(cooLatt, atomType, atomTypeName)
    if (allocated(forlatt)) deallocate(forLatt)

    isInit = .false.

  end subroutine geo_final

  !--------------------------------------------------------------------!
  !            information about currently loaded structure            !
  !--------------------------------------------------------------------!

  subroutine geo_print_info(verbose)

    implicit none

    integer, optional, intent(in) :: verbose

    integer                :: v, iat
    character(len=LINELEN) :: frmt

    if (.not. isInit) then
       write(*,*) 'No structure loaded.'
       return
    end if

    ! verbosity level
    if (present(verbose)) then
       v = min(0,verbose)
    else
       v = 0
    end if

    write(*,*) "input file      : ", trim(structureFile)
    write(*,*) "input format    : ", trim(structureFormat)
    write(*,*) "number of atoms : ", trim(io_adjustl(nAtoms))
    write(*,*) "atomic species  : ", trim(io_adjustl(nTypes))
    write(*,*) "periodic        : ", pbc
    if (hasEnergy) then
       write(*,*) "cohesive energy : ", cohesiveEnergy, ' eV'
       write(*,*) "total energy    : ", totalEnergy, ' eV'
    end if
    write(*,*)

    if (v > 0) then
       write(*,*) "Lattice Vectors"
       write(*,*)
       write(*,'(1x,3(F15.8,2x))') latticeVec(:,1)
       write(*,'(1x,3(F15.8,2x))') latticeVec(:,2)
       write(*,'(1x,3(F15.8,2x))') latticeVec(:,3)
       write(*,*)
    end if

    if (v > 1) then
       frmt = '(1x,A' // trim(io_adjustl(TYPELEN)) // ',3(F15.8,2x)'
       if (hasForces) then
          frmt = trim(frmt) // '2x,3(F15.8,2x))'
          write(*,*) "Lattice Coordinates (Angstrom) and Forces (eV/Angstrom)"
       else
          frmt = trim(frmt) // ')'
          write(*,*) "Lattice Coordinates (Angstrom)"
       end if
       write(*,*)
       do iat = 1, nAtoms
          if (hasForces) then
             write(*,frmt) trim(adjustl(atomTypeName(atomType(iat)))), &
                           cooLatt(:,iat), forLatt(:,iat)
          else
             write(*,frmt) trim(adjustl(atomTypeName(atomType(iat)))), &
                           cooLatt(:,iat)
          end if
       end do
    end if

  end subroutine geo_print_info

  !--------------------------------------------------------------------!
  !                  cartesian coordinates and forces                  !
  !--------------------------------------------------------------------!

  function cooCart(iatom) result(coo)

    implicit none

    integer,            intent(in) :: iatom
    double precision, dimension(3) :: coo

    if (.not. isInit) then
       coo(:) = 0.0d0
       return
    end if

    coo(:) = matmul(latticeVec, cooLatt(1:3,iatom))

  end function cooCart

  !--------------------------------------------------------------------!

  function forCart(iatom) result(force)

    implicit none

    integer,            intent(in) :: iatom
    double precision, dimension(3) :: force

    if (.not. isInit) then
       force(:) = 0.0d0
       return
    end if

    force(:) = matmul(latticeVec, forLatt(1:3,iatom))

  end function forCart

  !--------------------------------------------------------------------!
  !                        write out structure                         !
  !--------------------------------------------------------------------!

  subroutine geo_write(cooLatt, atomType, typeName, latticeVec,   &
                       pbc, forCart, cohesiveEnergy, totalEnergy, &
                       comment, fileformat, file, unit)

    implicit none

    !------------------------------------------------------------------!
    ! cooLatt(i,j)    i-th lattice coordinate of the j-th atom,        !
    ! atomType(i)     the atom type ID of the i-th atom                !
    ! typeName(i)     name of the i-th atomic species                  !
    ! latticeVec(i,j) i-th Cartesian component of the j-th lattice vec.!
    ! pbc             if .true. the structure is periodic              !
    !---------------------------- optional ----------------------------!
    ! forCart(i,j)    i-th Cartesian component of the force acting on  !
    !                 the j-th atom                                    !
    ! cohesiveEnergy  value of the cohesive energy                     !
    ! totalEnergy     value of the total energy                        !
    ! comment         character string containing a comment            !
    ! frmt            output format (currently: 'xyz' or 'xsf')        !
    ! file            name of the output file (will be overwritten)    !
    ! unit            name of the output unit (formatted access)       !
    ! if neither file nor unit is given, output will be written to     !
    ! standard out                                                     !
    !------------------------------------------------------------------!

    double precision, dimension(:,:),           intent(in) :: cooLatt
    integer,          dimension(:),             intent(in) :: atomType
    character(len=*), dimension(:),             intent(in) :: typeName
    double precision, dimension(3,3),           intent(in) :: latticeVec
    logical,                                    intent(in) :: pbc
    double precision, dimension(:,:), optional, intent(in) :: forCart
    double precision,                 optional, intent(in) :: cohesiveEnergy
    double precision,                 optional, intent(in) :: totalEnergy
    character(len=*),                 optional, intent(in) :: comment
    character(len=*),                 optional, intent(in) :: fileformat
    character(len=*),                 optional, intent(in) :: file
    integer,                          optional, intent(in) :: unit

    character(len=3)       :: frmt
    integer                :: u
    double precision       :: E_coh, E_tot
    character(len=LINELEN) :: cmt
    integer                :: nAtoms, nTypes
    logical                :: dummy

    if (present(fileformat)) then
       frmt = trim(fileformat)
    else
       frmt = 'xyz'
    end  if

    E_coh = 0.0d0
    E_tot = 0.0d0
    if (present(cohesiveEnergy)) E_coh = cohesiveEnergy
    if (present(totalEnergy))    E_tot = totalEnergy

    if (present(comment)) then
       cmt = trim(comment)
    else
       cmt = ' '
    end if

    nAtoms = size(cooLatt(1,:))
    nTypes = size(typeName(:))

    if (present(unit)) then
       u = unit
    else if (present(file)) then
       u = io_unit()
       open(u, file=trim(file), status='replace', action='write')
    else
       u = STDOUT
    end if

    select case(frmt)
    case('xyz')
       if (present(forCart)) then
          call geo_write_xyz(u, nAtoms, nTypes, cooLatt, atomType,    &
                             typeName, latticeVec, E_coh, E_tot, cmt, &
                             forces=forCart)
       else
          call geo_write_xyz(u, nAtoms, nTypes, cooLatt, atomType,  &
                             typeName, latticeVec, E_coh, E_tot, cmt)
       end if
    case default
       write(0,*) "Error: invalid file format in 'geo_write': ", trim(frmt)
    end select

    if (present(file)) close(u)

    ! do something with currently unused arguments to silence gfortran
    dummy = pbc

  end subroutine geo_write


  !====================================================================!
  !                                                                    !
  !                        auxiliary procedures                        !
  !                                                                    !
  !====================================================================!

  !--------------------------------------------------------------------!
  !                retrieve type number from type name                 !
  !--------------------------------------------------------------------!

  function geo_itype_of_name(name, typeName) result(itype)

    implicit none

    character(len=*),               intent(in) :: name
    character(len=*), dimension(:), intent(in) :: typeName
    integer                                    :: itype

    integer :: itype1, ntypes

    ntypes = size(typeName)

    itype = 0
    do itype1 = 1, nTypes
       if (trim(typeName(itype1)) == trim(name)) then
          itype = itype1
          exit
       end if
    end do

  end function geo_itype_of_name

  !--------------------------------------------------------------------!
  !    change atom type IDs to match input data and structure info     !
  !--------------------------------------------------------------------!

  subroutine geo_match_atom_types(nTypes_ref, typeName_ref, nTypes_orig, &
                                  typeName_orig, nAtoms, atomType_orig, &
                                  typeName, atomType)

    implicit none

    integer,                                        intent(in)    :: nTypes_ref
    character(len=TYPELEN), dimension(nTypes_ref),  intent(in)    :: typeName_ref
    integer,                                        intent(in)    :: nTypes_orig
    character(len=TYPELEN), dimension(nTypes_orig), intent(in)    :: typeName_orig
    integer,                                        intent(in)    :: nAtoms
    integer,                dimension(nAtoms),      intent(in)    :: atomType_orig
    character(len=TYPELEN), dimension(nTypes_ref),  intent(out)   :: typeName
    integer,                dimension(nAtoms),      intent(out)   :: atomType

    integer                    :: iatom, itype
    character(len=TYPELEN)     :: name

    if (nTypes_orig > nTypes_ref) then
       write(0,*) "Error: more atomic types present than known in `geo_match_atom_types'."
       stop
    end if

    do iatom = 1, nAtoms
       name  = typeName_orig(atomType_orig(iatom))
       itype = geo_itype_of_name(name, typeName_ref)
       atomType(iatom) = itype
    end do
    typeName(:) = typeName_ref(:)

  end subroutine geo_match_atom_types

  !--------------------------------------------------------------------!
  !                        convert type number                         !
  !                                                                    !
  ! For the case of two different atom type indices:                   !
  ! (1) the atom types of the parametrization                          !
  ! (2) the atom types of the currently loaded structure               !
  ! This routine allows to convert between the two.                    !
  ! 0 will be returned, if the species is not found.                   !
  !--------------------------------------------------------------------!

  function geo_type_conv(it1, nT1, name1, nT2, name2) result(it2)

    implicit none

    integer,                          intent(in) :: it1
    integer,                          intent(in) :: nT1
    character(len=*), dimension(nT1), intent(in) :: name1
    integer,                          intent(in) :: nT2
    character(len=*), dimension(nT2), intent(in) :: name2
    integer                                      :: it2

    integer :: itype

    it2 = 0
    search : do itype = 1, nT2
       if (trim(name2(itype)) == trim(name1(it1))) then
          it2 = itype
          exit search
       end if
    end do search

  end function geo_type_conv


  !--------------------------------------------------------------------!
  !                            cell volume                             !
  !--------------------------------------------------------------------!

  function geo_cell_volume(avec) result(V)

    implicit none

    double precision, dimension(3,3), intent(in) :: avec
    double precision                             :: V

    V = avec(1,1)*avec(2,2)*avec(3,3) &
      + avec(2,1)*avec(3,2)*avec(1,3) &
      + avec(3,1)*avec(1,2)*avec(2,3) &
      - avec(3,1)*avec(2,2)*avec(1,3) &
      - avec(1,1)*avec(3,2)*avec(2,3) &
      - avec(2,1)*avec(1,2)*avec(3,3)

  end function geo_cell_volume

  !--------------------------------------------------------------------!
  !           calculation of the reciprocal lattice vectors            !
  !                                                                    !
  ! if (cryst == .true.) the crystallographic reciprocal lattice will  !
  ! be returned, i.e., the vectors are not scaled by 2*PI              !
  !--------------------------------------------------------------------!

  function geo_recip_lattice(avec, cryst) result(bvec)

    implicit none

    double precision, dimension(3,3), intent(in) :: avec
    double precision, dimension(3,3)             :: bvec
    logical, optional,                intent(in) :: cryst

    double precision :: V

    bvec(1,1:3) =  vproduct(avec(1:3,2), avec(1:3,3))
    bvec(2,1:3) = -vproduct(avec(1:3,1), avec(1:3,3))
    bvec(3,1:3) =  vproduct(avec(1:3,1), avec(1:3,2))

    V = geo_cell_volume(avec)
    bvec(:,:) = bvec(:,:)/V

    if (present(cryst)) then
       if (.not. cryst) bvec(:,:) = bvec(:,:)*2.0d0*PI
    else
       bvec(:,:) = bvec(:,:)*2.0d0*PI
    end if

  end function geo_recip_lattice

  !------------------------------------------------------------------!
  !                       vector/cross product                       !
  !------------------------------------------------------------------!

  function vproduct(a,b) result(c)

    implicit none

    double precision, dimension(3), intent(in) :: a, b
    double precision, dimension(3)             :: c

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function vproduct

  !--------------------------------------------------------------------!
  !        deduce bounds from coordinates of isolated structure        !
  ! Returns input coordinates shifted by 'shift'.                      !
  !--------------------------------------------------------------------!

  subroutine geo_get_bounds(cooCart, scal, shift)

    implicit none

    double precision, dimension(:,:), intent(inout) :: cooCart
    double precision, dimension(3,3), intent(out)   :: scal
    double precision, dimension(3),   intent(out)   :: shift

    double precision :: x_min, x_max
    double precision :: y_min, y_max
    double precision :: z_min, z_max

    integer :: iat

    ! To avoid numerical problems we make the bounding box
    ! slightly larger than necessary and also consider this in
    ! the shift of the coordinates. This should guarantee that all
    ! scaled coordinates are in [0,1).
    double precision, parameter :: EPS = 2.0d-6

    x_min = minval(cooCart(1,:))
    x_max = maxval(cooCart(1,:))
    y_min = minval(cooCart(2,:))
    y_max = maxval(cooCart(2,:))
    z_min = minval(cooCart(3,:))
    z_max = maxval(cooCart(3,:))

    ! orthogonal bounding box
    scal(:,1) = (/ x_max - x_min + EPS, 0.0d0, 0.0d0  /)
    scal(:,2) = (/ 0.0d0, y_max - y_min + EPS, 0.0d0  /)
    scal(:,3) = (/ 0.0d0, 0.0d0, z_max - z_min + EPS  /)

    ! origin of the bounding box
    shift(1) = x_min - 0.5d0*EPS
    shift(2) = y_min - 0.5d0*EPS
    shift(3) = z_min - 0.5d0*EPS

    ! shift coordinates to bounding box
    do iat = 1, size(cooCart(1,:))
       cooCart(1:3,iat) = cooCart(1:3,iat) - shift(1:3)
    end do

  end subroutine geo_get_bounds

  !--------------------------------------------------------------------!
  !           update bounding box (e.g. during optimization)           !
  !--------------------------------------------------------------------!

  subroutine geo_update_bounds(cooLatt, scal, rscal, shift)

    implicit none

    double precision, dimension(:,:), intent(inout) :: cooLatt
    double precision, dimension(3,3), intent(inout) :: scal
    double precision, dimension(3,3), intent(inout) :: rscal
    double precision, dimension(3),   intent(inout) :: shift

    double precision, dimension(3) :: new_scal
    double precision, dimension(3) :: new_shift

    double precision :: x_min, x_max
    double precision :: y_min, y_max
    double precision :: z_min, z_max

    integer :: iat

    double precision, parameter :: EPS = 2.0d-6

    if (.not. (any(cooLatt < 0.0d0) .or. (any(cooLatt > 1.0d0)))) return

    x_min = minval(cooLatt(1,:))
    x_max = maxval(cooLatt(1,:))
    y_min = minval(cooLatt(2,:))
    y_max = maxval(cooLatt(2,:))
    z_min = minval(cooLatt(3,:))
    z_max = maxval(cooLatt(3,:))

    new_scal(1) = x_max - x_min + EPS/scal(1,1)
    new_scal(2) = y_max - y_min + EPS/scal(2,2)
    new_scal(3) = z_max - z_min + EPS/scal(3,3)

    new_shift(1) = x_min - 0.5d0*EPS/scal(1,1)
    new_shift(2) = y_min - 0.5d0*EPS/scal(2,2)
    new_shift(3) = z_min - 0.5d0*EPS/scal(3,3)

    ! update scaled coordinates
    do iat = 1, size(cooLatt(1,:))
       cooLatt(1:3,iat) = (cooLatt(1:3,iat) - new_shift(1:3))/new_scal(1:3)
    end do

    ! update shift (origin of bounding box)
    shift(1) = shift(1) + new_shift(1)*scal(1,1)
    shift(2) = shift(2) + new_shift(2)*scal(2,2)
    shift(3) = shift(3) + new_shift(3)*scal(3,3)

    ! update diagonal scaling matrix (lattice vectors)
    scal(1,1) = scal(1,1)*new_scal(1)
    scal(2,2) = scal(2,2)*new_scal(2)
    scal(3,3) = scal(3,3)*new_scal(3)

    ! inverse scaling factors (reciprocal lattice vectors)
    rscal(1,1) = 1.0d0/scal(1,1)
    rscal(2,2) = 1.0d0/scal(2,2)
    rscal(3,3) = 1.0d0/scal(3,3)

  end subroutine geo_update_bounds


  !====================================================================!
  !                                                                    !
  !                 Input / Output of specific formats                 !
  !                                                                    !
  !====================================================================!


  !--------------------------------------------------------------------!
  !                             XSF format                             !
  !                                                                    !
  ! The length unit in the XSF format is Angstrom, so no need to       !
  ! convert anything.  Energies and forces are assumed to be in eV and !
  ! eV/Angstrom.                                                       !
  !--------------------------------------------------------------------!

  subroutine geo_init_xsf(file)

    use xsflib,  only: xsf_init,                        &
                       xsf_final,                       &
                       latticeVec_in   => latticeVec,   &
                       nAtoms_in       => nAtoms,       &
                       nTypes_in       => nTypes,       &
                       atomType_in     => atomType,     &
                       hasForces_in    => hasForces,    &
                       atomTypeName_in => atomTypeName, &
                       forCart_in      => forCart,      &
                       cooCart_in      => cooCart,      &
                       pbc_in          => pbc,          &
                       hasOutput,                       &
                       E_coh, E_tot

    implicit none

    character(len=*), intent(in) :: file

    call xsf_init(file)

    pbc = pbc_in
    if (pbc) then
       latticeVec(:,:) = latticeVec_in(:,:)
       origin = (/ 0.0d0, 0.0d0, 0.0d0 /)
    else
       call geo_get_bounds(cooCart_in, latticeVec, origin)
    end if
    recLattVec(:,:) = geo_recip_lattice(latticeVec)

    nAtoms    = nAtoms_in
    nTypes    = nTypes_in
    hasForces = hasForces_in

    if ((nAtoms == 0) .or. (nTypes == 0)) then
       write(0,*) "Error: invalid input file: ", trim(file)
       write(0,*) "       nAtoms = ", trim(io_adjustl(nAtoms))
       write(0,*) "       nTypes = ", trim(io_adjustl(nTypes))
       call xsf_final()
       stop
    end if

    if (hasForces) then
       allocate(cooLatt(3,nAtoms), forLatt(3,nAtoms), &
                atomType(nAtoms), atomTypeName(nTypes))
    else
       allocate(cooLatt(3,nAtoms), atomType(nAtoms),  &
                atomTypeName(nTypes))
    end if

    atomType(:)     = atomType_in(:)
    atomTypeName(:) = atomTypeName_in(:)

    ! convert cartesian coordinates to lattice coordinates:
    cooLatt(1:3,1:nAtoms) = matmul(recLattVec, cooCart_in)/(2.0d0*PI)

    ! if forces are available, convert them too:
    if (hasForces) then
       forLatt(1:3,1:nAtoms) = matmul(recLattVec, forCart_in)/(2.0d0*PI)
    end if

    ! store cohesive energy, if available
    if (hasOutput) then
       cohesiveEnergy = E_coh
       totalEnergy    = E_tot
       hasEnergy      = .true.
    else
       cohesiveEnergy = 0.0d0
       totalEnergy    = 0.0d0
       hasEnergy      = .false.
    end if

    call xsf_final

  end subroutine geo_init_xsf


  !--------------------------------------------------------------------!
  !                  simple XYZ format (output only)                   !
  !--------------------------------------------------------------------!

  subroutine geo_write_xyz(u, nAtoms, nTypes, cooLatt, atomType,    &
                           typeName, latticeVec, E_coh, E_tot, cmt, &
                           forces)

    implicit none

    integer,                                         intent(in) :: u
    integer,                                         intent(in) :: nAtoms
    integer,                                         intent(in) :: nTypes
    double precision, dimension(3,nAtoms),           intent(in) :: cooLatt
    integer,          dimension(nAtoms),             intent(in) :: atomType
    character(len=*), dimension(nTypes),             intent(in) :: typeName
    double precision, dimension(3,3),                intent(in) :: latticeVec
    double precision,                                intent(in) :: E_coh, E_tot
    character(len=*),                                intent(in) :: cmt
    double precision, dimension(3,nAtoms), optional, intent(in) :: forces

    integer                        :: iat
    character(len=2)               :: symbol
    double precision, dimension(3) :: cooCart
    logical                        :: dummy

    write(u,*) nAtoms
    write(u,*) trim(cmt)
    do iat = 1, nAtoms
       symbol = typeName(atomType(iat))
       cooCart(1:3) = matmul(latticeVec, cooLatt(1:3,iat))
       write(u,'(1x,A2,3(2x,F12.6))') symbol, cooCart(1:3)
    end do

    ! do something with surrently unused arguments to silence gfortran
    if (E_coh > E_tot) dummy = .true.
    if (present(forces)) dummy = .true.

  end subroutine geo_write_xyz

end module geometry
