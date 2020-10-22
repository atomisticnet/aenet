!----------------------------------------------------------------------
!      lclist.f90 - implementation of a linked cell neighbor list
!----------------------------------------------------------------------
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
!----------------------------------------------------------------------

module lclist

  !--------------------------------------------------------------------!
  ! A simple, but universal neighbourlist implementation without any   !
  ! restrictions on the lattice vectors or on the number of atoms in   !
  ! the system.  It is equally well suited for small unit cells, where !
  ! |lattice vectors| << R_cut, and for unit cells that are much       !
  ! larger than the interaction radius.                                !
  !                                                                    !
  ! If the simulation cell is large wrt. the interaction cut-off, it   !
  ! is partitioned into cells and the information about the atoms in   !
  ! each cell is stored in a linked cell list.                         !
  !--------------------------------------------------------------------!
  ! 2011-11-05 Alexander Urban (AU)                                    !
  !--------------------------------------------------------------------!

  use sortlib, only: argsort

  implicit none
  save

  public  :: lcl_init,               &
             lcl_final,              &
             lcl_print_info,         &
             lcl_nmax_cell,          &
             lcl_nmax_nblist,        &
             lcl_nmax_nbdist,        &
             lcl_nblist,             &
             lcl_nbdist,             &
             lcl_nbdist_cart

  private :: cell_assign_atoms,      &
             cell_multiples,         &
             cell_get,               &
             cell_add,               &
             wrap_cooLatt,           &
             cell_volume,            &
             inner_sphere,           &
             vproduct,               &
             translation_vectors


  double precision, parameter, private :: PI = 3.1415926535897931d0

  !--------------------------------------------------------------------!

  logical,                                         private :: pbc
  double precision, dimension(3,3),                private :: latticeVec
  double precision, dimension(3,3),                private :: gridVec
  integer,                                         private :: nAtoms
  integer,          dimension(:),   pointer,       private :: atomType
  double precision, dimension(:,:), pointer,       private :: cooLatt

  double precision,                                private :: r_min
  double precision,                                private :: r_max
  integer,          dimension(3),                  private :: nCells
  integer,          dimension(:,:,:), allocatable, private :: cell
  integer,          dimension(:),     allocatable, private :: atomList
  integer,          dimension(:,:),   allocatable, private :: cellList

  integer,                                         private :: nCvecs
  integer,          dimension(:,:),   allocatable, private :: Cvec

  integer,                                         private :: nTvecs
  integer,          dimension(:,:),   allocatable, private :: Tvec

  integer,                                         private :: nmax_cell
  integer,                                         private :: nmax_nblist
  integer,                                         private :: nmax_nbdist

  logical,                                         private :: isInit = .false.

contains

  subroutine lcl_init(r_min_in, r_max_in, latticeVec_in, &
                      nAtoms_in, atomType_in, cooLatt_in, pbc_in)

    implicit none

    double precision,                                 intent(in) :: r_min_in
    double precision,                                 intent(in) :: r_max_in
    double precision, dimension(3,3),                 intent(in) :: latticeVec_in
    integer,                                          intent(in) :: nAtoms_in
    integer,          dimension(nAtoms_in),   target, intent(in) :: atomType_in
    double precision, dimension(3,nAtoms_in), target, intent(in) :: cooLatt_in
    logical, optional,                                intent(in) :: pbc_in

    if (isInit) then
       write(0,*) "Error: module already initialized in `lclist_init'."
       return
    end if

    r_min           = r_min_in
    r_max           = r_max_in
    nAtoms          = nAtoms_in
    latticeVec(:,:) = latticeVec_in(:,:)
    atomType        => atomType_in(:)
    cooLatt         => cooLatt_in(:,:)
    if (present(pbc_in)) then
       pbc = pbc_in
    else
       pbc = .true.
    end if

    ! get optimal divison of the unit cell:
    call cell_multiples(r_max, latticeVec, gridVec, nCells)
    allocate(cell(nCells(3),nCells(2),nCells(1)), atomList(nAtoms), &
             cellList(3,nAtoms))

    ! assign atoms to cells:
    call cell_assign_atoms(nCells, nAtoms, cooLatt, cell, atomList, cellList)

    ! set up half start of vectors pointing to cells within range of r_max:
    nCvecs = 0
    call translation_vectors(r_max, gridVec, nCvecs, Cvec, nc=nCells, pbc=pbc)
    allocate(Cvec(3,nCvecs))
    call translation_vectors(r_max, gridVec, nCVecs, cVec, nc=nCells, pbc=pbc)

    ! set up half start of translation vectors pointing to periodic
    ! images of the unit cell within range of r_max:
    nTvecs = 0
    if (pbc) then
       call translation_vectors(r_max, latticeVec, nTvecs, Tvec)
       allocate(Tvec(3,nTvecs))
       call translation_vectors(r_max, latticeVec, nTvecs, Tvec)
    end if

    isInit = .true.

    ! max. number of neighbours and max. number of atoms per cell:
    nmax_cell   = lcl_nmax_cell()
    nmax_nblist = lcl_nmax_nblist()
    nmax_nbdist = lcl_nmax_nbdist(r_min, r_max)

  end subroutine lcl_init

  !--------------------------------------------------------------------!

  subroutine lcl_final()

    implicit none

    if (.not. isInit) return

    if (allocated(cell)) deallocate(cell, atomList, cellList)
    if (allocated(Cvec)) deallocate(Cvec)
    if (allocated(Tvec)) deallocate(Tvec)

    cooLatt  => null()
    atomType => null()

    isInit = .false.

  end subroutine lcl_final

  !--------------------------------------------------------------------!
  !       print information about the linked cell list to stdout       !
  !--------------------------------------------------------------------!

  subroutine lcl_print_info()

    implicit none

    write(*,*) 'Linked Cell List (Neighbourlist)'
    write(*,*) '--------------------------------'
    write(*,*)

    if (.not. isInit) then
       write(*,*) "Module not initialized."
       write(*,*)
       return
    end if

    write(*,'(1x,"Number of atoms             : ",I10)')      nAtoms
    write(*,'(1x,"Minimal distance            : ",ES10.3)')   r_min
    write(*,'(1x,"Cut-off radius              : ",ES10.3)')   r_max
    write(*,'(1x,"Number of cells             : ",3(I4,2x))') nCells
    write(*,'(1x,"Cells within cut-off        : ",I10)')      nCvecs + 1
    write(*,'(1x,"Translation vectors         : ",I10)')      nTvecs + 1
    write(*,'(1x,"Max. atoms per cell         : ",I10)')      nmax_cell
    write(*,'(1x,"Max. possible neighbour IDs : ",I10)')      nmax_nblist
    write(*,'(1x,"Max. real neighbours        : ",I10)')      nmax_nbdist
    write(*,*)

  end subroutine lcl_print_info

  !--------------------------------------------------------------------!
  !        max. number of neighbour candidates within cut-off          !
  !--------------------------------------------------------------------!

  function lcl_nmax_nblist() result(nmax)

    implicit none

    integer          :: nmax

    nmax = min(lcl_nmax_cell()*(nCvecs + 1), nAtoms-1)

  end function lcl_nmax_nblist

  !--------------------------------------------------------------------!
  !          max. number of (real) neighbours within cut-off           !
  !--------------------------------------------------------------------!

  function lcl_nmax_nbdist(rmin, rmax) result(nmax)

    implicit none

    double precision, intent(in) :: rmin, rmax
    integer                      :: nmax

    double precision :: V_atom, V_cut

    V_atom = (0.5d0*rmin)**3
    V_cut  = (rmax+0.5d0*rmin)**3

    ! max number of atoms assuming close packing
    ! pi/(s*sqrt(2)) ~ 0.7405
    nmax = ceiling(V_cut/V_atom*0.7405d0)

  end function lcl_nmax_nbdist

  !--------------------------------------------------------------------!
  !           estimated max. average number of atoms per cell          !
  !--------------------------------------------------------------------!

  function lcl_nmax_cell() result(nmax)

    implicit none

    integer          :: nmax
    double precision :: V_cell, V_atom

    if (.not. isInit) then
       write(0,*) "Error: module not initialized in `max_atoms_per_cell'."
       stop
    end if

    V_cell = cell_volume(gridVec)
    ! effective atom volume in a densly packed crystal structure:
    ! 3*sqrt(2)/pi * 4/3*pi * r^3 = 4*sqrt(2) * r^3
    V_atom = sqrt(32.0d0)*(0.5d0*r_min)**3
    ! wrong?! assuming the also partial spheres in the cell.
    ! V_atom = 8.0d0/27.0d0*(0.5d0*r_min)**3
    nmax   = ceiling(V_cell/V_atom)

  end function lcl_nmax_cell

  !--------------------------------------------------------------------!
  !                      retrieve neighbour list                       !
  !                                                                    !
  ! This routine returns the complete list of atom numbers of possible !
  ! neighbours.  It is not checked, if the returned atoms really are   !
  ! within the cut-off.                                                !
  !--------------------------------------------------------------------!

  subroutine lcl_nblist(iatom, nnb, nblist, stat)

    implicit none

    integer,                 intent(in)    :: iatom
    integer,                 intent(inout) :: nnb
    integer, dimension(nnb), intent(out)   :: nblist
    integer,       optional, intent(out)   :: stat

    integer, dimension(3) :: ic0, ic
    integer               :: inb, iat
    integer               :: iv

    if (.not. isInit) then
       write(0,*) "Error: module not initialized in `nblist'."
       stop
    end if

    if ((iatom < 1) .or. (iatom > nAtoms)) then
       write(0,*) "Error: invalid atom number in `nblist' : ", iatom
       stop
    end if

    if ((nnb < nmax_nblist) .and. (.not. present(stat))) then
       write(0,*) "Error: array nblist too small in `nblist'."
       stop
    end if

    if (present(stat)) stat = 0
    if (nnb == 0) return

    ! iatom is in cell ic0:
    ic0(1:3) = cellList(1:3,iatom)

    ! (1) other atoms in the same cell:

    inb = 0
    iat = cell(ic0(3),ic0(2),ic0(1))
    if (iat /= 0) then
       if (iat /= iatom) then
          inb = inb + 1
          if (inb > nnb) then
             if (present(stat)) stat = 1
             return
          end if
          nblist(inb) = iat
       end if
       do while(atomList(iat) /= 0)
          iat = atomList(iat)
          if (iat == iatom) cycle
          inb = inb + 1
          if (inb > nnb) then
             if (present(stat)) stat = 2
             return
          end if
          nblist(inb) = iat
       end do
    end if

    ! (2) atoms in neighbouring cells:

    ivloop : do iv = 1, nCvecs

       ! get cell for cell vector
       ic(1:3) = cell_from_Cvec(ic0, Cvec(1:3,iv), nCells, pbc)
       if (any(ic < 1) .or. any(ic > nCells)) cycle ivloop

       ! all atoms in that cell:
       iat = cell(ic(3),ic(2),ic(1))
       if (iat /= 0) then
          if (.not. any(nblist(1:inb) == iat)) then
             inb = inb + 1
             if (inb > nnb) then
                if (present(stat)) stat = 3
                return
             end if
             nblist(inb) = iat
          end if
          do while(atomList(iat) /= 0)
             iat = atomList(iat)
             if (.not. any(nblist(1:inb) == iat)) then
                inb = inb + 1
                if (inb > nnb) then
                   if (present(stat)) stat = 4
                   return
                end if
                nblist(inb) = iat
             end if
          end do
       end if

    end do ivloop

    nnb = inb

  end subroutine lcl_nblist

  !--------------------------------------------------------------------!
  !      retrieve real PBC coordinates of the neighbouring atoms       !
  !--------------------------------------------------------------------!

  subroutine lcl_nbdist(iatom, nnb, nbcoo, nbdist, r_cut, itype_opt, stat)

    implicit none

    integer,                            intent(in)    :: iatom
    integer,                            intent(inout) :: nnb
    double precision, dimension(3,nnb), intent(out)   :: nbcoo
    double precision, dimension(nnb),   intent(out)   :: nbdist
    double precision, optional,         intent(in)    :: r_cut
    integer,          optional,         intent(in)    :: itype_opt
    integer,          optional,         intent(out)   :: stat

    integer,          dimension(nmax_nblist)          :: nblist_loc
    integer                                           :: nnb_tot, inb, nnb2
    integer                                           :: iat, iT
    integer                                           :: itype
    double precision                                  :: Rc, Rc2, dist2
    double precision, dimension(3)                    :: coo2, cart

    if (.not. isInit) then
       write(0,*) "Error: module not initialized in `nbdist'."
       stop
    end if

    if ((iatom < 1) .or. (iatom > nAtoms)) then
       write(0,*) "Error: invalid atom number in `nbdist' : ", iatom
       stop
    end if

    if ((nnb < nmax_nbdist) .and. (.not. present(stat))) then
       write(0,*) "Error: array nbcoo too small in `nbdist'. (1)"
       stop
    end if

    if (present(itype_opt)) then
       itype = itype_opt
    else
       itype = -1
    end if

    if (present(stat)) stat = 0
    if (present(r_cut)) then
       Rc = r_cut
    else
       Rc = r_max
    end if
    Rc2 = Rc*Rc
    nnb_tot = 0

    ! (1) get neighbour list:

    nnb2 = nmax_nblist
    call lcl_nblist(iatom, nnb2, nblist_loc)

    ! (2) check distance do the periodic images of the central atom:

    nnb_tot = 0
    if ( (.not. present(itype_opt)) .or. (atomType(iatom) == itype)) then
       do iT = 1, nTvecs
          cart(1:3) = matmul(latticeVec, dble(Tvec(1:3,iT)))
          dist2 = sum(cart*cart)
          if (dist2 <= Rc2) then
             nnb_tot = nnb_tot + 1
             if (nnb_tot > nnb) then
                if (present(stat)) then
                   stat = 1
                   return
                else
                   write(0,*) "Error: array nbcoo too small in `nbdist'. (2)"
                   stop
                end if
             end if ! nbcoo too small
             nbcoo(1,nnb_tot) = cooLatt(1,iatom) + dble(Tvec(1,iT))
             nbcoo(2,nnb_tot) = cooLatt(2,iatom) + dble(Tvec(2,iT))
             nbcoo(3,nnb_tot) = cooLatt(3,iatom) + dble(Tvec(3,iT))
             nbdist(nnb_tot)  = sqrt(dist2)
          end if
       end do
    end if ! correct atom type

    ! (3) calculate PBC distances for all atoms in the neighborlist:

    do inb = 1, nnb2
       iat = nblist_loc(inb)
       if (present(itype_opt) .and. (atomType(iat) /= itype)) cycle
       ! in home unit cell:
       cart(1:3) = cooLatt(1:3,iat) - cooLatt(1:3,iatom)
       cart(1:3) = matmul(latticeVec, cart(1:3))
       dist2 = sum(cart*cart)
       if (dist2 <= Rc2) then
          nnb_tot = nnb_tot + 1
          if (nnb_tot > nnb) then
             if (present(stat)) then
                stat = 2
                return
             else
                write(0,*) "Error: array nbcoo too small in `nbdist'. (3)"
                stop
             end if
          end if ! nbcoo too small
          nbcoo(1:3,nnb_tot) = cooLatt(1:3,iat)
          nbdist(nnb_tot)    = sqrt(dist2)
       end if ! within cut-off
       ! in periodic images:
       do iT = 1, nTvecs
          coo2(1) = cooLatt(1,iat) + dble(Tvec(1,iT))
          coo2(2) = cooLatt(2,iat) + dble(Tvec(2,iT))
          coo2(3) = cooLatt(3,iat) + dble(Tvec(3,iT))
          cart(1:3) = coo2(1:3) - cooLatt(1:3,iatom)
          cart(1:3) = matmul(latticeVec, cart(1:3))
          dist2 = sum(cart*cart)
          if (dist2 <= Rc2) then
             nnb_tot = nnb_tot + 1
             if (nnb_tot > nnb) then
                if (present(stat)) then
                   stat = 3
                   return
                else
                   write(0,*) "Error: array nbcoo too small in `ndist'. (4)"
                   stop
                end if
             end if ! nbcoo too small
             nbcoo(1:3,nnb_tot) = coo2(1:3)
             nbdist(nnb_tot)    = sqrt(dist2)
          end if ! within cut-off
       end do
    end do

    ! (4) return number of neighbours found:

    nnb = nnb_tot

  end subroutine lcl_nbdist

  !--------------------------------------------------------------------!
  !    same as `lcl_nbdist', but cartesian coordinates are returned    !
  !--------------------------------------------------------------------!

  subroutine lcl_nbdist_cart(iatom, nnb, nbcoo, nbdist, r_cut, itype_opt, &
                             nblist, nbtype, stat)

    implicit none

    integer,                                    intent(in)    :: iatom
    integer,                                    intent(inout) :: nnb
    double precision, dimension(3,nnb),         intent(out)   :: nbcoo
    double precision, dimension(nnb),           intent(out)   :: nbdist
    double precision,                 optional, intent(in)    :: r_cut
    integer,                          optional, intent(in)    :: itype_opt
    integer,          dimension(nnb), optional, intent(out)   :: nblist
    integer,          dimension(nnb), optional, intent(out)   :: nbtype
    integer,                          optional, intent(out)   :: stat

    double precision,      parameter :: EPS = 1.0d-3

    integer,          dimension(nmax_nblist) :: nblist_loc
    integer                                  :: nnb_tot, inb, nnb2
    integer                                  :: iat, iT
    integer                                  :: itype
    double precision                         :: Rc, Rc2, dist2
    double precision, dimension(3)           :: coo1, coo2, cart
    integer                                  :: nblist_stat

    if (.not. isInit) then
       write(0,*) "Error: module not initialized in `nbdist'."
       stop
    end if

    if ((iatom < 1) .or. (iatom > nAtoms)) then
       write(0,*) "Error: invalid atom number in `nbdist' : ", iatom
       stop
    end if

    if ((nnb < nmax_nbdist) .and. (.not. present(stat))) then
       write(0,*) "Error: array nbcoo too small in `nbdist_cart'. (1)"
       stop
    end if

    if (present(itype_opt)) then
       itype = itype_opt
    else
       itype = -1
    end if

    if (present(stat)) stat = 0
    if (present(r_cut)) then
       Rc = r_cut
    else
       Rc = r_max
    end if
    Rc2 = (Rc + EPS)*(Rc + EPS)
    nnb_tot = 0

    ! (1) get neighbour list:

    nnb2 = nmax_nblist
    call lcl_nblist(iatom, nnb2, nblist_loc, stat=nblist_stat)
    if (nblist_stat /= 0) then
       write(0,*) "Error: Neighbor list returned with status ", nblist_stat
       if (present(stat)) then
          stat = 4
          return
       else
          stop
       end if
    end if

    ! (2) check distance do the periodic images of the central atom:

    nnb_tot = 0
    if ( (.not. present(itype_opt)) .or. (atomType(iatom) == itype)) then
       do iT = 1, nTvecs
          cart(1:3) = matmul(latticeVec, dble(Tvec(1:3,iT)))
          dist2 = sum(cart*cart)
          cutoff : if (dist2 <= Rc2) then
             nnb_tot = nnb_tot + 1
             if (nnb_tot > nnb) then
                if (present(stat)) then
                   stat = 1
                   return
                else
                   write(0,*) "Error: array nbcoo too small in `nbdist_cart'. (2)"
                   stop
                end if
             end if ! nbcoo too small
             if (present(nblist)) nblist(nnb_tot) = iatom
             if (present(nbtype)) nbtype(nnb_tot) = atomType(iatom)
             nbcoo(1,nnb_tot) = cooLatt(1,iatom) + dble(Tvec(1,iT))
             nbcoo(2,nnb_tot) = cooLatt(2,iatom) + dble(Tvec(2,iT))
             nbcoo(3,nnb_tot) = cooLatt(3,iatom) + dble(Tvec(3,iT))
             nbcoo(:,nnb_tot) = matmul(latticeVec, nbcoo(:,nnb_tot))
             nbdist(nnb_tot)  = sqrt(dist2)
          end if cutoff
       end do
    end if ! correct atom type

    ! (3) calculate PBC distances for all atoms in the neighbourlist:

    ! cartesian coordinates of central atom:
    coo1(1:3) = matmul(latticeVec, cooLatt(1:3,iatom))

    do inb = 1, nnb2
       iat = nblist_loc(inb)
       if (present(itype_opt) .and. (atomType(iat) /= itype)) cycle
       ! in home unit cell:
       coo2(1:3) = matmul(latticeVec, cooLatt(1:3,iat))
       cart(1:3) = coo2(1:3) - coo1(1:3)
       dist2 = sum(cart*cart)
       if (dist2 <= Rc2) then
          nnb_tot = nnb_tot + 1
          if (nnb_tot > nnb) then
             if (present(stat)) then
                stat = 2
                return
             else
                write(0,*) "Error: array nbcoo too small in `nbdist_cart'. (3)"
                stop
             end if
          end if ! nbcoo too small
          if(present(nblist)) nblist(nnb_tot) = iat
          if(present(nbtype)) nbtype(nnb_tot) = atomType(iat)
          nbcoo(1:3,nnb_tot) = coo2(1:3)
          nbdist(nnb_tot)    = sqrt(dist2)
       end if ! within cut-off
       ! in periodic images:
       do iT = 1, nTvecs
          coo2(1) = cooLatt(1,iat) + dble(Tvec(1,iT))
          coo2(2) = cooLatt(2,iat) + dble(Tvec(2,iT))
          coo2(3) = cooLatt(3,iat) + dble(Tvec(3,iT))
          coo2(:) = matmul(latticeVec, coo2(1:3))
          cart(:) = coo2(1:3) - coo1(1:3)
          dist2 = sum(cart*cart)
          if (dist2 <= Rc2) then
             nnb_tot = nnb_tot + 1
             if (nnb_tot > nnb) then
                if (present(stat)) then
                   stat = 3
                   return
                else
                   write(0,*) "Error: array nbcoo too small in `nbdist_cart'. (4)"
                   stop
                end if
             end if ! nbcoo too small
             if(present(nblist)) nblist(nnb_tot) = iat
             if(present(nbtype)) nbtype(nnb_tot) = atomType(iat)
             nbcoo(1:3,nnb_tot) = coo2(1:3)
             nbdist(nnb_tot)    = sqrt(dist2)
          end if ! within cut-off
       end do
    end do

    ! (4) return number of neighbours found:

    nnb = nnb_tot

  end subroutine lcl_nbdist_cart


  !====================================================================!
  !                                                                    !
  !                         private procedures                         !
  !                                                                    !
  !====================================================================!



  !--------------------------------------------------------------------!
  !                       assign atoms to cells                        !
  !--------------------------------------------------------------------!

  subroutine cell_assign_atoms(nc, nAtoms, cooLatt, cell, atomList, &
                               cellList)

    implicit none

    integer,          dimension(3),                 intent(in)    :: nc
    integer,                                        intent(in)    :: nAtoms
    double precision, dimension(3,nAtoms),          intent(inout) :: cooLatt
    integer,          dimension(nc(3),nc(2),nc(1)), intent(out)   :: cell
    integer,          dimension(nAtoms),            intent(out)   :: atomList
    integer,          dimension(3,nAtoms),          intent(out)   :: cellList

    integer               :: iatom
    integer, dimension(3) :: ic

    cell(:,:,:)   = 0
    atomList(:)   = 0
    cellList(:,:) = 0

    do iatom = 1, nAtoms
       ! coordinates of isolated structures will be checked,
       ! but will not be wrapped
       call wrap_cooLatt(cooLatt(1:3,iatom))
       ic(1:3) = cell_get(cooLatt(1:3,iatom), nc)
       call cell_add(iatom, ic, nAtoms, nc, cell, atomList, cellList)
    end do

  end subroutine cell_assign_atoms

  !--------------------------------------------------------------------!
  !                determine the optimal cell division                 !
  !--------------------------------------------------------------------!

  subroutine cell_multiples(r_max, latticeVec, gridVec, nCells)

    implicit none

    double precision,                 intent(in)  :: r_max
    double precision, dimension(3,3), intent(in)  :: latticeVec
    double precision, dimension(3,3), intent(out) :: gridVec
    integer,          dimension(3),   intent(out) :: nCells

    double precision, dimension(3)   :: v
    double precision, dimension(3)   :: h
    double precision                 :: h_min
    double precision                 :: r

    ! compute the three heights:

    v  = vproduct(latticeVec(:,1),latticeVec(:,2))
    v  = v/sqrt(sum(v*v))
    h(3) = abs(sum(v*latticeVec(:,3)))

    v  = vproduct(latticeVec(:,1),latticeVec(:,3))
    v  = v/sqrt(sum(v*v))
    h(2) = abs(sum(v*latticeVec(:,2)))

    v  = vproduct(latticeVec(:,2),latticeVec(:,3))
    v  = v/sqrt(sum(v*v))
    h(1) = abs(sum(v*latticeVec(:,1)))

    ! shortest height, use minimal value for low-dimensional structures:

    h_min = max(minval(h), 0.5d0)

    ! compute multiples of h_min wrt. each height, but guarantee that
    ! at least one cell is used in each direction:

    nCells(1) = max(floor(h(1)/h_min), 1)
    nCells(2) = max(floor(h(2)/h_min), 1)
    nCells(3) = max(floor(h(3)/h_min), 1)

    ! compute largest sphere that fits into such a cell:

    gridVec(:,1) = latticeVec(:,1)/dble(nCells(1))
    gridVec(:,2) = latticeVec(:,2)/dble(nCells(2))
    gridVec(:,3) = latticeVec(:,3)/dble(nCells(3))
    call inner_sphere(gridVec, r)

    ! reduce cell size if that sphere is unnecessary large:

    if (r > 0.5d0*r_max) then
       nCells(1) = floor(nCells(1)*(2.0d0*r/r_max))
       nCells(2) = floor(nCells(2)*(2.0d0*r/r_max))
       nCells(3) = floor(nCells(3)*(2.0d0*r/r_max))
       gridVec(:,1) = latticeVec(:,1)/dble(nCells(1))
       gridVec(:,2) = latticeVec(:,2)/dble(nCells(2))
       gridVec(:,3) = latticeVec(:,3)/dble(nCells(3))
    end if

  end subroutine cell_multiples

  !--------------------------------------------------------------------!
  ! Remove redundant translation vectors that point to the same cells  !
  !--------------------------------------------------------------------!

  subroutine delete_redundant_cell_vecs(Cvec, nCvecs)

    implicit none

    integer, dimension(:,:), intent(inout) :: Cvec
    integer,                 intent(inout) :: nCvecs

    integer, dimension(3,nCvecs) :: Cvec_keep
    integer, dimension(3,nCvecs) :: cells
    integer, dimension(3)        :: ic
    integer                      :: i, nc

    if (.not. pbc) return

    nc = 0
    do i = 1, nCvecs
       ic(1:3) = cell_from_Cvec([1,1,1], Cvec(1:3,i), nCells, .true.)
       write(0,*) i, ic
       if (.not. (all(ic == 1) .or. in_list(ic, cells(1:3,1:nc)))) then
          nc = nc + 1
          Cvec_keep(1:3,nc) = Cvec(1:3,i)
          cells(1:3,nc) = ic(1:3)
       end if
    end do

    write(0,*) "Original number of cell vectors: ", nCvecs
    write(0,*) "Reduced number of cell vectors: ", nc

    Cvec(1:3,1:nc) = Cvec_keep(1:3,1:nc)
    nCvecs = nc

  end subroutine delete_redundant_cell_vecs



  !====================================================================!
  !                                                                    !
  !                 operations on the linked cell list                 !
  !                                                                    !
  !====================================================================!

  !--------------------------------------------------------------------!
  !             get cell for specific lattice coordinates              !
  !--------------------------------------------------------------------!

  function cell_get(coo, nc) result(ic)

    implicit none

    double precision, dimension(3), intent(in) :: coo
    integer,          dimension(3), intent(in) :: nc
    integer,          dimension(3)             :: ic

    ic(1) = min(nint(coo(1)*dble(nc(1)) + 0.5d0), nc(1))
    ic(2) = min(nint(coo(2)*dble(nc(2)) + 0.5d0), nc(2))
    ic(3) = min(nint(coo(3)*dble(nc(3)) + 0.5d0), nc(3))

  end function cell_get

  !--------------------------------------------------------------------!
  !                get cell relative to reference cell                 !
  ! Returns [0,0,0] if the cell is out of bounds (isolated structures  !
  ! only).                                                             !
  !--------------------------------------------------------------------!

  function cell_from_Cvec(ref_cell, Cvec_loc, nCells_loc, pbc) result(ic)

    implicit none

    !------------------------------------------------------------------!
    ! ref_cell        reference cell with                              !
    !                       1 <= ref_cell(i) <= nCells_loc(i)          !
    ! Cvec_loc(i)     i-th component of a cell translation vector      !
    ! nCells_loc(i)   total number of cells in direction i             !
    ! pbc             if .true., impose periodic poundary conditions   !
    !------------------------------------------------------------------!

    integer, dimension(3), intent(in) :: ref_cell
    integer, dimension(3), intent(in) :: Cvec_loc
    integer, dimension(3), intent(in) :: nCells_loc
    logical,               intent(in) :: pbc

    integer, dimension(3) :: ic

    if (pbc) then
       ic(1) = modulo(ref_cell(1) + Cvec_loc(1) - 1, nCells_loc(1)) + 1
       ic(2) = modulo(ref_cell(2) + Cvec_loc(2) - 1, nCells_loc(2)) + 1
       ic(3) = modulo(ref_cell(3) + Cvec_loc(3) - 1, nCells_loc(3)) + 1
    else
       ic(1:3) = ref_cell(1:3) + Cvec_loc(1:3)
       if (any(ic > nCells_loc) .or. any(ic < 1)) ic(1:3) = 0
    end if

  end function cell_from_Cvec

  !--------------------------------------------------------------------!
  !                          add atom to cell                          !
  !--------------------------------------------------------------------!

  subroutine cell_add(iatom, ic, nAtoms, nc, cell, atomList, cellList)

    implicit none

    integer,                               intent(in)    :: iatom
    integer, dimension(3),                 intent(in)    :: ic
    integer,                               intent(in)    :: nAtoms
    integer, dimension(3),                 intent(in)    :: nc
    integer, dimension(nc(3),nc(2),nc(1)), intent(inout) :: cell
    integer, dimension(nAtoms),            intent(inout) :: atomList
    integer, dimension(3,nAtoms),          intent(inout) :: cellList

    atomList(iatom) = cell(ic(3),ic(2),ic(1))
    cell(ic(3),ic(2),ic(1)) = iatom

    cellList(1:3,iatom) = ic(1:3)

  end subroutine cell_add

  !--------------------------------------------------------------------!
  !                       remove atom from cell                        !
  !--------------------------------------------------------------------!

!!$  subroutine cell_del(iatom, ic, nAtoms, nc, cell, atomList)
!!$
!!$    implicit none
!!$
!!$    integer,                               intent(in)    :: iatom
!!$    integer, dimension(3),                 intent(in)    :: ic
!!$    integer,                               intent(in)    :: nAtoms
!!$    integer, dimension(3),                 intent(in)    :: nc
!!$    integer, dimension(nc(3),nc(2),nc(1)), intent(inout) :: cell
!!$    integer, dimension(nAtoms),            intent(inout) :: atomList
!!$
!!$    integer :: i
!!$
!!$    if (cell(ic(3),ic(2),ic(1)) == iatom) then
!!$       cell(ic(3),ic(2),ic(1)) = atomList(iatom)
!!$    else
!!$       i = cell(ic(3),ic(2),ic(1))
!!$       do while(atomList(i) /= iatom)
!!$          if (i == 0) then
!!$             write(0,*) "Warning: atom not in cell in `cell_del'."
!!$             return
!!$          end if
!!$          i = atomList(i)
!!$       end do
!!$       atomList(i) = atomlist(iatom)
!!$    end if
!!$
!!$  end subroutine cell_del



  !====================================================================!
  !                                                                    !
  !                        auxilliary functions                        !
  !                                                                    !
  !====================================================================!



  !--------------------------------------------------------------------!
  !          wrap lattice coordinates to the interval [0,1[
  !--------------------------------------------------------------------!

  subroutine wrap_cooLatt(coo)

    implicit none

    double precision, dimension(3), intent(inout) :: coo

    if (.not. pbc) then
       if (any(coo < 0.0d0) .or. any(coo > 1.0d0)) then
          write(0,*) "Warning: atoms outside of bounding box in LC list!"
       end if
       return
    end if

    do while(coo(1) < 0.0d0)
       coo(1) = coo(1) + 1.0d0
    end do
    do while(coo(1) >= 1.0d0)
       coo(1) = coo(1) - 1.0d0
    end do

    do while(coo(2) < 0.0d0)
       coo(2) = coo(2) + 1.0d0
    end do
    do while(coo(2) >= 1.0d0)
       coo(2) = coo(2) - 1.0d0
    end do

    do while(coo(3) < 0.0d0)
       coo(3) = coo(3) + 1.0d0
    end do
    do while(coo(3) >= 1.0d0)
       coo(3) = coo(3) - 1.0d0
    end do

  end subroutine wrap_cooLatt

  !--------------------------------------------------------------------!
  !               volume of the cell spanned by a1,a2,a3               !
  !--------------------------------------------------------------------!

  function cell_volume(avec) result(V)

    implicit none

    double precision, dimension(3,3), intent(in) :: avec
    double precision                             :: V

    V = avec(1,1)*avec(2,2)*avec(3,3) &
      + avec(2,1)*avec(3,2)*avec(1,3) &
      + avec(3,1)*avec(1,2)*avec(2,3) &
      - avec(3,1)*avec(2,2)*avec(1,3) &
      - avec(1,1)*avec(3,2)*avec(2,3) &
      - avec(2,1)*avec(1,2)*avec(3,3)

    V = abs(V)

  end function cell_volume

  !--------------------------------------------------------------------!
  !          largest sphere within a cell spanned by a1,a2,a3          !
  !--------------------------------------------------------------------!

  subroutine inner_sphere(avec, r, c)

    implicit none

    double precision, dimension(3,3),           intent(in)  :: avec
    double precision,                           intent(out) :: r
    double precision, dimension(3),   optional, intent(out) :: c

    double precision, dimension(3) :: v
    double precision               :: h

    ! determine the shortest of the three heights:

    v = vproduct(avec(:,1),avec(:,2))
    v = v/sqrt(sum(v*v))
    h = abs(sum(v*avec(:,3)))
    r = h

    v = vproduct(avec(:,1),avec(:,3))
    v = v/sqrt(sum(v*v))
    h = abs(sum(v*avec(:,2)))
    r = min(r, h)

    v = vproduct(avec(:,2),avec(:,3))
    v = v/sqrt(sum(v*v))
    h = abs(sum(v*avec(:,1)))
    r = min(r, h)

    r = 0.5d0*r

    if (present(c)) c = 0.5d0*(avec(:,1) + avec(:,2) + avec(:,3))

  end subroutine inner_sphere

  !--------------------------------------------------------------------!
  !                         vector/cross product                       !
  !--------------------------------------------------------------------!

  function vproduct(a,b) result(c)

    implicit none

    double precision, dimension(3), intent(in) :: a, b
    double precision, dimension(3)             :: c

    c(1) = a(2)*b(3) - a(3)*b(2)
    c(2) = a(3)*b(1) - a(1)*b(3)
    c(3) = a(1)*b(2) - a(2)*b(1)

  end function vproduct

  !--------------------------------------------------------------------!
  !       smallest supercell to contain a sphere of given radius       !
  !--------------------------------------------------------------------!

  function sphere_supercell(avec, r) result(n)

    implicit none

    double precision, dimension(3, 3), intent(in) :: avec
    double precision,                  intent(in) :: r
    integer,          dimension(3)                :: n

    double precision                 :: V
    double precision, dimension(3,3) :: bvec
    double precision, dimension(3)   :: blen

    double precision, parameter :: EPS = 0.001d0

    V = cell_volume(avec)

    bvec(1,:) = vproduct(avec(2,:), avec(3,:))/V
    bvec(2,:) = -vproduct(avec(1,:), avec(3,:))/V
    bvec(3,:) = vproduct(avec(1,:), avec(2,:))/V

    blen(1) = sqrt(sum(bvec(1,1:3)*bvec(1,1:3)))
    blen(2) = sqrt(sum(bvec(2,1:3)*bvec(2,1:3)))
    blen(3) = sqrt(sum(bvec(3,1:3)*bvec(3,1:3)))

    n(1) = ceiling(r*blen(1) + EPS)
    n(2) = ceiling(r*blen(2) + EPS)
    n(3) = ceiling(r*blen(3) + EPS)

  end function sphere_supercell

  !--------------------------------------------------------------------!
  !        distance of a cell from the home cell at T = [0,0,0]        !
  !--------------------------------------------------------------------!

  function cell_distance_squared(T, avec) result(dist)

    implicit none

    integer,          dimension(3),   intent(in) :: T
    double precision, dimension(3,3), intent(in) :: avec

    double precision :: dist

    double precision, dimension(3) :: dT, v1, v2, normal, h1, h2, d
    double precision               :: d1, d2

    d = [0.0d0, 0.0d0, 0.0d0]

    ! first lattice direction
    v1 = matmul(avec, dble([0, 1, 0]))
    v2 = matmul(avec, dble([0, 0, 1]))
    normal = vproduct(v1, v2)
    normal = normal/sum(normal*normal)
    ! (1) positive segment
    dT = matmul(avec, dble(T - [1, 0, 0]))
    d1 = sum(normal*dT)
    h1 = d1*normal
    ! (2) negative segment
    dT = matmul(avec, dble(T + [1, 0, 0]))
    d2 = -sum(normal*dT)
    h2 = -d2*normal
    if (d1 > 0.0d0) then
       d = h1
    else if (d2 > 0.0d0) then
       d = h2
    else
       d = 0.0d0
    end if

    ! second lattice direction
    v1 = matmul(avec, dble([0, 0, 1]))
    v2 = matmul(avec, dble([1, 0, 0]))
    normal = vproduct(v1, v2)
    normal = normal/sum(normal*normal)
    ! (1) positive segment
    dT = matmul(avec, dble(T - [0, 1, 0]))
    d1 = sum(normal*dT)
    h1 = d1*normal
    ! (2) negative segment
    dT = matmul(avec, dble(T + [0, 1, 0]))
    d2 = -sum(normal*dT)
    h2 = -d2*normal
    if (d1 > 0.0d0) then
       d = d + h1
    else if (d2 > 0.0d0) then
       d = d + h2
    end if

    ! third lattice direction
    v1 = matmul(avec, dble([1, 0, 0]))
    v2 = matmul(avec, dble([0, 1, 0]))
    normal = vproduct(v1, v2)
    normal = normal/sum(normal*normal)
    ! (1) positive segment
    dT = matmul(avec, dble(T - [0, 0, 1]))
    d1 = sum(normal*dT)
    h1 = d1*normal
    ! (2) negative segment
    dT = matmul(avec, dble(T + [0, 0, 1]))
    d2 = -sum(normal*dT)
    h2 = -d2*normal
    if (d1 > 0.0d0) then
       d = d + h1
    else if (d2 > 0.0d0) then
       d = d + h2
    end if

    dist = sum(d*d)

  end function cell_distance_squared


  !--------------------------------------------------------------------!
  !                     check if vector is in list                     !
  !--------------------------------------------------------------------!

  function in_list(vec, list) result(in)

    implicit none

    integer, dimension(:),   intent(in) :: vec
    integer, dimension(:,:), intent(in) :: list
    logical                             :: in

    integer :: d, n, i

    d = size(vec)
    if (d /= size(list(:,1))) then
       write(0,*) "Error: Incompatible vector/list dimensions."
       stop
    end if

    n = size(list(1,:))

    in = .false.
    do i = 1, n
       if (all(vec(:) == list(:,i))) then
          in = .true.
          return
       end if
    end do

  end function in_list


  !--------------------------------------------------------------------!
  !                         Cells within range                         !
  !                                                                    !
  ! The algorithm is as follows:                                       !
  !                                                                    !
  ! First, all T vectors T_i in a sphere with radius Rc around (0,0,0) !
  ! are determined.  Every T vector points to a corner that is shared  !
  ! by 8 cells, so in a second step all vectors                        !
  !                                                                    !
  !               T_i' = T_i + (-1..0,-1..0,-1..0)                     !
  !                                                                    !
  ! are added to the list. It is not sufficient to consider the star   !
  ! around (0,0,0), but instead all corners of the home unit cell have !
  ! to be considered.  Thus, in the final step, all vectors            !
  !                                                                    !
  !               T_i'' = T_i' + (0..1,0..1,0..1)                      !
  !                                                                    !
  ! shifted to the corners of the home unit cell are included.         !
  !--------------------------------------------------------------------!

  subroutine translation_vectors(Rc, avec, nT, T, Tnorm, nc, pbc)

    implicit none

    !------------------------------------------------------------------!
    ! Rc         : cut-off radius                                      !
    ! avec(i,j)  : i-th component of lattice/cell vector j             !
    ! nT         : on entry: dimension of array T --> dimension(3,nT)  !
    !              on exit : number of T vectors found (can be > nT)   !
    !              --> the routine can be used just to calculate the   !
    !                  number of T vectors by setting nT = 0.          !
    ! T(i,j)     : (output) i-th component of the j-th T vector        !
    ! Tnorm(j)   : norm/length of the j-th T vector                    !
    ! nc(i)      : optional input; max. number of cells in direction i !
    !              This is used when setting up the vector star for    !
    !              the cell grid.                                      !
    ! pbc        : optional input; if 'nc' is also specified, only     !
    !              non-wrapping translation vectors (assuming periodic !
    !              boundary conditions) are returned.  This means, it  !
    !              is asserted that no two vectors point to the same   !
    !              cell.
    !------------------------------------------------------------------!

    double precision,                          intent(in)    :: Rc
    double precision, dimension(3,3),          intent(in)    :: avec
    integer,                                   intent(inout) :: nT
    integer,          dimension(3,nT),         intent(out)   :: T
    double precision, dimension(nT), optional, intent(out)   :: Tnorm
    integer,          dimension(3),  optional, intent(in)    :: nc
    logical,                         optional, intent(in)    :: pbc

    double precision,      parameter :: EPS = 1.0d-3

    integer                        :: i1, i2, i3
    double precision               :: Rc2
    double precision               :: cell_dist
    double precision, dimension(3) :: v1, v2
    double precision               :: vnorm
    integer,          dimension(3) :: n
    logical                        :: is_pbc, redundant

    integer                                               :: iT
    integer                                               :: nT1
    integer,          dimension(3)                        :: T_new
    integer,          dimension(:,:), allocatable, target :: T1_tgt1, T1_tgt2
    integer,          dimension(:,:), pointer             :: T1
    integer,          dimension(3)                        :: ic
    integer,          dimension(:,:), allocatable, target :: cell_tgt1, cell_tgt2
    integer,          dimension(:,:), pointer             :: cell
    double precision, dimension(:),   allocatable, target :: Tnorm1_tgt1, Tnorm1_tgt2
    double precision, dimension(:),   pointer             :: Tnorm1

    if (present(pbc)) then
       is_pbc = pbc
    else
       is_pbc = .true.
    end if

    Rc2 = (Rc+EPS)*(Rc+EPS)
    iT  = 0

    ! supercell that can contain a sphere of radius Rc centered at [0,0,0]
    n = sphere_supercell(avec, Rc+EPS)

    ! bounds for T vectors based on that supercell; add +1 in each direction
    ! so that the supercell is large enough from every corner of the home
    ! unit cell
    n(1) = (n(1)/2 + 1) + 1
    n(2) = (n(2)/2 + 1) + 1
    n(3) = (n(3)/2 + 1) + 1

    ! if the number of cells is limited, adjust bounds
    if (present(nc)) then
       n(1) = min(n(1), nc(1) - 1)
       n(2) = min(n(2), nc(2) - 1)
       n(3) = min(n(3), nc(3) - 1)
    end if

    ! allocate initial buffers
    nT1 = max(100, nT)
    allocate(T1_tgt1(3,nT1), Tnorm1_tgt1(nT1), cell_tgt1(3,nT1))
    T1 => T1_tgt1
    Tnorm1 => Tnorm1_tgt1
    cell => cell_tgt1

    redundant = .false.

    ! loop over all T vectors in that supercell
    loop3 : do i3 = -n(3), n(3)
    loop2 : do i2 = -n(2), n(2)
    loop1 : do i1 = -n(1), n(1)
       T_new(1:3) = [i1, i2, i3]
       if (all(T_new == 0)) cycle loop1
       if (present(nc) .and. is_pbc) then
          ! check if cell vector points to same cell as other vector
          ic(1:3) = cell_from_Cvec([1,1,1], T_new, nc, is_pbc)
          redundant = in_list(ic, cell(1:3,1:iT))
       end if
       ! cell_dist = cell_distance_squared(T_new, avec)
       ! cutoff : if (cell_dist < Rc2) then
          ! only consider T_new if the vector is not stored already
          if (.not. (redundant .or. in_list(T_new, T1(1:3,1:iT)))) then
             iT = iT + 1
             ! increase buffer size, if necessary
             if (iT > nT1) then
                nT1 = 2*nT1
                if (allocated(T1_tgt1)) then
                   allocate(T1_tgt2(3,nT1), Tnorm1_tgt2(nT1), cell_tgt2(3,nT1))
                   T1_tgt2(1:3,1:iT-1) = T1_tgt1(1:3,1:iT-1)
                   Tnorm1_tgt2(1:iT-1) = Tnorm1_tgt1(1:iT-1)
                   cell_tgt2(1:3,1:iT-1) = cell_tgt1(1:3,1:iT-1)
                   T1 => T1_tgt2
                   Tnorm1 => Tnorm1_tgt2
                   cell => cell_tgt2
                   deallocate(T1_tgt1, Tnorm1_tgt1, cell_tgt1)
                else if (allocated(T1_tgt2)) then
                   allocate(T1_tgt1(3,nT1), Tnorm1_tgt1(nT1), cell_tgt1(3,nT1))
                   T1_tgt1(1:3,1:iT-1) = T1_tgt2(1:3,1:iT-1)
                   Tnorm1_tgt1(1:iT-1) = Tnorm1_tgt2(1:iT-1)
                   cell_tgt1(1:3,1:iT-1) = cell_tgt2(1:3,1:iT-1)
                   T1 => T1_tgt1
                   Tnorm1 => Tnorm1_tgt1
                   cell => cell_tgt1
                   deallocate(T1_tgt2, Tnorm1_tgt2, cell_tgt2)
                else
                   write(0,*) "Error: unable to allocate memory."
                   stop
                end if
             end if
             ! store T vector and its norm
             T1(1:3,iT) = T_new(1:3)
             v2 = matmul(avec, dble(T_new))
             Tnorm1(iT) = sqrt(sum(v2*v2))
             if (present(nc)) then
                cell(1:3,iT) = ic(1:3)
             end if
          end if
       ! end if cutoff
    end do loop1
    end do loop2
    end do loop3

    if (iT <= nT) then
       T(:,1:iT) = T1(:,1:iT)
       if (present(Tnorm)) Tnorm(1:iT) = Tnorm1(1:iT)
    end if
    nT = iT

    if (allocated(T1_tgt1)) deallocate(T1_tgt1, Tnorm1_tgt1, cell_tgt1)
    if (allocated(T1_tgt2)) deallocate(T1_tgt2, Tnorm1_tgt2, cell_tgt2)

  end subroutine translation_vectors

end module lclist
