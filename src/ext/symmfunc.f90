!-----------------------------------------------------------------------
!                         symmetry functions
!           --> translation/rotation invariant coordinates
!
! This module implements the symmetry function basis by Behler following
! the original publication: J. Behler, J. Chem. Phys. 134 (2011) 074106.
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
! 2011-04-25 Alexander Urban (AU)
! 2011-07-20 AU --- prefactors for angular functions added
! 2013-05-07 AU --- earlier evaluation of these prefactors
! 2013-09-30 AU --- correct wrong behavior of angular functions
!-----------------------------------------------------------------------

module symmfunc

  implicit none

  double precision, parameter, private :: PI    = 3.14159265358979d0
  double precision, parameter, private :: PI2   = 2.0d0*PI
  double precision, parameter, private :: EPS   = 1.0d-12
  integer,          parameter, private :: NG    = 5

  !--------------------------------------------------------------------!
  ! sf_nG_type(i)    : number of symmetry functions for species i      !
  !                                                                    !
  ! sf_nTypes        : number of different atom/point types            !
  ! sf_nPairs        : number type pairs (= nTypes!)                   !
  ! sf_nGmax         : max. num. symmetry functs. per type pair/triple !
  ! sf_idx2(i,j)     : pair index for species i-j                      !
  ! sf_idx3(i,j,k)   : triple index for species i-j-k                  !
  ! sf_nG(i,pij)     : num. symmetry function type i for type pair pij !
  ! sf_pG?(i,iG,pij) : parameter i for the iG-th symm.-func. of type ? !
  !                    for type pair pij                               !
  !                    Parameters:  G1:  Rc                            !
  !                                 G2:  Rc, Rs, eta                   !
  !                                 G3:  Rc, kappa                     !
  !                                 G4:  Rc, lambda, zeta, eta         !
  !                                 G5:  Rc, lambda, zeta, eta         !
  !--------------------------------------------------------------------!

  integer,          dimension(:),     allocatable, public  :: sf_nG_type

  integer,                                         private :: sf_nTypes
  integer,                                         private :: sf_nPairs
  integer,                                         private :: sf_nTriples
  integer,                                         private :: sf_nGMax
  integer,          dimension(:,:),   allocatable, private :: sf_idx2
  integer,          dimension(:,:,:), allocatable, private :: sf_idx3
  integer,          dimension(:),     allocatable, private :: sf_iG02
  integer,          dimension(:),     allocatable, private :: sf_iG03
  integer,          dimension(:,:),   allocatable, private :: sf_nG2
  integer,          dimension(:,:),   allocatable, private :: sf_nG3
  double precision, dimension(:,:),   allocatable, private :: sf_pG1
  double precision, dimension(:,:,:), allocatable, private :: sf_pG2
  double precision, dimension(:,:,:), allocatable, private :: sf_pG3
  double precision, dimension(:,:,:), allocatable, private :: sf_pG4
  double precision, dimension(:,:,:), allocatable, private :: sf_pG5

  !--------------------------------------------------------------------!

  integer, public :: stdout = 6
  integer, public :: stderr = 0

  !--------------------------------------------------------------------!

  logical, private :: isInit  = .false.
  integer, private :: memSize = 0

contains !-------------------------------------------------------------!

  subroutine sf_init(ntypes, nGmax)

    implicit none

    integer, intent(in) :: ntypes
    integer, intent(in) :: nGmax

    integer :: i, j, k, idx2, idx3

    sf_nTypes = ntypes
    sf_nGMax  = nGmax

    ! initialize type pair/triple index:
    ! note: for angular symmetry functions A-B-C = A-C-B
    allocate(sf_idx2(ntypes,ntypes), sf_idx3(ntypes,ntypes,ntypes))
    idx2 = 0
    idx3 = 0
    do i = 1, ntypes
       do j = 1, ntypes
          idx2 = idx2 + 1
          sf_idx2(i,j) = idx2
          idx3 = idx3 + 1
          sf_idx3(i,j,j) = idx3
          do k = j+1, ntypes
             idx3 = idx3 + 1
             sf_idx3(i,j,k) = idx3
             sf_idx3(i,k,j) = idx3
          end do
       end do
    end do
    sf_nPairs = idx2
    sf_nTriples = idx3

    allocate(sf_nG_type(ntypes),            &
             sf_nG2(NG, sf_nPairs),         &
             sf_nG3(NG, sf_nTriples),       &
             sf_iG02(sf_nPairs),            &
             sf_iG03(sf_nTriples),          &
             sf_pG1(nGmax, sf_nPairs),      &
             sf_pG2(3, nGmax, sf_nPairs),   &
             sf_pG3(2, nGmax, sf_nPairs),   &
             sf_pG4(4, nGmax, sf_nTriples), &
             sf_pG5(4, nGmax, sf_nTriples)  )
    sf_nG_type(:) = 0
    sf_nG2(:,:) = 0
    sf_nG3(:,:) = 0

    isInit = .true.

  end subroutine sf_init

  !--------------------------------------------------------------------!

  subroutine sf_final()

    implicit none

    if (isInit) then
       deallocate(sf_idx2, sf_idx3, sf_nG2, sf_nG3, sf_pG1, sf_pG2, &
                  sf_pG3,  sf_pG4, sf_pG5, sf_nG_type, sf_iG02, sf_iG03)
       isInit = .false.
    end if

  end subroutine sf_final

  !--------------------------------------------------------------------!
  !                      parameter initialization                      !
  !--------------------------------------------------------------------!

  subroutine sf_add_rad(funct, type1, type2, Rc, Rs, eta, kappa)

    implicit none

    integer,                    intent(in) :: funct
    integer,                    intent(in) :: type1, type2
    double precision, optional, intent(in) :: Rc, Rs, eta, kappa

    integer :: iG, pij

    pij = sf_idx2(type1,type2)

    if (sf_nG2(funct,pij) >= sf_nGMax) then
       write(stderr,*) "Error: Memory exceeded in `sf_add()'."
       write(stderr,*) sf_nG2(funct,pij), sf_nGMax
       stop
    end if

    if (.not. present(Rc)) then
       write(stderr,*) "Error: Rc undefined in `sf_add()'"
       stop
    end if

    iG = sf_nG2(funct,pij) + 1
    select case(funct)
    case(1)
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG2(funct,pij) = iG
       sf_pG1(iG,pij)   = Rc
    case(2)
       if (.not. (present(Rs) .and. present(eta))) then
          write(stderr,*) "Error: needed for G2: Rs, eta"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG2(funct,pij) = iG
       sf_pG2(1,iG,pij) = Rc
       sf_pG2(2,iG,pij) = Rs
       sf_pG2(3,iG,pij) = eta
    case(3)
       if (.not. present(kappa)) then
          write(stderr,*) "Error: needed for G3: kappa"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG2(funct,pij) = iG
       sf_pG3(1,iG,pij) = Rc
       sf_pG3(2,iG,pij) = kappa
    case default
       write(stderr,*) "Error: Not a radial symmetry function type: ", funct
       stop
    end select

    call sf_build_index()

  end subroutine sf_add_rad

  !--------------------------------------------------------------------!

  subroutine sf_add_ang(funct, type1, type2, type3, Rc, eta, lambda, zeta)

    implicit none

    integer,                    intent(in) :: funct
    integer,                    intent(in) :: type1, type2, type3
    double precision, optional, intent(in) :: Rc, eta, lambda, zeta

    integer :: iG, tijk

    tijk = sf_idx3(type1, type2, type3)

    if (sf_nG3(funct,tijk) >= sf_nGMax) then
       write(stderr,*) "Error: Memory exceeded in `sf_add()'."
       write(stderr,*) sf_nG3(funct,tijk), sf_nGMax
       stop
    end if

    if (.not. present(Rc)) then
       write(stderr,*) "Error: Rc undefined in `sf_add()'"
       stop
    end if

    iG = sf_nG3(funct,tijk) + 1
    select case(funct)
    case(4)
       if (.not. (present(lambda) .and. present(zeta) &
           .and. present(eta))) then
          write(stderr,*) "Error: needed for G4: eta, lambda, zeta"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG3(funct,tijk) = iG
       sf_pG4(1,iG,tijk) = Rc
       sf_pG4(2,iG,tijk) = lambda
       sf_pG4(3,iG,tijk) = zeta
       sf_pG4(4,iG,tijk) = eta
    case(5)
       if (.not. (present(lambda) .and. present(zeta) &
           .and. present(eta))) then
          write(stderr,*) "Error: needed for G5: eta, lambda, zeta"
          stop
       end if
       sf_nG_type(type1) = sf_nG_type(type1) + 1
       sf_nG3(funct,tijk) = iG
       sf_pG5(1,iG,tijk) = Rc
       sf_pG5(2,iG,tijk) = lambda
       sf_pG5(3,iG,tijk) = zeta
       sf_pG5(4,iG,tijk) = eta
    case default
       write(stderr,*) "Error: Not an angular symmetry function type: ", funct
       stop
    end select

    call sf_build_index()

  end subroutine sf_add_ang

  !--------------------------------------------------------------------!
  ! The functions are organized in memory not in the same order as     !
  ! they are added, but rather depending on their type and species.    !
  ! The procedure `sf_build_index()' generates the indices for radial  !
  ! and angular functions.                                             !
  !--------------------------------------------------------------------!

  subroutine sf_build_index()

    implicit none

    integer :: pij, tijk, i
    integer :: itype1, itype2, itype3

    do itype1 = 1, sf_nTypes
       i = 0
       do itype2 = 1, sf_nTypes
          pij = sf_idx2(itype1, itype2)
          sf_iG02(pij) = i
          i = i + sum(sf_nG2(:,pij))
          do itype3 = itype2, sf_nTypes
             tijk = sf_idx3(itype1, itype2, itype3)
             sf_iG03(tijk) = i
             i = i + sum(sf_nG3(:,tijk))
          end do
       end do
    end do

  end subroutine sf_build_index

  !--------------------------------------------------------------------!
  !               evaluate finger print for single atom                !
  !--------------------------------------------------------------------!

  subroutine sf_fingerprint(type1, Xi, nx, X, typej, n, G, dGi, dGj)

    implicit none

    !------------------------------------------------------------------!
    ! type1     : type index of the central particle                   !
    ! Xi(1:3)   : cartesian coordinates of the central particle        !
    ! nx        : number of other particles                            !
    ! X(1:3,j)  : cartesian coordinates of the j-th particle           !
    ! typej(j)  : atom type index of j-th atom                         !
    ! n         : number of symmetry functions                         !
    ! G(j)      : (out) value of the j-th symmetry function            !
    ! dGi(1:3,j): (out, optional) derivative of the j-th symmetry      !
    !             function wrt. the central particle's coordinates     !
    ! dGj(:,j,k): (out, optional) derivative of the j-th symmetry      !
    !             function wrt. the coordinates of the k-th particle   !
    ! Note: dGi and dGj have to be present/absent simultaneously       !
    !------------------------------------------------------------------!

    integer,                                       intent(in)  :: type1
    double precision, dimension(3),                intent(in)  :: Xi
    integer,                                       intent(in)  :: nx
    double precision, dimension(3,nx),             intent(in)  :: X
    integer,          dimension(nx),               intent(in)  :: typej
    integer,                                       intent(in)  :: n
    double precision, dimension(n),                intent(out) :: G
    double precision, dimension(3,n),    optional, intent(out) :: dGi
    double precision, dimension(3,n,nx), optional, intent(out) :: dGj

    integer          :: j, k, pij, tijk
    integer          :: iG_ang0, iG
    double precision :: Rij, Rik, Rjk, cost
    logical          :: deriv

    double precision, dimension(3) :: R1, R2, R3

    if (present(dGi) .and. present(dGj)) then
       deriv = .true.
    else
       deriv = .false.
    end if

    if (.not. isInit) return

    if (n < sf_nG_type(type1)) then
       write(stderr,*) "Error: number of symmetry functions exceeded."
       stop
    end if

    G(1:n)  = 0.0d0
    iG_ang0 = 0
    if (deriv) then
       dGi(1:3,1:n) = 0.0d0
       dGj(1:3,1:n,1:nx) = 0.0d0
    end if
    jloop : do j = 1, nx
       pij = sf_idx2(type1, typej(j))

       R1  = X(:,j) - Xi(:)
       Rij = sqrt(sum(R1*R1))
       R1  = R1/Rij

       iG = sf_iG02(pij)

       ! radial symmetry functions:
       if (deriv) then
          call sf_G1_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j))
          call sf_G2_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j))
          call sf_G3_update(pij, R1(1:3), Rij, n, G, iG, dGi, dGj(1:3,1:n,j))
       else
          call sf_G1_update(pij, R1(1:3), Rij, n, G, iG)
          call sf_G2_update(pij, R1(1:3), Rij, n, G, iG)
          call sf_G3_update(pij, R1(1:3), Rij, n, G, iG)
       end if

       ! if (iG_ang0 == 0) iG_ang0 = sum(sf_nG2(:,:))
       kloop : do k = j+1, nx
          tijk = sf_idx3(type1, typej(j), typej(k))

          R2  = X(:,k) - Xi(:)
          Rik = sqrt(sum(R2*R2))
          if (Rik <= 1.0d-8) then
             write(stderr, *) "Warning: an atom in the neighbor list is " &
                  // "identical to the central atom (sf_fingerprint)"
          else
             R2  = R2/Rik
          end if

          R3  = X(:,k) - X(:,j)
          Rjk = sqrt(sum(R3*R3))
          if (Rjk <= 1.0d-8) then
             write(stderr, *) "Warning: redundant atoms in neighbor " &
                              // "list (sf_fingerprint)"
          else
             R3  = R3/Rjk
          end if

          ! cos(theta_ijk)
          cost = sum(R1*R2)
          cost = max(cost,-1.0d0)
          cost = min(cost, 1.0d0)

          iG = sf_iG03(tijk)

          ! angular symmetry functions:
          if (deriv) then
             call sf_G4_update(tijk, R1(1:3), R2(1:3), R3(1:3), Rij, Rik, Rjk, cost, &
                               n, G, iG, dGi, dGj(1:3,1:n,j), dGj(1:3,1:n,k))
             call sf_G5_update(tijk, R1(1:3), R2(1:3), Rij, Rik, cost, &
                               n, G, iG, dGi, dGj(1:3,1:n,j), dGj(1:3,1:n,k))
          else
             call sf_G4_update(tijk, R1(1:3), R2(1:3), R3(1:3), Rij, Rik, Rjk, &
                               cost, n, G, iG)
             call sf_G5_update(tijk, R1(1:3), R2(1:3), Rij, Rik, cost, &
                               n, G, iG)
          end if

       end do kloop
    end do jloop

  end subroutine sf_fingerprint



  !==================== radial symmetry functions =====================!



  !--------------------------------------------------------------------!
  !            radial symm. function 1 (eq. 5 in Ref. [1])             !
  !                                                                    !
  ! here: vecRij = vec{R}_ij/||vec{R})_ij||                            !
  !--------------------------------------------------------------------!

  subroutine sf_G1_ij(Rij, Rc, Gij, dGij)

    implicit none

    double precision, intent(in)  :: Rij
    double precision, intent(in)  :: Rc
    double precision, intent(out) :: Gij
    double precision, intent(out) :: dGij

    double precision :: fc, dfc

    call sf_cut(Rij, Rc, fc, dfc)

    Gij  = fc
    dGij = dfc

  end subroutine sf_G1_ij

  !--------------------------------------------------------------------!

  subroutine sf_G1_update(pij, vecRij, Rij, n, G, iG, dGi, dGj)

    implicit none

    integer,                                    intent(in)    :: pij
    double precision, dimension(3),             intent(in)    :: vecRij
    double precision,                           intent(in)    :: Rij
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj

    integer          :: iG1
    double precision :: G1, dG1
    double precision :: Rc

    do iG1 = 1, sf_nG2(1,pij)
       iG    = iG + 1
       Rc    = sf_pG1(iG1,pij)
       call sf_G1_ij(Rij, Rc, G1, dG1)
       G(iG) = G(iG) + G1
       if (present(dGi)) dGi(1:3,iG) = dGi(1:3,iG) - dG1*vecRij(1:3)
       if (present(dGj)) dGj(1:3,iG) = dGj(1:3,iG) + dG1*vecRij(1:3)
    end do

  end subroutine sf_G1_update

  !--------------------------------------------------------------------!
  !            radial symm. function 2 (eq. 6 in Ref. [1])             !
  !--------------------------------------------------------------------!

  subroutine sf_G2_ij(Rij, Rc, Rs, eta, Gij, dGij)

    implicit none

    double precision, intent(in)  :: Rij
    double precision, intent(in)  :: Rc
    double precision, intent(in)  :: Rs
    double precision, intent(in)  :: eta
    double precision, intent(out) :: Gij
    double precision, intent(out) :: dGij

    double precision :: fc, dfc
    double precision :: fexp, arg

    call sf_cut(Rij, Rc, fc, dfc)

    arg  = Rij - Rs
    fexp = exp(-eta*arg*arg)

    Gij  = fexp*fc
    dGij = fexp*( dfc - 2.0d0*eta*arg*fc  )

  end subroutine sf_G2_ij

  !--------------------------------------------------------------------!

  subroutine sf_G2_update(pij, vecRij, Rij, n, G, iG, dGi, dGj)

    implicit none

    integer,                                    intent(in)    :: pij
    double precision, dimension(3),             intent(in)    :: vecRij
    double precision,                           intent(in)    :: Rij
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj

    integer          :: iG2
    double precision :: Rc, Rs, eta
    double precision :: G2, dG2

    do iG2 = 1, sf_nG2(2,pij)
       iG    = iG + 1
       Rc    = sf_pG2(1,iG2,pij)
       Rs    = sf_pG2(2,iG2,pij)
       eta   = sf_pG2(3,iG2,pij)
       call sf_G2_ij(Rij, Rc, Rs, eta, G2, dG2)
       G(iG) = G(iG) + G2
       if (present(dGi)) dGi(1:3,iG) = dGi(1:3,iG) - dG2*vecRij(1:3)
       if (present(dGj)) dGj(1:3,iG) = dGj(1:3,iG) + dG2*vecRij(1:3)
    end do

  end subroutine sf_G2_update

  !--------------------------------------------------------------------!
  !            radial symm. function 3 (eq. 6 in Ref. [1])             !
  !--------------------------------------------------------------------!

  subroutine sf_G3_ij(Rij, Rc, kappa, Gij, dGij)

    implicit none

    double precision, intent(in)  :: Rij
    double precision, intent(in)  :: Rc
    double precision, intent(in)  :: kappa
    double precision, intent(out) :: Gij
    double precision, intent(out) :: dGij

    double precision :: fc, dfc
    double precision :: fcos, fsin, arg

    call sf_cut(Rij, Rc, fc, dfc)

    arg  = kappa*Rij
    fcos = cos(arg)
    fsin = sin(arg)

    Gij  = fcos*fc
    dGij = fcos*dfc - kappa*fsin*fc

  end subroutine sf_G3_ij

  !--------------------------------------------------------------------!

  subroutine sf_G3_update(pij, vecRij, Rij, n, G, iG, dGi, dGj)

    implicit none

    integer,                                    intent(in)    :: pij
    double precision, dimension(3),             intent(in)    :: vecRij
    double precision,                           intent(in)    :: Rij
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj

    integer          :: iG3
    double precision :: Rc, kappa
    double precision :: G3, dG3

    do iG3 = 1, sf_nG2(3,pij)
       iG    = iG + 1
       Rc    = sf_pG3(1,iG3,pij)
       kappa = sf_pG3(2,iG3,pij)
       call sf_G3_ij(Rij, Rc, kappa, G3, dG3)
       G(iG) = G(iG) + G3
       if (present(dGi)) dGi(1:3,iG) = dGi(1:3,iG) - dG3*vecRij(1:3)
       if (present(dGj)) dGj(1:3,iG) = dGj(1:3,iG) + dG3*vecRij(1:3)
    end do

  end subroutine sf_G3_update



  !==================== angular symmetry functions ====================!

  !--------------------------------------------------------------------!
  !          first angular symm. function (eq. 8 in Ref. [1])          !
  !--------------------------------------------------------------------!

  !--------------------------------------------------------------------!
  !                        angle dependend part                        !
  ! This is function F_1(R_ij,R_ik) in the documentation.              !
  !--------------------------------------------------------------------!

  subroutine sf_F1_ijk(vecRij, vecRik, Rij, Rik, cost, lambda, &
                       zeta, Fijk, dFijk_dRj, dFijk_dRk)

    implicit none

    double precision, dimension(3), intent(in)  :: vecRij, vecRik
    double precision,               intent(in)  :: Rij, Rik
    double precision,               intent(in)  :: cost
    double precision,               intent(in)  :: lambda
    double precision,               intent(in)  :: zeta
    double precision,               intent(out) :: Fijk
    double precision, dimension(3), intent(out) :: dFijk_dRj
    double precision, dimension(3), intent(out) :: dFijk_dRk

    double precision :: arg, prefactor

    if (abs(cost) > 1.0d0) write(*,*) "cos(theta) = ", cost

    arg = 0.5d0*(1.0d0 + lambda*cost)
    prefactor = 0.5d0*zeta*lambda*arg**(zeta-1.0d0)

    Fijk  = arg**zeta

    dFijk_dRj = prefactor*( -cost*( vecRij/Rij ) &
                            + vecRik/Rij )

    dFijk_dRk = prefactor*( -cost*( vecRik/Rik ) &
                            + vecRij/Rik )

  end subroutine sf_F1_ijk

  !--------------------------------------------------------------------!
  !                      distance dependend part                       !
  ! This is function F_2(R) in the documentation.                      !
  ! --> very similar to sf_G2_ij(), may be combined in future.         !
  !--------------------------------------------------------------------!

  subroutine sf_F2_ij(Rij, Rc, eta, Fij, dFij)

    implicit none

    double precision, intent(in)  :: Rij
    double precision, intent(in)  :: Rc
    double precision, intent(in)  :: eta
    double precision, intent(out) :: Fij
    double precision, intent(out) :: dFij

    double precision :: fc, dfc
    double precision :: fexp

    call sf_cut(Rij, Rc, fc, dfc)

    fexp = exp(-eta*Rij*Rij)

    Fij  = fexp*fc
    dFij = fexp*( dfc - 2.0d0*eta*Rij*fc  )

  end subroutine sf_F2_ij

  !--------------------------------------------------------------------!

  subroutine sf_G4_update(tijk, vecRij, vecRik, vecRjk, Rij, Rik, Rjk, &
                          cost, n, G, iG, dGi, dGj, dGk)

    implicit none

    integer,                                    intent(in)    :: tijk
    double precision, dimension(3),             intent(in)    :: vecRij, vecRik, vecRjk
    double precision,                           intent(in)    :: Rij, Rik, Rjk
    double precision,                           intent(in)    :: cost
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj
    double precision, dimension(3,n), optional, intent(inout) :: dGk

    integer                        :: iG4
    double precision               :: Rc, lambda, zeta, eta
    double precision               :: G4
    double precision, dimension(3) :: dG4, dF1j, dF1k

    double precision :: F1, F2ij, dF2ij, F2ik, dF2ik, F2jk, dF2jk

    do iG4 = 1, sf_nG3(4,tijk)
       iG     = iG + 1
       Rc     = sf_pG4(1,iG4,tijk)
       lambda = sf_pG4(2,iG4,tijk)
       zeta   = sf_pG4(3,iG4,tijk)
       eta    = sf_pG4(4,iG4,tijk)

       call sf_F1_ijk(vecRij, vecRik, Rij, Rik, cost, lambda, &
                      zeta, F1, dF1j, dF1k)
       call sf_F2_ij(Rij, Rc, eta, F2ij, dF2ij)
       call sf_F2_ij(Rik, Rc, eta, F2ik, dF2ik)
       call sf_F2_ij(Rjk, Rc, eta, F2jk, dF2jk)

       G4    = F1*F2ij*F2ik*F2jk
       ! factor of 2 for k<j in sf_fingerprint()
       G(iG) = G(iG) + 2.0d0*G4

       if (present(dGi)) then
          dG4(1:3) = F2jk*( -(dF1j(1:3)+dF1k(1:3))*F2ij*F2ik &
                            - vecRij(1:3)*F1*dF2ij*F2ik      &
                            - vecRik(1:3)*F1*F2ij*dF2ik )
          dGi(1:3,iG) = dGi(1:3,iG) + 2.0d0*dG4(1:3)
       end if
       if (present(dGj)) then
          dG4(1:3) =  dF1j(1:3)*F2ij*F2ik*F2jk           &
                      + vecRij(1:3)*F1*dF2ij*F2ik*F2jk   &
                      - vecRjk(1:3)*F1*F2ij*F2ik*dF2jk
          dGj(1:3,iG) = dGj(1:3,iG) + 2.0d0*dG4(1:3)
       end if
       if (present(dGk)) then
          dG4(1:3) =  dF1k(1:3)*F2ij*F2ik*F2jk           &
                      + vecRik(1:3)*F1*F2ij*dF2ik*F2jk   &
                      + vecRjk(1:3)*F1*F2ij*F2ik*dF2jk
          dGk(1:3,iG) = dGk(1:3,iG) + 2.0d0*dG4(1:3)
       end if

    end do

  end subroutine sf_G4_update

  !--------------------------------------------------------------------!
  !         second angular symm. function (eq. 9 in Ref. [1])          !
  !--------------------------------------------------------------------!

  subroutine sf_G5_update(tijk, vecRij, vecRik, Rij, Rik, cost, n, &
                          G, iG, dGi, dGj, dGk)

    implicit none

    integer,                                    intent(in)    :: tijk
    double precision, dimension(3),             intent(in)    :: vecRij, vecRik
    double precision,                           intent(in)    :: Rij, Rik
    double precision,                           intent(in)    :: cost
    integer,                                    intent(in)    :: n
    double precision, dimension(n),             intent(inout) :: G
    integer,                                    intent(inout) :: iG
    double precision, dimension(3,n), optional, intent(inout) :: dGi
    double precision, dimension(3,n), optional, intent(inout) :: dGj
    double precision, dimension(3,n), optional, intent(inout) :: dGk

    integer                        :: iG5
    double precision               :: Rc, lambda, zeta, eta
    double precision               :: G5
    double precision, dimension(3) :: dG5, dF1j, dF1k

    double precision :: F1, F2ij, dF2ij, F2ik, dF2ik

    do iG5 = 1, sf_nG3(5,tijk)
       iG     = iG + 1
       Rc     = sf_pG5(1,iG5,tijk)
       lambda = sf_pG5(2,iG5,tijk)
       zeta   = sf_pG5(3,iG5,tijk)
       eta    = sf_pG5(4,iG5,tijk)

       call sf_F1_ijk(vecRij, vecRik, Rij, Rik, cost, lambda, &
                      zeta, F1, dF1j, dF1k)
       call sf_F2_ij(Rij, Rc, eta, F2ij, dF2ij)
       call sf_F2_ij(Rik, Rc, eta, F2ik, dF2ik)

       G5    = F1*F2ij*F2ik
       ! factor of two for k<j in sf_fingerprint()
       G(iG) = G(iG) + 2.0d0*G5

       if (present(dGi)) then
          dG5(1:3) = -(dF1j(1:3)+dF1k(1:3))*F2ij*F2ik        &
                   - vecRij(1:3)*F1*dF2ij*F2ik &
                   - vecRik(1:3)*F1*F2ij*dF2ik
          dGi(1:3,iG) = dGi(1:3,iG) + 2.0d0*dG5(1:3)
       end if
       if (present(dGj)) then
          dG5(1:3) = dF1j(1:3)*F2ij*F2ik        &
                   + vecRij(1:3)*F1*dF2ij*F2ik
          dGj(1:3,iG) = dGj(1:3,iG) + 2.0d0*dG5(1:3)
       end if
       if (present(dGk)) then
          dG5(1:3) = dF1k(1:3)*F2ij*F2ik        &
                   + vecRik(1:3)*F1*F2ij*dF2ik
          dGk(1:3,iG) = dGk(1:3,iG) + 2.0d0*dG5(1:3)
       end if

    end do

  end subroutine sf_G5_update


  !======================= auxiliary functions ========================!



  !--------------------------------------------------------------------!
  !                      cosine cut-off function                       !
  !                        (eq. 4 in Ref. [1])                         !
  !--------------------------------------------------------------------!

  subroutine sf_cut(Rij, Rc, fc, dfc)

    implicit none

    double precision, intent(in)  :: Rij, Rc
    double precision, intent(out) :: fc, dfc

    if (Rij >= Rc) then
       fc  = 0.0d0
       dfc = 0.0d0
    else
       fc  =  0.5d0*(cos(PI*Rij/Rc) + 1.0d0)
       dfc = -0.5d0*PI/Rc*sin(PI*Rij/Rc)
    end if

  end subroutine sf_cut

end module symmfunc
