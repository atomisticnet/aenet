!-----------------------------------------------------------------------
! aenet.f90 -- public interface to the aenet routines (aenetLib)
!
! This module also provides the C bindings for libaenet.so/.a .
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
! 2014-08-31 Nongnuch Artrith, Alexander Urban
!-----------------------------------------------------------------------
module aenet

  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_char, c_bool, &
                                         c_ptr, c_f_pointer, c_null_char

  use aeio,        only: TYPELEN, PATHLEN

  use feedforward, only: ff_eval, ff_deriv

  use geometry,    only: geo_type_conv,       &
                         geo_recip_lattice

  use io,          only: io_cstring_len,      &
                         io_cstring2f

  use lclist,      only: lcl_nmax_nbdist,     &
                         lcl_init,            &
                         lcl_final,           &
                         lcl_nbdist_cart,     &
                         lcl_print_info

  use potential,   only: NNPot,               &
                         load_NNPot,          &
                         del_NNPot,           &
                         pot_init,            &
                         pot_final,           &
                         pot_print_info,      &
                         pot_get_range

  use sfbasis,     only: FingerprintBasis,       &
                         new_SFBasis,            &
                         del_SFBasis,            &
                         save_SFBasis,           &
                         load_SFBasis,           &
                         save_SFBasis_ASCII,     &
                         load_SFBasis_ASCII,     &
                         sfb_print_info,         &
                         sfb_set_typeid,         &
                         sfb_set_typespin,       &
                         sfb_eval,               &
                         sfb_reconstruct_radial

  use sfsetup,     only: stp_init,            &
                         stp_final,           &
                         stp_nsf_max,         &
                         stp_eval

  use timing,      only: tng_timing, tng_timing2

  implicit none
  private
  save

  public :: aenet_init,                     &
            aenet_final,                    &
            aenet_all_loaded,               &
            aenet_atomic_energy,            &
            aenet_atomic_energy_and_forces, &
            aenet_convert_atom_types,       &
            aenet_free_atom_energy,         &
            aenet_load_potential,           &
            aenet_print_info,               &
            aenet_nbl_init,                 &
            aenet_nbl_final,                &
            aenet_nbl_neighbors,            &
            aenet_sfb_init,                 &
            aenet_sfb_final,                &
            aenet_sfb_eval,                 &
            aenet_sfb_reconstruct_radial

  !---------------------------- constants -----------------------------!

  ! return status
  integer(kind=c_int), bind(C, name='AENET_OK'),         public :: AENET_OK = 0_c_int
  integer(kind=c_int), bind(C, name='AENET_ERR_INIT'),   public :: AENET_ERR_INIT = 1_c_int
  integer(kind=c_int), bind(C, name='AENET_ERR_MALLOC'), public :: AENET_ERR_MALLOC = 2_c_int
  integer(kind=c_int), bind(C, name='AENET_ERR_IO'),     public :: AENET_ERR_IO = 3_c_int
  integer(kind=c_int), bind(C, name='AENET_ERR_TYPE'),   public :: AENET_ERR_TYPE = 4_c_int

  integer(kind=c_int), bind(C, name='AENET_TYPELEN'),    public :: AENET_TYPELEN = TYPELEN
  integer(kind=c_int), bind(C, name='AENET_PATHLEN'),    public :: AENET_PATHLEN = PATHLEN

  logical(kind=c_bool), bind(C, name='AENET_TRUE'),      public :: AENET_TRUE = .true.
  logical(kind=c_bool), bind(C, name='AENET_FALSE'),     public :: AENET_FALSE = .false.

  !---------------------------- variables -----------------------------!

  integer(kind=c_int), bind(C, name='aenet_nsf_max'), public :: aenet_nsf_max
  integer(kind=c_int), bind(C, name='aenet_nnb_max'), public :: aenet_nnb_max
  real(kind=c_double), bind(C, name='aenet_Rc_min'),  public :: aenet_Rc_min
  real(kind=c_double), bind(C, name='aenet_Rc_max'),  public :: aenet_Rc_max

  !----------------------------- private ------------------------------!

  logical, private :: aenet_is_init = .false.
  logical, private :: aenet_is_loaded = .false.

  integer,                                             private :: aenet_ntypes
  integer,                                             private :: aenet_nvalues_max
  integer,                                             private :: aenet_nweights_max
  character(len=TYPELEN), dimension(:),   allocatable, private :: aenet_atom_types
  type(NNPot),            dimension(:),   allocatable, private :: aenet_pot
  type(FingerprintBasis),                              private :: aenet_sfb

  real(kind=c_double),    dimension(:,:), allocatable, target, private :: aenet_coo_latt

contains

  !--------------------------------------------------------------------!
  !                  Initialization and Finalization                   !
  !--------------------------------------------------------------------!

  subroutine aenet_init(atom_types, stat)

    implicit none

    character(len=*), dimension(:), intent(in)  :: atom_types
    integer,                        intent(out) :: stat

    integer :: ok

    stat = AENET_OK
    if (aenet_is_init) then
       stat = AENET_ERR_INIT
       return
    end if

    aenet_ntypes = size(atom_types)
    allocate(aenet_pot(aenet_ntypes),        &
             aenet_atom_types(aenet_ntypes), &
             stat=ok)
    if (ok /= 0) then
       aenet_ntypes = 0
       stat = AENET_ERR_MALLOC
       return
    end if
    aenet_atom_types = atom_types
    aenet_nvalues_max = 0
    aenet_nweights_max = 0
    aenet_is_init = .true.

  end subroutine aenet_init

  subroutine aenet_init_C(ntypes, atom_types, stat) bind(C, name='aenet_init')

    implicit none

    integer(kind=c_int), value,             intent(in)  :: ntypes
    type(c_ptr), dimension(ntypes), target, intent(in)  :: atom_types
    integer(kind=c_int),                    intent(out) :: stat

    character,                dimension(:), pointer :: fstringp
    character(len=TYPELEN+1), dimension(ntypes)     :: fatom_types

    integer :: i, slen

    do i = 1, ntypes
       call c_f_pointer(atom_types(i), fstringp, [TYPELEN+1])
       fatom_types(i) = transfer(fstringp(1:TYPELEN+1), fatom_types(i))
       slen = index(fatom_types(i), c_null_char) - 1
       fatom_types(i)(slen+1:TYPELEN+1) = ' '
    end do

    call aenet_init(fatom_types, stat)

  end subroutine aenet_init_C

  !--------------------------------------------------------------------!

  subroutine aenet_final(stat) bind(C)

    implicit none

    integer(kind=c_int), intent(out) :: stat

    integer :: itype
    integer :: ok

    stat = AENET_OK
    if (aenet_is_init) then
       call stp_final(aenet_ntypes, aenet_pot(1:aenet_ntypes)%stp)
       do itype = 1, aenet_ntypes
          call del_NNPot(aenet_pot(itype))
       end do
       deallocate(aenet_pot, aenet_atom_types, stat=ok)
       if (ok /= 0) then
          stat = AENET_ERR_MALLOC
          return
       end if
       aenet_ntypes = 0
       aenet_is_init = .false.
       aenet_is_loaded = .false.
    end if

  end subroutine aenet_final

  !--------------------------------------------------------------------!
  !                               output                               !
  !--------------------------------------------------------------------!

  subroutine aenet_print_info() bind(C)

    implicit none

    integer :: ipot

    if (.not. aenet_is_init) then
       write(*,*) "Nothing to print. AenetLib is not initialized."
    else
       do ipot = 1, aenet_ntypes
          if (aenet_pot(ipot)%init) then
             call pot_print_info(aenet_pot(ipot))
          end if
       end do
    end if

  end subroutine aenet_print_info

  !--------------------------------------------------------------------!
  !                        Load ANN potentials                         !
  !--------------------------------------------------------------------!

  subroutine aenet_load_potential(type_id, filename, stat)

    implicit none

    integer,          intent(in)  :: type_id
    character(len=*), intent(in)  :: filename
    integer,          intent(out) :: stat

    logical :: fexists

    stat = AENET_OK
    if (.not. aenet_is_init) then
       stat = AENET_ERR_INIT
       return
    end if

    if ((type_id <= 0) .or. (type_id > aenet_ntypes)) then
       stat = AENET_ERR_TYPE
       return
    end if

    inquire(file=trim(filename), exist=fexists)
    if (.not. fexists) then
       stat = AENET_ERR_IO
       return
    end if

    aenet_pot(type_id) = load_NNPot(aenet_atom_types, filename)
    aenet_nvalues_max = max(aenet_nvalues_max, aenet_pot(type_id)%net%nvalues)
    aenet_nweights_max = max(aenet_nweights_max, aenet_pot(type_id)%net%Wsize)

    ! when all potentials are loaded, determine array sizes
    if (aenet_all_loaded()) then
       call pot_get_range(aenet_ntypes, aenet_pot, aenet_Rc_min, aenet_Rc_max)
       if (stat /= 0) return
       aenet_nnb_max = lcl_nmax_nbdist(aenet_Rc_min, aenet_Rc_max)
       call stp_init(aenet_ntypes, aenet_pot(1:aenet_ntypes)%stp, aenet_nnb_max)
       aenet_nsf_max = stp_nsf_max()
       aenet_is_loaded = .true.
    end if

  end subroutine aenet_load_potential

  subroutine aenet_load_potential_C(type_id, filename, stat) &
       bind(C, name='aenet_load_potential')

    implicit none

    integer(kind=c_int), value,           intent(in)  :: type_id
    character(kind=c_char), dimension(*), intent(in)  :: filename
    integer(kind=c_int),                  intent(out) :: stat

    integer :: slen
    character(len=1024) :: ffilename

    slen = io_cstring_len(filename)
    if (slen > 1024) then
       stat = aenet_ERR_MALLOC
       return
    end if
    ffilename = io_cstring2f(filename, slen)

    call aenet_load_potential(type_id, trim(ffilename), stat)

  end subroutine aenet_load_potential_C

  function aenet_all_loaded() result(all_loaded) bind(C)

    implicit none

    logical(kind=c_bool) :: all_loaded
    integer :: ipot

    all_loaded = .true.
    do ipot = 1, aenet_ntypes
       if (.not. aenet_pot(ipot)%init) then
          all_loaded = .false.
          return
       end if
    end do

  end function aenet_all_loaded

  !--------------------------------------------------------------------!
  !                            information                             !
  !--------------------------------------------------------------------!

  function aenet_free_atom_energy(type_id) result(E_atom) bind(C)

    implicit none

    integer(kind=c_int), intent(in) :: type_id
    real(kind=c_double)             :: E_atom

    E_atom = aenet_pot(type_id)%E_atom

  end function aenet_free_atom_energy

  !--------------------------------------------------------------------!
  !                             Evaluation                             !
  !                                                                    !
  ! Attention: all routines require synchronized atom type IDs, i.e.,  !
  !            the IDs passed to the evaluation routines must be       !
  !            compatible with the ANN potential type IDs.             !
  !                                                                    !
  ! Notes:     * Coordinates are Cartesian.                            !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine aenet_atomic_energy(coo_i, type_i, n_j, coo_j, type_j, &
                                 E_i, stat) bind(C)

    implicit none

    real(kind=c_double), dimension(3),     intent(in)  :: coo_i
    integer(kind=c_int), value,            intent(in)  :: type_i
    integer(kind=c_int), value,            intent(in)  :: n_j
    real(kind=c_double), dimension(3,n_j), intent(in)  :: coo_j
    integer(kind=c_int), dimension(n_j),   intent(in)  :: type_j
    real(kind=c_double),                   intent(out) :: E_i
    integer(kind=c_int),                   intent(out) :: stat

    double precision, dimension(1)              :: E_i_arr
    double precision, dimension(:), allocatable :: sfval
    double precision, dimension(:), allocatable :: values
    double precision, dimension(:), allocatable :: derivs
    integer                                     :: nsf, nvalues
    integer                                     :: ok

    stat = AENET_OK
    if (.not. (aenet_is_init .and. aenet_is_loaded)) then
       stat = aenet_ERR_INIT
       return
    end if

    nvalues = aenet_pot(type_i)%net%nvalues
    allocate(sfval(aenet_nsf_max), &
             values(nvalues), &
             derivs(nvalues), stat=ok)
    if (ok /= 0) then
       stat = AENET_ERR_MALLOC
       return
    end if

    nsf = aenet_pot(type_i)%stp%nsf
    call stp_eval(type_i, coo_i, n_j, coo_j, type_j, &
                  aenet_pot(type_i)%stp, sfval=sfval, scaled=.true.)
    call ff_eval(aenet_pot(type_i)%net, nsf, sfval, 1, values, derivs, E_i_arr)
    E_i = aenet_pot(type_i)%E_scale*E_i_arr(1) + aenet_pot(type_i)%E_shift

    E_i = E_i + aenet_pot(type_i)%E_atom

    deallocate(sfval, values, derivs, stat=ok)
    if (ok /= 0) then
       stat = AENET_ERR_MALLOC
       return
    end if

  end subroutine aenet_atomic_energy

  subroutine aenet_atomic_energy_and_forces( &
       coo_i, type_i, index_i, n_j, coo_j, type_j, index_j, natoms, &
       E_i, F, stat) bind(C)

    implicit none

    real(kind=c_double), dimension(3),        intent(in)    :: coo_i
    integer(kind=c_int), value,               intent(in)    :: type_i
    integer(kind=c_int), value,               intent(in)    :: index_i
    integer(kind=c_int), value,               intent(in)    :: n_j
    real(kind=c_double), dimension(3,n_j),    intent(in)    :: coo_j
    integer(kind=c_int), dimension(n_j),      intent(in)    :: type_j
    integer(kind=c_int), dimension(n_j),      intent(in)    :: index_j
    integer(kind=c_int), value,               intent(in)    :: natoms
    real(kind=c_double),                      intent(out)   :: E_i
    real(kind=c_double), dimension(3,natoms), intent(inout) :: F
    integer(kind=c_int),                      intent(out)   :: stat

    double precision, dimension(1)                  :: E_i_arr
    double precision, dimension(:),     allocatable :: sfval
    double precision, dimension(:,:),   allocatable :: sfderiv_i
    double precision, dimension(:,:,:), allocatable :: sfderiv_j
    double precision, dimension(:),     allocatable :: values
    double precision, dimension(:),     allocatable :: derivs
    double precision, dimension(:),     allocatable :: jacobian
    integer                                         :: nvalues, nweights
    integer                                         :: nsf, j
    double precision, dimension(:),     allocatable :: dE_dG
    integer                                         :: ok

    stat = AENET_OK
    if (.not. (aenet_is_init .and. aenet_is_loaded)) then
       stat = aenet_ERR_INIT
       return
    end if

    nvalues = aenet_pot(type_i)%net%nvalues
    nweights = aenet_pot(type_i)%net%Wsize
    allocate(dE_dG(aenet_nsf_max),        &
             sfval(aenet_nsf_max),        &
             sfderiv_i(3, aenet_nsf_max), &
             sfderiv_j(3, aenet_nsf_max, aenet_nnb_max), &
             values(nvalues), &
             derivs(nvalues), &
             jacobian(nweights), stat=ok)
    if (ok /= 0) then
       stat = AENET_ERR_MALLOC
       return
    end if

    nsf = aenet_pot(type_i)%stp%nsf
    call stp_eval(type_i, coo_i, n_j, coo_j, type_j, &
                  aenet_pot(type_i)%stp, sfval=sfval, &
                  sfderiv_i=sfderiv_i, sfderiv_j=sfderiv_j, scaled=.true.)

    call ff_eval(aenet_pot(type_i)%net, nsf, sfval, 1, values, derivs, E_i_arr)
    call ff_deriv(aenet_pot(type_i)%net, nsf, 1, derivs, jacobian, dE_dG(1:nsf))

    E_i = aenet_pot(type_i)%E_scale*E_i_arr(1) + aenet_pot(type_i)%E_shift
    E_i = E_i + aenet_pot(type_i)%E_atom

    F(1:3, index_i) = F(1:3, index_i) - aenet_pot(type_i)%E_scale &
                    * matmul(sfderiv_i(1:3,1:nsf), dE_dG(1:nsf))

    do j = 1, n_j
       F(1:3, index_j(j)) = F(1:3, index_j(j)) - aenet_pot(type_i)%E_scale &
                          * matmul(sfderiv_j(1:3,1:nsf,j), dE_dG(1:nsf))
    end do

    deallocate(dE_dG, sfval, sfderiv_i, sfderiv_j, &
               values, derivs, jacobian, stat=ok)
    if (ok /= 0) then
       stat = AENET_ERR_MALLOC
       return
    end if

  end subroutine aenet_atomic_energy_and_forces

  !--------------------------------------------------------------------!
  !                       convert atom type IDs                        !
  !--------------------------------------------------------------------!

  subroutine aenet_convert_atom_types(&
       atom_types_in, type_id_in, type_id_out, stat)

    implicit none

    character(len=*), dimension(:), intent(in)  :: atom_types_in
    integer,          dimension(:), intent(in)  :: type_id_in
    integer,          dimension(:), intent(out) :: type_id_out
    integer,                        intent(out) :: stat

    integer :: iat,it,  nTypes_in, natoms_in

    stat = AENET_OK
    if (.not. aenet_is_init) then
       stat = aenet_ERR_INIT
       return
    end if

    ntypes_in = size(atom_types_in)
    natoms_in = size(type_id_in)

    do iat = 1, natoms_in
       it = geo_type_conv(type_id_in(iat), ntypes_in, atom_types_in, &
                          aenet_ntypes, aenet_atom_types)
       if (it == 0) then
          stat = aenet_ERR_TYPE
          return
       end if
       type_id_out(iat) = it
    end do

  end subroutine aenet_convert_atom_types

  subroutine aenet_convert_atom_types_C(ntypes_in, atom_types_in, &
                         natoms_in, type_id_in, type_id_out, stat &
                         ) bind(C, name="aenet_convert_atom_types")

    implicit none

    integer(kind=c_int), value,                        intent(in)  :: ntypes_in
    type(c_ptr),         dimension(ntypes_in), target, intent(in)  :: atom_types_in
    integer(kind=c_int), value,                        intent(in)  :: natoms_in
    integer(kind=c_int), dimension(natoms_in),         intent(in)  :: type_id_in
    integer(kind=c_int), dimension(natoms_in),         intent(out) :: type_id_out
    integer(kind=c_int),                               intent(out) :: stat

    character,                dimension(:), pointer :: fstringp
    character(len=TYPELEN+1), dimension(ntypes_in)  :: fatom_types

    integer :: i, slen

    do i = 1, ntypes_in
       call c_f_pointer(atom_types_in(i), fstringp, [TYPELEN+1])
       fatom_types(i) = transfer(fstringp(1:TYPELEN+1), fatom_types(i))
       slen = index(fatom_types(i), c_null_char) - 1
       fatom_types(i)(slen+1:TYPELEN+1) = ' '
    end do

    call aenet_convert_atom_types(&
         fatom_types, type_id_in, type_id_out, stat)

  end subroutine aenet_convert_atom_types_C


  !====================================================================!
  !                                                                    !
  !                           Neighbor List                            !
  !                                                                    !
  ! Interface to aenet's linked-cell neighbor list to determin the     !
  ! atoms within the local atomic environment.                         !
  !                                                                    !
  ! Most molecular dynamics codes provide an implementation of         !
  ! neighbor lists of some kind that is adapted to simulations where   !
  ! the atomic coordinates change little between steps, e.g., Verlet   !
  ! lists.  Those approaches are superior to linked-cell list provided !
  ! here and should be used whenever possible.                         !
  !                                                                    !
  !====================================================================!

  subroutine aenet_nbl_init(latt_vec, natoms, atom_types, coords, &
                            cartesian, pbc) bind(C)

    implicit none

    real(kind=c_double), dimension(3,3),       intent(in)    :: latt_vec
    integer(kind=c_int), value,                intent(in)    :: natoms
    integer(kind=c_int), dimension(natoms),    intent(in)    :: atom_types
    real(kind=c_double), dimension(3, natoms), intent(inout) :: coords
    logical(kind=c_bool), value,               intent(in)    :: cartesian
    logical(kind=c_bool), value,               intent(in)    :: pbc

    real(kind=c_double), dimension(3,3)       :: rec_latt_vec
    logical                                   :: pbc_f

    allocate(aenet_coo_latt(3, natoms))

    if (cartesian) then
       rec_latt_vec = geo_recip_lattice(latt_vec, cryst=.true.)
       aenet_coo_latt = matmul(rec_latt_vec, coords)
    else
       aenet_coo_latt = coords
    end if

    ! C to Fortran conversion
    if (pbc) then
       pbc_f = .true.
    else
       pbc_f = .false.
    end if

    call lcl_init(aenet_Rc_min, aenet_Rc_max, latt_vec, natoms, &
                  atom_types, aenet_coo_latt, pbc_f)

    coords = matmul(latt_vec, aenet_coo_latt)

  end subroutine aenet_nbl_init

  !--------------------------------------------------------------------!

  subroutine aenet_nbl_final() bind(C)

    implicit none

    call lcl_final()
    deallocate(aenet_coo_latt)

  end subroutine aenet_nbl_final

  !--------------------------------------------------------------------!

  subroutine aenet_nbl_neighbors(iatom, nnb, nbcoo, nbdist, nblist, &
                                 nbtype) bind(C)

    implicit none

    integer(kind=c_int), value,            intent(in)    :: iatom
    integer(kind=c_int),                   intent(inout) :: nnb
    real(kind=c_double), dimension(3,nnb), intent(out)   :: nbcoo
    real(kind=c_double), dimension(nnb),   intent(out)   :: nbdist
    integer(kind=c_int), dimension(nnb),   intent(out)   :: nblist
    integer(kind=c_int), dimension(nnb),   intent(out)   :: nbtype

    integer :: nnb_here

    call lcl_nbdist_cart(iatom, nnb, nbcoo, nbdist, aenet_Rc_max, &
                         nblist=nblist, nbtype=nbtype)

  end subroutine aenet_nbl_neighbors

  !====================================================================!
  !                                                                    !
  !                    structural fingerprint basis                    !
  !                                                                    !
  ! The procedures here are for use without ANN potentials.  The ANN   !
  ! potential code handles structural fingerprint setups separately in !
  ! sfsetup.f90.  At some point, this might become an independent      !
  ! library.                                                           !
  !                                                                    !
  !====================================================================!

  subroutine aenet_sfb_init(atom_types, radial_order, angular_order, &
                            radial_Rc, angular_Rc, stat)

    implicit none

    character(len=*), dimension(:), intent(in)  :: atom_types
    integer,                        intent(in)  :: radial_order
    integer,                        intent(in)  :: angular_order
    double precision,               intent(in)  :: radial_Rc
    double precision,               intent(in)  :: angular_Rc
    integer,                        intent(out) :: stat

    integer                                     :: num_types

    stat = AENET_OK
    if (aenet_sfb%initialized) then
       stat = AENET_ERR_INIT
       return
    end if

    num_types = size(atom_types)

    aenet_sfb = new_SFBasis(num_types, atom_types, radial_order, &
                            angular_order, radial_Rc, angular_Rc)

  end subroutine aenet_sfb_init

  !--------------------------------------------------------------------!

  subroutine aenet_sfb_init_C(ntypes, atom_types, radial_order, &
                              angular_order, radial_Rc, angular_Rc, &
                              stat) bind(C, name='aenet_sfb_init')

    implicit none

    integer(kind=c_int), value,             intent(in)  :: ntypes
    type(c_ptr), dimension(ntypes), target, intent(in)  :: atom_types
    integer(kind=c_int), value,             intent(in)  :: radial_order
    integer(kind=c_int), value,             intent(in)  :: angular_order
    real(kind=c_double), value,             intent(in)  :: radial_Rc
    real(kind=c_double), value,             intent(in)  :: angular_Rc
    integer(kind=c_int),                    intent(out) :: stat

    character,                dimension(:), pointer :: fstringp
    character(len=TYPELEN+1), dimension(ntypes)     :: fatom_types

    integer, parameter :: typelen_local = 2
    integer :: i, slen

    do i = 1, ntypes
       call c_f_pointer(atom_types(i), fstringp, [typelen_local+1])
       fatom_types(i) = transfer(fstringp(1:typelen_local+1), fatom_types(i))
       slen = index(fatom_types(i), c_null_char) - 1
       fatom_types(i)(slen+1:typelen_local+1) = ' '
    end do

    call aenet_sfb_init(fatom_types, radial_order, angular_order, &
                        radial_Rc, angular_Rc, stat)

  end subroutine aenet_sfb_init_C

  !--------------------------------------------------------------------!

  subroutine aenet_sfb_final(stat) bind(C)

    implicit none

    integer(kind=c_int), intent(out) :: stat

    stat = AENET_OK
    if (.not. aenet_sfb%initialized) then
       stat = AENET_ERR_INIT
       return
    end if

    call del_SFBasis(aenet_sfb)

  end subroutine aenet_sfb_final

  !--------------------------------------------------------------------!

  function aenet_sfb_nvalues() result(nvalues) bind(C)

    implicit none

    integer(kind=c_int) :: nvalues

    if (.not. aenet_sfb%initialized) then
       nvalues = 0
    else
       nvalues = aenet_sfb%N
    end if

  end function aenet_sfb_nvalues

  !--------------------------------------------------------------------!

  subroutine aenet_sfb_eval(itype0, coo0, nat, itype1, coo1, nv, &
                            values, stat) bind(C)

    implicit none

    integer(kind=c_int), value,            intent(in)    :: itype0
    real(kind=c_double), dimension(3),     intent(in)    :: coo0
    integer(kind=c_int), value,            intent(in)    :: nat
    integer(kind=c_int), dimension(nat),   intent(in)    :: itype1
    real(kind=c_double), dimension(3,nat), intent(in)    :: coo1
    integer(kind=c_int), value,            intent(in)    :: nv
    real(kind=c_double), dimension(nv),    intent(inout) :: values
    integer(kind=c_int),                   intent(out)   :: stat

    stat = AENET_OK
    if (.not. aenet_sfb%initialized) then
       stat = AENET_ERR_INIT
       return
    end if

    call sfb_eval(aenet_sfb, itype0, coo0, nat, itype1, coo1, nv, values)

  end subroutine aenet_sfb_eval

  !--------------------------------------------------------------------!

  subroutine aenet_sfb_reconstruct_radial(nv, values, nx, x, y, stat) bind(C)

    implicit none

    integer(kind=c_int), value,            intent(in)    :: nv
    real(kind=c_double), dimension(nv),    intent(in)    :: values
    integer(kind=c_int), value,            intent(in)    :: nx
    real(kind=c_double), dimension(nx),    intent(out)   :: x
    real(kind=c_double), dimension(nx),    intent(out)   :: y
    integer(kind=c_int),                   intent(out)   :: stat

    double precision, dimension(:), allocatable :: coeff

    stat = AENET_OK
    if (.not. aenet_sfb%initialized) then
       stat = AENET_ERR_INIT
       return
    end if

    allocate(coeff(aenet_sfb%r_N))
    coeff = values(aenet_sfb%r_i1:aenet_sfb%r_f1)
    call sfb_reconstruct_radial(aenet_sfb, coeff, nx, x, y)
    deallocate(coeff)

  end subroutine aenet_sfb_reconstruct_radial

  !--------------------------------------------------------------------!
  !                         internal routines                          !
  !--------------------------------------------------------------------!

end module aenet
