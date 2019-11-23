!-----------------------------------------------------------------------
!    parallel.f90 - F90 layer for the MPI message passing interface
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
!
! This module mainly provides generic interfaces to the FORTRAN 77 MPI
! bindings.  The MPI send, receive and broadcast operations are
! simplified for the case of a rank 0 master process and the single
! communicator MPI_COMM_WORLD.
!
! The prodedures for sending and receiving derived types defined in
! this module only send primitive data types.  This is probably not
! very efficient, but should be very portable.
!
!-----------------------------------------------------------------------
! 2011-11-10 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------


module parallel

  use aeio,        only: aeio_header

  use feedforward, only: Network,     &
                         new_Network

  use input,       only: InputData

  use io,          only: io_adjustl

  use potential,   only: NNPot

  use sfsetup,     only: Setup,       &
                         new_Setup

  use trainset,    only: TrnSet,      &
                         new_TrnSet_info

  implicit none
  save

#ifdef PARALLEL
  include 'mpif.h'
#endif

  public::  pp_init,              &
            pp_final,             &
            pp_print_info,        &
            pp_barrier,           &
            pp_bcast,             &
            pp_send,              &
            pp_recv,              &
            pp_sum,               &
            pp_sum2d,             &
            pp_send_Network,      &
            pp_recv_Network,      &
            pp_bcast_Network,     &
            pp_send_Setup,        &
            pp_recv_Setup,        &
            pp_bcast_Setup,       &
            pp_send_InputData,    &
            pp_recv_InputData,    &
            pp_bcast_InputData,   &
            pp_send_TrnSet_info,  &
            pp_recv_TrnSet_info,  &
            pp_bcast_TrnSet_info, &
            pp_send_NNPot,        &
            pp_recv_NNPot,        &
            pp_bcast_NNPot,       &
            pp_bcast_coo,         &
            pp_bcast_latt,        &
            pp_sum_weights

  !------------------------------ public ------------------------------!

  logical, public  :: ppMaster = .true.
  integer, public  :: ppRank   = 0
  integer, public  :: ppSize   = 1

  !----------------------------- private ------------------------------!

  integer, private :: ierr
  logical, private :: isInit     = .false.
  logical, private :: isParallel = .false.

#ifdef PARALLEL
  integer, dimension(MPI_STATUS_SIZE) :: status
#endif

  !--------------------------------------------------------------------!
  !                    generic BROADCAST interface                     !
  !                                                                    !
  ! pp_bcast(val[, root])                                              !
  ! pp_bcast(val(1:n), n[, root])                                      !
  !                                                                    !
  ! Broadcast either single values or arrays.  The default root        !
  ! process is 0.  The following implementations are available:        !
  !                                                                    !
  ! pp_bcast_i1(val)       -- single integer value                     !
  ! pp_bcast_in(val,n)     -- integer array of length n                !
  ! pp_bcast_d1(val)       -- single double precision value            !
  ! pp_bcast_dn(val,n)     -- double precision array of length n       !
  ! pp_bcast_l1(val)       -- single logical value                     !
  ! pp_bcast_ln(val)       -- array of logical values                  !
  ! pp_bcast_c1(val)       -- single character string                  !
  ! pp_bcast_cn(val)       -- array of character strings               !
  !--------------------------------------------------------------------!

  interface pp_bcast
     module procedure pp_bcast_i1, pp_bcast_in, pp_bcast_d1,      &
                      pp_bcast_dn, pp_bcast_l1, pp_bcast_ln,      &
                      pp_bcast_c1, pp_bcast_cn
  end interface

  !--------------------------------------------------------------------!
  !                       generic SEND interface                       !
  !                                                                    !
  ! pp_send(val, dest)                                                 !
  ! pp_send(val(1:n), n, dest)                                         !
  !                                                                    !
  ! Send a single value or an array do process `dest'.                 !
  !                                                                    !
  ! pp_send_i1(val, dest)    -- single integer value                   !
  ! pp_send_in(val,n, dest)  -- integer array of length n              !
  ! pp_send_d1(val, dest)    -- single double precision value          !
  ! pp_send_dn(val,n, dest)  -- double precision array of length n     !
  ! pp_send_l1(val, dest)    -- single logical value                   !
  ! pp_send_c1(val, dest)    -- single character string                !
  ! pp_send_cn(val, dest)    -- array of character strings             !
  !--------------------------------------------------------------------!

  interface pp_send
     module procedure pp_send_i1, pp_send_in, pp_send_d1, &
                      pp_send_dn, pp_send_l1, pp_send_c1, &
                      pp_send_cn
  end interface

  !--------------------------------------------------------------------!
  !                     generic RECEIVE interface                      !
  !                                                                    !
  ! pp_recv(val[, src])                                                !
  ! pp_recv(val(1:n), n[, src])                                        !
  !                                                                    !
  ! Receive a single value or array from source `src'.  The default    !
  ! source is MPI_ANY_SOURCE.                                          !
  !                                                                    !
  ! pp_recv_i1(val)    -- single integer value                         !
  ! pp_recv_in(val,n)  -- integer array of length n                    !
  ! pp_recv_d1(val)    -- single double precision value                !
  ! pp_recv_dn(val,n)  -- double precision array of length n           !
  ! pp_recv_l1(val)    -- single logical value                         !
  ! pp_recv_c1(val)    -- single character string                      !
  ! pp_recv_cn(val)    -- array of character strings                   !
  !--------------------------------------------------------------------!

  interface pp_recv
     module procedure pp_recv_i1, pp_recv_in, pp_recv_d1, &
                      pp_recv_dn, pp_recv_l1, pp_recv_c1, &
                      pp_recv_cn
  end interface

  !--------------------------------------------------------------------!
  !                       generic SUM interface                        !
  !                                                                    !
  ! pp_sum(val)                                                        !
  ! pp_sum(val(1:n), n)                                                !
  !                                                                    !
  ! Sum variables or arrays over all processes and return the result   !
  ! to every process (ALLREDUCE).                                      !
  !                                                                    !
  ! pp_sum_i1(val)    -- single integer value                          !
  ! pp_sum_in(val,n)  -- integer array of length n                     !
  ! pp_sum_d1(val)    -- single double precision value                 !
  ! pp_sum_dn(val,n)  -- double precision array of length n            !
  !--------------------------------------------------------------------!

  interface pp_sum
     module procedure pp_sum_i1, pp_sum_in, pp_sum_d1, &
                      pp_sum_dn
  end interface


contains !=============================================================!


  subroutine pp_init()

    implicit none

    if (isInit) return

#ifdef PARALLEL
    call MPI_Init(ierr)

    call MPI_Comm_size(MPI_COMM_WORLD, ppSize, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, ppRank, ierr)

    isParallel = .true.
#endif

    if (ppRank == 0) then
       ppMaster = .true.
    else
       ppMaster = .false.
    end if

    isInit = .true.

  end subroutine pp_init

  !--------------------------------------------------------------------!

  subroutine pp_final()

    implicit none

    if (.not. isInit) return

#ifdef PARALLEL
    call MPI_Finalize(ierr)
#endif

    isInit     = .false.
    isParallel = .false.

  end subroutine pp_final

  !--------------------------------------------------------------------!

  subroutine pp_print_info()

    implicit none

    if (.not. isInit)     return
    if (.not. ppMaster)   return
    if (.not. isParallel) return

    call aeio_header("Parallel run")
    write(*,*)

    write(*,*) 'Number of processes : ', trim(io_adjustl(ppSize))
    write(*,*)

  end subroutine pp_print_info

  !--------------------------------------------------------------------!
  !                       send/receive networks                        !
  !                                                                    !
  ! Actually only the information that's necessary to construct an     !
  ! identical network is communicated.                                 !
  !--------------------------------------------------------------------!

  subroutine pp_send_Network(net, dest)

    implicit none

    type(Network), intent(in) :: net
    integer,       intent(in) :: dest

    if (.not. isParallel) return

    ! send nlayer, nnodes, and f_a, so that the other process
    ! can allocate an equivalent network:

    call pp_send(net%nlayers, dest)
    call pp_send(net%nnodes, net%nlayers, dest)
    call pp_send(net%f_a, net%nlayers-1, dest)

    ! send network weights:

    call pp_send(net%W, net%Wsize, dest)

  end subroutine pp_send_Network

  !--------------------------------------------------------------------!

  function pp_recv_Network() result(net)

    implicit none

    type(Network)                       :: net

    integer                             :: nlayers
    integer, dimension(:), allocatable  :: nnodes
    integer, dimension(:), allocatable  :: f_a

    if (.not. isParallel) return

    ! receive nlayers, nnodes, f_a:

    call pp_recv(nlayers)
    allocate(nnodes(nlayers), f_a(2:nlayers))
    call pp_recv(nnodes(1:nlayers), nlayers)
    call pp_recv(f_a(2:nlayers), nlayers-1)

    ! allocate a network defined by nlayers, nnodes, f_a:

    net = new_Network(nnodes)
    net%f_a(:) = f_a(:)
    deallocate(nnodes, f_a)

    ! receive and store the network weights:

    call pp_recv(net%W, net%Wsize)

  end function pp_recv_Network

  !--------------------------------------------------------------------!

  subroutine pp_bcast_Network(net)

    implicit none

    type(Network), intent(inout) :: net

    integer :: dest

    if (.not. isParallel) return

    if (ppMaster) then
       do dest = 0, ppSize - 1
          if (dest == ppRank) cycle
          call pp_send_Network(net, dest)
       end do
    else
       net = pp_recv_Network()
    end if

  end subroutine pp_bcast_Network

  !--------------------------------------------------------------------!
  !                                                                    !
  !                 send/receive basis function set-up                 !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine pp_send_Setup(stp, dest)

    implicit none

    type(Setup), intent(in) :: stp
    integer,     intent(in) :: dest

    integer :: isf, nsfp, nenv

    if (.not. isParallel) return

    call pp_send(stp%nsf, dest)
    call pp_send(stp%nenv, dest)
    call pp_send(stp%ntypes_global, dest)

    call pp_send(stp%neval, dest)
    call pp_send(stp%atomtype, dest)
    call pp_send(stp%Rc_min, dest)
    call pp_send(stp%Rc_max, dest)
    call pp_send(stp%sftype, dest)

    call pp_send(stp%envtypes, stp%nenv, dest)
    call pp_send(stp%gtype, stp%nenv, dest)
    call pp_send(stp%ltype, stp%ntypes_global, dest)

    call pp_send(stp%sf, stp%nsf, dest)
    nsfp = stp%nsfparam
    call pp_send(nsfp, dest)
    nenv = size(stp%sfenv(:,isf))
    do isf = 1, stp%nsf
       call pp_send(stp%sfparam(1:nsfp,isf), nsfp, dest)
       call pp_send(stp%sfenv(1:nenv,isf), nenv, dest)
    end do
    call pp_send(stp%sfval_min(1:stp%nsf), stp%nsf, dest)
    call pp_send(stp%sfval_max(1:stp%nsf), stp%nsf, dest)
    call pp_send(stp%sfval_avg(1:stp%nsf), stp%nsf, dest)
    call pp_send(stp%sfval_cov(1:stp%nsf), stp%nsf, dest)

  end subroutine pp_send_Setup

  !--------------------------------------------------------------------!

  function pp_recv_Setup() result(stp)

    implicit none

    type(Setup) :: stp
    integer     :: isf, nsf, nenv, ntypes_global, nsfp

    if (.not. isParallel) return

    call pp_recv(nsf)
    call pp_recv(nenv)
    call pp_recv(ntypes_global)
    stp = new_Setup(nsf, nenv, ntypes_global)

    call pp_recv(stp%neval)
    call pp_recv(stp%atomtype)
    call pp_recv(stp%Rc_min)
    call pp_recv(stp%Rc_max)
    call pp_recv(stp%sftype)

    call pp_recv(stp%envtypes, stp%nenv)
    call pp_recv(stp%gtype, stp%nenv)
    call pp_recv(stp%ltype, stp%ntypes_global)

    call pp_recv(stp%sf, nsf)
    call pp_recv(nsfp)
    nenv = size(stp%sfenv(:,isf))
    do isf = 1, stp%nsf
       call pp_recv(stp%sfparam(1:nsfp,isf), nsfp)
       call pp_recv(stp%sfenv(1:nenv,isf), nenv)
    end do
    call pp_recv(stp%sfval_min(1:stp%nsf), nsf)
    call pp_recv(stp%sfval_max(1:stp%nsf), nsf)
    call pp_recv(stp%sfval_avg(1:stp%nsf), nsf)
    call pp_recv(stp%sfval_cov(1:stp%nsf), nsf)

  end function pp_recv_Setup

  !--------------------------------------------------------------------!

  subroutine pp_bcast_Setup(stp)

    implicit none

    type(Setup), intent(inout) :: stp

    integer :: dest

    if (.not. isParallel) return

    if (ppMaster) then
       do dest = 0, ppSize - 1
          if (dest == ppRank) cycle
          call pp_send_Setup(stp, dest)
       end do
    else
       stp = pp_recv_Setup()
    end if

  end subroutine pp_bcast_Setup

  !--------------------------------------------------------------------!
  !                                                                    !
  !                  send/receive training set info                    !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine pp_send_TrnSet_info(ts, dest)

    implicit none

    type(TrnSet), intent(in) :: ts
    integer,      intent(in) :: dest

    if (.not. isParallel) return

    call pp_send(ts%init, dest)
    if (.not. ts%init) return

    call pp_send(ts%nTypes, dest)

    call pp_send(ts%normalized, dest)
    call pp_send(ts%scale, dest)
    call pp_send(ts%shift, dest)
    call pp_send(ts%nAtomsTot, dest)
    call pp_send(ts%nStrucs, dest)
    call pp_send(ts%iStruc, dest)
    call pp_send(ts%E_min, dest)
    call pp_send(ts%E_max, dest)
    call pp_send(ts%E_av, dest)

    call pp_send(ts%typeName, ts%nTypes, dest)
    call pp_send(ts%E_atom, ts%nTypes, dest)

  end subroutine pp_send_TrnSet_info

  !--------------------------------------------------------------------!

  function pp_recv_TrnSet_info() result(ts)

    implicit none

    type(TrnSet) :: ts
    integer      :: nTypes

    if (.not. isParallel) return

    call pp_recv(ts%init)
    if (.not. ts%init) return

    call pp_recv(nTypes)
    ts = new_TrnSet_info(nTypes)
    ts%file = ''
    ts%unit = -1
    ts%mode = 'info'

    call pp_recv(ts%normalized)
    call pp_recv(ts%scale)
    call pp_recv(ts%shift)
    call pp_recv(ts%nAtomsTot)
    call pp_recv(ts%nStrucs)
    call pp_recv(ts%iStruc)
    call pp_recv(ts%E_min)
    call pp_recv(ts%E_max)
    call pp_recv(ts%E_av)

    call pp_recv(ts%typeName, ts%nTypes)
    call pp_recv(ts%E_atom, ts%nTypes)

  end function pp_recv_TrnSet_info

  !--------------------------------------------------------------------!

  subroutine pp_bcast_TrnSet_info(ts)

    implicit none

    type(TrnSet), intent(inout) :: ts

    integer :: dest

    if (.not. isParallel) return

    if (ppMaster) then
       do dest = 0, ppSize - 1
          if (dest == ppRank) cycle
          call pp_send_TrnSet_info(ts, dest)
       end do
    else
       ts = pp_recv_TrnSet_info()
    end if

  end subroutine pp_bcast_TrnSet_info

  !--------------------------------------------------------------------!
  !                                                                    !
  !                 send/receive NN potential objects                  !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine pp_send_NNPot(pot, dest)

    implicit none

    type(NNPot), intent(in) :: pot
    integer,     intent(in) :: dest

    if (.not. isParallel) return

    call pp_send(pot%init, dest)
    if (.not. pot%init) return

    call pp_send(pot%file, dest)
    call pp_send(pot%unit, dest)
    call pp_send(pot%typeName, dest)
    call pp_send(pot%E_atom, dest)
    call pp_send(pot%E_scale, dest)
    call pp_send(pot%E_shift, dest)

    call pp_send_Network(pot%net, dest)
    call pp_send_Setup(pot%stp, dest)
    call pp_send_TrnSet_info(pot%ts, dest)

  end subroutine pp_send_NNPot

  !--------------------------------------------------------------------!

  function pp_recv_NNPot() result(pot)

    implicit none

    type(NNPot) :: pot

    if (.not. isParallel) return

    call pp_recv(pot%init)
    if (.not. pot%init) return

    call pp_recv(pot%file)
    call pp_recv(pot%unit)
    call pp_recv(pot%typeName)
    call pp_recv(pot%E_atom)
    call pp_recv(pot%E_scale)
    call pp_recv(pot%E_shift)
    pot%net = pp_recv_Network()
    pot%stp = pp_recv_Setup()
    pot%ts  = pp_recv_TrnSet_info()

  end function pp_recv_NNPot

  !--------------------------------------------------------------------!

  subroutine pp_bcast_NNPot(pot)

    implicit none

    type(NNPot), intent(inout) :: pot

    integer :: dest

    if (.not. isParallel) return

    if (ppMaster) then
       do dest = 0, ppSize - 1
          if (dest == ppRank) cycle
          call pp_send_NNPot(pot, dest)
       end do
    else
       pot = pp_recv_NNPot()
    end if

  end subroutine pp_bcast_NNPot


  !--------------------------------------------------------------------!
  !                                                                    !
  !                   communicate input data objects                   !
  !                                                                    !
  !--------------------------------------------------------------------!

  subroutine pp_send_InputData(inp, dest)

    implicit none

    type(InputData), intent(in) :: inp
    integer,         intent(in) :: dest

    integer :: ig
    logical :: isalloc

    if (.not. isParallel) return

    call pp_send(inp%init, dest)
    if (.not. inp%init) return

    call pp_send(inp%mode, dest)
    call pp_send(inp%nTypes, dest)
    call pp_send(inp%do_debug, dest)

    isalloc = allocated(inp%typeName)
    call pp_send(isalloc, dest)
    if (isalloc) call pp_send(inp%typeName, inp%nTypes, dest)

    isalloc = allocated(inp%atomicEnergy)
    call pp_send(isalloc, dest)
    if (isalloc) call pp_send(inp%atomicEnergy, inp%nTypes, dest)

    isalloc = allocated(inp%netFile)
    call pp_send(isalloc, dest)
    if (isalloc) call pp_send(inp%netFile,inp%nTypes, dest)

    isalloc = allocated(inp%netArch)
    call pp_send(isalloc, dest)
    if (isalloc) call pp_send(inp%netArch,inp%nTypes, dest)

    call pp_send(inp%do_forces, dest)
    call pp_send(inp%do_timing, dest)
    call pp_send(inp%nStrucs, dest)

    isalloc = allocated(inp%strucFile)
    call pp_send(isalloc, dest)
    if (isalloc) call pp_send(inp%strucFile, inp%nStrucs, dest)

    call pp_send(inp%do_relax, dest)
    call pp_send(inp%relax_method, dest)
    call pp_send(inp%relax_steps, dest)
    call pp_send(inp%relax_F_conv, dest)
    call pp_send(inp%relax_E_conv, dest)

    call pp_send(inp%trn_file, dest)
    call pp_send(inp%trn_testset, dest)
    call pp_send(inp%trn_maxenergy, dest)
    call pp_send(inp%trn_steps, dest)
    call pp_send(inp%trn_sampling, dest)
    call pp_send(inp%trn_method, dest)
    call pp_send(inp%trn_nparams, dest)
    call pp_send(inp%do_save_energies, dest)
    isalloc = allocated(inp%trn_param)
    call pp_send(isalloc, dest)
    if (isalloc) call pp_send(inp%trn_param, inp%trn_nparams, dest)

    call pp_send(inp%T, dest)
    call pp_send(inp%T_final, dest)
    call pp_send(inp%nSteps, dest)
    call pp_send(inp%nSweeps, dest)
    call pp_send(inp%ensemble, dest)

    isalloc = allocated(inp%mu)
    call pp_send(isalloc, dest)
    if (isalloc) call pp_send(inp%mu, inp%nTypes, dest)

    isalloc = allocated(inp%activeType)
    call pp_send(isalloc, dest)
    if (isalloc) call pp_send(inp%activeType, inp%nTypes, dest)

    call pp_send(inp%mc_relax_final, dest)

    call pp_send(inp%mc_ngroups, dest)
    if (inp%mc_ngroups > 0) then
       call pp_send(inp%mc_group, inp%nTypes, dest)
       call pp_send(inp%mc_ntypes_group, inp%mc_ngroups, dest)
       do ig = 1, inp%mc_ngroups
          call pp_send(inp%mc_group_type(ig,:), inp%nTypes, dest)
       end do
    end if

    call pp_send(inp%mc_relax_final, dest)

  end subroutine pp_send_InputData

  !--------------------------------------------------------------------!

  function pp_recv_InputData() result(inp)

    implicit none

    type(InputData) :: inp

    integer :: ig
    logical :: isalloc

    if (.not. isParallel) return

    call pp_recv(inp%init)
    if (.not. inp%init) return

    call pp_recv(inp%mode)
    call pp_recv(inp%nTypes)
    call pp_recv(inp%do_debug)

    call pp_recv(isalloc)
    if (isalloc) then
       allocate(inp%typeName(inp%nTypes))
       call pp_recv(inp%typeName, inp%nTypes)
    end if

    call pp_recv(isalloc)
    if (isalloc) then
       allocate(inp%atomicEnergy(inp%nTypes))
       call pp_recv(inp%atomicEnergy, inp%nTypes)
    end if

    call pp_recv(isalloc)
    if (isalloc) then
       allocate(inp%netFile(inp%nTypes))
       call pp_recv(inp%netFile(1:inp%nTypes),inp%nTypes)
    end if

    call pp_recv(isalloc)
    if (isalloc) then
       allocate(inp%netArch(inp%nTypes))
       call pp_recv(inp%netArch(1:inp%nTypes),inp%nTypes)
    end if

    call pp_recv(inp%do_forces)
    call pp_recv(inp%do_timing)
    call pp_recv(inp%nStrucs)

    call pp_recv(isalloc)
    if (isalloc) then
       allocate(inp%strucFile(inp%nStrucs))
       call pp_recv(inp%strucFile, inp%nStrucs)
    end if

    call pp_recv(inp%do_relax)
    call pp_recv(inp%relax_method)
    call pp_recv(inp%relax_steps)
    call pp_recv(inp%relax_F_conv)
    call pp_recv(inp%relax_E_conv)

    call pp_recv(inp%trn_file)
    call pp_recv(inp%trn_testset)
    call pp_recv(inp%trn_maxenergy)
    call pp_recv(inp%trn_steps)
    call pp_recv(inp%trn_sampling)
    call pp_recv(inp%trn_method)
    call pp_recv(inp%trn_nparams)
    call pp_recv(inp%do_save_energies)
    call pp_recv(isalloc)
    if (isalloc) then
       allocate(inp%trn_param(inp%trn_nparams))
       call pp_recv(inp%trn_param, inp%trn_nparams)
    end if

    call pp_recv(inp%T)
    call pp_recv(inp%T_final)
    call pp_recv(inp%nSteps)
    call pp_recv(inp%nSweeps)
    call pp_recv(inp%ensemble)

    call pp_recv(isalloc)
    if (isalloc) then
       allocate(inp%mu(inp%nTypes))
       call pp_recv(inp%mu, inp%nTypes)
    end if

    call pp_recv(isalloc)
    if (isalloc) then
       allocate(inp%activeType(inp%nTypes))
       call pp_recv(inp%activeType, inp%nTypes)
    end if

    call pp_recv(inp%mc_relax_final)

    call pp_recv(inp%mc_ngroups)
    if (inp%mc_ngroups > 0) then
       allocate(inp%mc_group(inp%nTypes),                  &
                inp%mc_ntypes_group(inp%mc_ngroups),       &
                inp%mc_group_type(inp%mc_ngroups, inp%nTypes))
       call pp_recv(inp%mc_group, inp%nTypes)
       call pp_recv(inp%mc_ntypes_group, inp%mc_ngroups)
       do ig = 1, inp%mc_ngroups
          call pp_recv(inp%mc_group_type(ig,1:inp%nTypes), inp%nTypes)
       end do
    end if

    call pp_recv(inp%mc_relax_final)

  end function pp_recv_InputData

  !--------------------------------------------------------------------!

  subroutine pp_bcast_InputData(inp)

    implicit none

    type(InputData), intent(inout) :: inp

    integer :: dest

    if (.not. isParallel) return

    if (ppMaster) then
       do dest = 0, ppSize - 1
          if (dest == ppRank) cycle
          call pp_send_InputData(inp, dest)
       end do
    else
       inp = pp_recv_InputData()
    end if

  end subroutine pp_bcast_InputData

  !--------------------------------------------------------------------!
  !                  reduce / sum up network weights                   !
  !--------------------------------------------------------------------!

  subroutine pp_sum_weights(Dw, nw, ntypes)

    implicit none

    integer,                                intent(in) :: nw
    integer,                                intent(in) :: ntypes
    double precision, dimension(nw,ntypes), intent(inout) :: Dw

    integer :: itype

    if (.not. isParallel) return

    do itype = 1, ntypes
       call pp_sum(Dw(1:nw,itype), nw)
    end do

  end subroutine pp_sum_weights

  subroutine pp_sum2d(A, n, m)

    implicit none

    integer,                          intent(in)    :: n
    integer,                          intent(in)    :: m
    double precision, dimension(n,m), intent(inout) :: A

    integer :: i

    if (.not. isParallel) return

    do i = 1, m
       call pp_sum(A(1:n,i), n)
    end do

  end subroutine pp_sum2d

  !--------------------------------------------------------------------!
  !           broadcast 2D arrays of forces and coordinates            !
  !--------------------------------------------------------------------!

  subroutine pp_bcast_coo(coo, nat)

    implicit none

    integer,                            intent(in)    :: nat
    double precision, dimension(3,nat), intent(inout) :: coo

    integer :: iat

    if (.not. isParallel) return

    do iat = 1, nat
       call pp_bcast(coo(1:3,iat),3)
    end do

  end subroutine pp_bcast_coo

  !--------------------------------------------------------------------!
  !                     broadcast lattice vectors                      !
  !--------------------------------------------------------------------!

  subroutine pp_bcast_latt(avec)

    implicit none

    double precision, dimension(3,3), intent(inout) :: avec

    if (.not. isParallel) return

    call pp_bcast(avec(1:3,1),3)
    call pp_bcast(avec(1:3,2),3)
    call pp_bcast(avec(1:3,3),3)

  end subroutine pp_bcast_latt


  !====================================================================!
  !                                                                    !
  !                        basic MPI operations                        !
  !                                                                    !
  !====================================================================!


  !--------------------------------------------------------------------!
  !                            MPI barrier                             !
  !--------------------------------------------------------------------!

  subroutine pp_barrier()

    implicit none

    if (.not. isParallel) return

#ifdef PARALLEL
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
#endif

  end subroutine pp_barrier


  !====================================================================!
  !                                                                    !
  !         implementation of the generic BROADCAST interface          !
  !                                                                    !
  !====================================================================!


  !-------------------------- single integer --------------------------!

  subroutine pp_bcast_i1(val, root)

    implicit none

    integer,           intent(inout) :: val
    integer, optional, intent(in)    :: root

    if (.not. isParallel) return

#ifdef PARALLEL
    if (present(root)) then
       call MPI_Bcast(val, 1, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    else
       call MPI_Bcast(val, 1, MPI_INTEGER, 0,    MPI_COMM_WORLD, ierr)
    end if
#else
    write(*,*) val, root
#endif

  end subroutine pp_bcast_i1

  !-------------------------- integer array ---------------------------!

  subroutine pp_bcast_in(val, n, root)

    implicit none

    integer,               intent(in)    :: n
    integer, dimension(n), intent(inout) :: val
    integer, optional,     intent(in)    :: root

    if (.not. isParallel) return

#ifdef PARALLEL
    if (present(root)) then
       call MPI_Bcast(val, n, MPI_INTEGER, root, MPI_COMM_WORLD, ierr)
    else
       call MPI_Bcast(val, n, MPI_INTEGER, 0,    MPI_COMM_WORLD, ierr)
    end if
#else
    write(*,*) n, val, root
#endif

  end subroutine pp_bcast_in

  !--------------------- single double precision ----------------------!

  subroutine pp_bcast_d1(val, root)

    implicit none

    double precision,  intent(inout) :: val
    integer, optional, intent(in)    :: root

    if (.not. isParallel) return

#ifdef PARALLEL
    if (present(root)) then
       call MPI_Bcast(val, 1, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
    else
       call MPI_Bcast(val, 1, MPI_DOUBLE_PRECISION, 0,    MPI_COMM_WORLD, ierr)
    end if
#else
    write(*,*) val, root
#endif

  end subroutine pp_bcast_d1

  !---------------------- double precision array ----------------------!

  subroutine pp_bcast_dn(val, n, root)

    implicit none

    integer,                        intent(in)    :: n
    double precision, dimension(n), intent(inout) :: val
    integer, optional,              intent(in)    :: root

    if (.not. isParallel) return

#ifdef PARALLEL
    if (present(root)) then
       call MPI_Bcast(val, n, MPI_DOUBLE_PRECISION, root, MPI_COMM_WORLD, ierr)
    else
       call MPI_Bcast(val, n, MPI_DOUBLE_PRECISION, 0,    MPI_COMM_WORLD, ierr)
    end if
#else
    write(*,*) n, val, root
#endif

  end subroutine pp_bcast_dn

  !----------------------------- logical ------------------------------!

  subroutine pp_bcast_l1(val, root)

    implicit none

    logical,           intent(inout) :: val
    integer, optional, intent(in)    :: root

    if (.not. isParallel) return

#ifdef PARALLEL
    if (present(root)) then
       call MPI_Bcast(val, 1, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
    else
       call MPI_Bcast(val, 1, MPI_LOGICAL, 0,    MPI_COMM_WORLD, ierr)
    end if
#else
    write(*,*) val, root
#endif

  end subroutine pp_bcast_l1

  !------------------------ array of logicals -------------------------!

  subroutine pp_bcast_ln(val, n, root)

    implicit none

    integer,               intent(in)    :: n
    logical, dimension(n), intent(inout) :: val
    integer, optional,     intent(in)    :: root

    if (.not. isParallel) return

#ifdef PARALLEL
    if (present(root)) then
       call MPI_Bcast(val, n, MPI_LOGICAL, root, MPI_COMM_WORLD, ierr)
    else
       call MPI_Bcast(val, n, MPI_LOGICAL, 0,    MPI_COMM_WORLD, ierr)
    end if
#else
    write(*,*) n, val, root
#endif

  end subroutine pp_bcast_ln

  !---------------------------- character -----------------------------!

  subroutine pp_bcast_c1(val, root)

    implicit none

    character(len=*),  intent(inout) :: val
    integer, optional, intent(in)    :: root

    integer :: n

    if (.not. isParallel) return

    n = len(val)

#ifdef PARALLEL
    if (present(root)) then
       call MPI_Bcast(val, n, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
    else
       call MPI_Bcast(val, n, MPI_CHARACTER, 0,    MPI_COMM_WORLD, ierr)
    end if
#else
    n = root
    write(*,*) val, root
#endif

  end subroutine pp_bcast_c1

  !------------------------- array of strings -------------------------!

  subroutine pp_bcast_cn(val, n, root)

    implicit none

    integer,                        intent(in)    :: n
    character(len=*), dimension(n), intent(inout) :: val
    integer, optional,              intent(in)    :: root

    integer :: m

    if (.not. isParallel) return

    m = n*len(val)

#ifdef PARALLEL
    if (present(root)) then
       call MPI_Bcast(val, m, MPI_CHARACTER, root, MPI_COMM_WORLD, ierr)
    else
       call MPI_Bcast(val, m, MPI_CHARACTER, 0,    MPI_COMM_WORLD, ierr)
    end if
#else
    write(*,*) n, val, root
#endif

  end subroutine pp_bcast_cn



  !====================================================================!
  !                                                                    !
  !            implementation of the generic SEND interface            !
  !                                                                    !
  !====================================================================!


  !-------------------------- single integer --------------------------!

  subroutine pp_send_i1(val, dest)

    implicit none

    integer, intent(in) :: val
    integer, intent(in) :: dest

    integer,  parameter :: tag = 1

    if (.not. isParallel) return

#ifdef PARALLEL
    call MPI_Send(val, 1, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, ierr)
#else
    write(*,*) val, dest
#endif

  end subroutine pp_send_i1

  !-------------------------- integer array ---------------------------!

  subroutine pp_send_in(val, n, dest)

    implicit none

    integer,               intent(in) :: n
    integer, dimension(n), intent(in) :: val
    integer,               intent(in) :: dest

    integer,  parameter :: tag = 1

    if (.not. isParallel) return

#ifdef PARALLEL
    call MPI_Send(val, n, MPI_INTEGER, dest, tag, MPI_COMM_WORLD, ierr)
#else
    write(*,*) n, val, dest
#endif

  end subroutine pp_send_in

  !--------------------- single double precision ----------------------!

  subroutine pp_send_d1(val, dest)

    implicit none

    double precision, intent(in) :: val
    integer,          intent(in) :: dest

    integer,  parameter :: tag = 1

    if (.not. isParallel) return

#ifdef PARALLEL
    call MPI_Send(val, 1, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
#else
    write(*,*) val, dest
#endif

  end subroutine pp_send_d1

  !---------------------- double precision array ----------------------!

  subroutine pp_send_dn(val, n, dest)

    implicit none

    integer,                        intent(in) :: n
    double precision, dimension(n), intent(in) :: val
    integer,                        intent(in) :: dest

    integer,  parameter :: tag = 1

    if (.not. isParallel) return

#ifdef PARALLEL
    call MPI_Send(val, n, MPI_DOUBLE_PRECISION, dest, tag, MPI_COMM_WORLD, ierr)
#else
    write(*,*) n, val, dest
#endif

  end subroutine pp_send_dn

  !----------------------------- logical ------------------------------!

  subroutine pp_send_l1(val, dest)

    implicit none

    logical, intent(in) :: val
    integer, intent(in) :: dest

    integer,  parameter :: tag = 1

    if (.not. isParallel) return

#ifdef PARALLEL
    call MPI_Send(val, 1, MPI_LOGICAL, dest, tag, MPI_COMM_WORLD, ierr)
#else
    write(*,*) val, dest
#endif

  end subroutine pp_send_l1

  !---------------------------- character -----------------------------!

  subroutine pp_send_c1(val, dest)

    implicit none

    character(len=*), intent(in) :: val
    integer,          intent(in) :: dest

    integer,  parameter :: tag = 1

    integer :: n

    if (.not. isParallel) return

#ifdef PARALLEL
    n = len(val)
    call MPI_Send(val, n, MPI_CHARACTER, dest, tag, MPI_COMM_WORLD, ierr)
#else
    n = dest
    write(*,*) val
#endif

  end subroutine pp_send_c1

  !--------------------------- N characters ---------------------------!

  subroutine pp_send_cn(val, n, dest)

    implicit none

    integer,                        intent(in) :: n
    character(len=*), dimension(n), intent(in) :: val
    integer,                        intent(in) :: dest

    integer,  parameter :: tag = 1

    integer :: m

    if (.not. isParallel) return

#ifdef PARALLEL
    m = len(val(1))
    call MPI_Send(val, m*n, MPI_CHARACTER, dest, tag, MPI_COMM_WORLD, ierr)
#else
    m = n*dest
    write(*,*) val
#endif

  end subroutine pp_send_cn


  !====================================================================!
  !                                                                    !
  !          implementation of the generic RECEIVE interface           !
  !                                                                    !
  !====================================================================!


  !-------------------------- single integer --------------------------!

  subroutine pp_recv_i1(val, src)

    implicit none

    integer,                intent(out) :: val
    integer,      optional, intent(in)  :: src

    integer,                  parameter :: tag = 1

    if (.not. isParallel) return

#ifdef PARALLEL
    if (present(src)) then
       call MPI_Recv(val, 1, MPI_INTEGER, src, tag,       &
                     MPI_COMM_WORLD, status, ierr)
    else
       call MPI_Recv(val, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                     tag, MPI_COMM_WORLD, status, ierr)
    end if
#else
    val = src
#endif

  end subroutine pp_recv_i1

  !-------------------------- integer array ---------------------------!

  subroutine pp_recv_in(val, n, src)

    implicit none

    integer,               intent(in)   :: n
    integer, dimension(n), intent(out)  :: val
    integer, optional,     intent(in)   :: src

    integer,                  parameter :: tag = 1

    if (.not. isParallel) return

#ifdef PARALLEL
    if (present(src)) then
       call MPI_Recv(val, n, MPI_INTEGER, src, tag,       &
                     MPI_COMM_WORLD, status, ierr)
    else
       call MPI_Recv(val, n, MPI_INTEGER, MPI_ANY_SOURCE, &
                     tag, MPI_COMM_WORLD, status, ierr)
    end if
#else
    val = n*src
#endif

  end subroutine pp_recv_in

  !--------------------- single double precision ----------------------!

  subroutine pp_recv_d1(val, src)

    implicit none

    double precision,       intent(out) :: val
    integer,      optional, intent(in)  :: src

    integer,                  parameter :: tag = 1

    if (.not. isParallel) return

#ifdef PARALLEL
    if (present(src)) then
       call MPI_Recv(val, 1, MPI_DOUBLE_PRECISION, src, tag,       &
                     MPI_COMM_WORLD, status, ierr)
    else
       call MPI_Recv(val, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                     tag, MPI_COMM_WORLD, status, ierr)
    end if
#else
    val = 1.0*src
#endif

  end subroutine pp_recv_d1

  !---------------------- double precision array ----------------------!

  subroutine pp_recv_dn(val, n, src)

    implicit none

    integer,                        intent(in)   :: n
    double precision, dimension(n), intent(out)  :: val
    integer, optional,              intent(in)   :: src

    integer,                           parameter :: tag = 1

    if (.not. isParallel) return

#ifdef PARALLEL
    if (present(src)) then
       call MPI_Recv(val, n, MPI_DOUBLE_PRECISION, src, tag,       &
                     MPI_COMM_WORLD, status, ierr)
    else
       call MPI_Recv(val, n, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
                     tag, MPI_COMM_WORLD, status, ierr)
    end if
#else
    val = 1.0d0*n*src
#endif

  end subroutine pp_recv_dn

  !----------------------------- logical ------------------------------!

  subroutine pp_recv_l1(val, src)

    implicit none

    logical,                intent(out) :: val
    integer,      optional, intent(in)  :: src

    integer,                  parameter :: tag = 1

    if (.not. isParallel) return

#ifdef PARALLEL
    if (present(src)) then
       call MPI_Recv(val, 1, MPI_LOGICAL, src, tag,       &
                     MPI_COMM_WORLD, status, ierr)
    else
       call MPI_Recv(val, 1, MPI_LOGICAL, MPI_ANY_SOURCE, &
                     tag, MPI_COMM_WORLD, status, ierr)
    end if
#else
    val = (src == 1)
#endif

  end subroutine pp_recv_l1

  !---------------------------- character -----------------------------!

  subroutine pp_recv_c1(val, src)

    implicit none

    character(len=*),       intent(out) :: val
    integer,      optional, intent(in)  :: src

    integer,                  parameter :: tag = 1

    integer :: n

    if (.not. isParallel) return

#ifdef PARALLEL
    n = len(val)
    if (present(src)) then
       call MPI_Recv(val, n, MPI_CHARACTER, src, tag,       &
                     MPI_COMM_WORLD, status, ierr)
    else
       call MPI_Recv(val, n, MPI_CHARACTER, MPI_ANY_SOURCE, &
                     tag, MPI_COMM_WORLD, status, ierr)
    end if
#else
    n = src
    val = ' '
#endif

  end subroutine pp_recv_c1

  !----------------------- N character strings ------------------------!

  subroutine pp_recv_cn(val, n, src)

    implicit none

    integer,                        intent(in)  :: n
    character(len=*), dimension(n), intent(out) :: val
    integer, optional,              intent(in)  :: src

    integer,                  parameter :: tag = 1

    integer :: m

    if (.not. isParallel) return

#ifdef PARALLEL
    m = len(val(1))
    if (present(src)) then
       call MPI_Recv(val, n*m, MPI_CHARACTER, src, tag,       &
                     MPI_COMM_WORLD, status, ierr)
    else
       call MPI_Recv(val, n*m, MPI_CHARACTER, MPI_ANY_SOURCE, &
                     tag, MPI_COMM_WORLD, status, ierr)
    end if
#else
    m = src
    val = ' '
#endif

  end subroutine pp_recv_cn

  !====================================================================!
  !                                                                    !
  !           implementation of a generic MPI sum interface            !
  !                                                                    !
  !====================================================================!


  !-------------------------- single integer --------------------------!

  subroutine pp_sum_i1(val)

    implicit none

    integer, intent(inout) :: val
    integer                :: buff

    if (.not. isParallel) return

    buff = val
#ifdef PARALLEL
    call MPI_allreduce(buff, val, 1, MPI_INTEGER, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
#endif

  end subroutine pp_sum_i1

  !-------------------------- integer array ---------------------------!

  subroutine pp_sum_in(val,n)

    implicit none

    integer,               intent(in)    :: n
    integer, dimension(n), intent(inout) :: val
    integer, dimension(n)                :: buff

    if (.not. isParallel) return

    buff(1:n) = val(1:n)
#ifdef PARALLEL
    call MPI_allreduce(buff, val, n, MPI_INTEGER, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
#endif

  end subroutine pp_sum_in

  !------------------------- double precision -------------------------!

  subroutine pp_sum_d1(val)

    implicit none

    double precision, intent(inout) :: val
    double precision                :: buff

    if (.not. isParallel) return

    buff = val
#ifdef PARALLEL
    call MPI_allreduce(buff, val, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
#endif

  end subroutine pp_sum_d1

  !---------------------- double precision array ----------------------!

  subroutine pp_sum_dn(val,n)

    implicit none

    integer,                        intent(in)    :: n
    double precision, dimension(n), intent(inout) :: val
    double precision, dimension(n)                :: buff

    if (.not. isParallel) return

    buff(1:n) = val(1:n)
#ifdef PARALLEL
    call MPI_allreduce(buff, val, n, MPI_DOUBLE_PRECISION, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
#endif

  end subroutine pp_sum_dn


end module parallel
