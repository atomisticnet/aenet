!-----------------------------------------------------------------------
!      generate.f90 - generate training sets for use with train.x
!-----------------------------------------------------------------------
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2018 Nongnuch Artrith and Alexander Urban
!+
!+ This Source Code Form is subject to the terms of the Mozilla Public
!+ License, v. 2.0. If a copy of the MPL was not distributed with this
!+ file, You can obtain one at http://mozilla.org/MPL/2.0/.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!+ Mozilla Public License, v. 2.0, for more details.
!-----------------------------------------------------------------------
! 2011-10-19 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

program generate

  use aeio,     only: aeio_readline,        &
                      aeio_header,          &
                      aeio_timestamp,       &
                      aeio_print_copyright, &
                      PATHLEN, LINELEN

  use geometry, only: geo_init,          &
                      geo_final,         &
                      geo_itype_of_name, &
                      geo_type_conv,     &
                      pbc,               &
                      latticeVec,        &
                      nAtoms,            &
                      nTypes,            &
                      atomType,          &
                      atomTypeName,      &
                      cooLatt,           &
                      cooCart,           &
                      forCart,           &
                      hasEnergy,         &
                      hasForces,         &
                      cohesiveEnergy,    &
                      totalEnergy

  use input,    only: InputData,         &
                      read_InpGenerate,  &
                      del_InputData

  use io,       only: io_adjustl,        &
                      io_center,         &
                      io_lower,          &
                      io_readnext,       &
                      io_unit

  use lclist,   only: lcl_init,          &
                      lcl_final,         &
                      lcl_print_info,    &
                      lcl_nmax_nbdist,   &
                      lcl_nbdist_cart

  use sfsetup,  only: Setup,                 &
                      read_Setup_parameters, &
                      save_Setup,            &
                      del_Setup,             &
                      stp_init,              &
                      stp_final,             &
                      stp_get_range,         &
                      stp_print_info,        &
                      stp_eval,              &
                      nsf_max

  use timing,   only: tng_init,          &
                      tng_final,         &
                      tng_timing,        &
                      tng_timing2,       &
                      tng_timing3,       &
                      tng_dump

  use trainset, only: TrnSet,                 &
                      new_TrnSet,             &
                      close_TrnSet,           &
                      ts_print_info,          &
                      ts_write_header,        &
                      ts_write_sf_info,       &
                      ts_write_atom_info,     &
                      ts_write_structure_info

  implicit none

  !--------------------------------------------------------------------!
  ! stp(i)         structural fingerprint basis setup for atom type i  !
  ! r_min, r_max   lower and upper bound for atomic interactions       !
  ! ts             training set reference                              !
  !                                                                    !
  ! nnb_max, nnb   max. and actual number of neighboring atoms         !
  ! nbcoo(i,j)     i-th component of the coordinates of the j-th       !
  !                neighboring atom                                    !
  ! nbdist(i)      distance of the i-th neighbor                       !
  !                                                                    !
  ! sfval(i)         value of the i-th basis function                  !
  ! sfderiv_i(i,j)   i-th component of the derivative of the j-th SF   !
  !                  with respect to the central atom                  !
  !                  sfderiv_i(3,nsf_max)                              !
  ! sfderiv_j(i,j,k) i-th component of the derivative of the j-th SF   !
  !                  with respect to the coordinates of atom k         !
  !                  sfderiv_j(3,nsf_max,nnb_max)                      !
  !                                                                    !
  ! E_coh          cohesive energy                                     !
  ! nFiles_inv     = 1/inp%nStrucs                                     !
  !                                                                    !
  ! inFile         name of the input file for the generate.x program   !
  ! cooFile        name of the currently active structure file         !
  ! keyword        the last keyword read from the input file           !
  !                                                                    !
  ! do_debug       if .true., additional files containing debugging    !
  !                info will be created                                !
  !                                                                    !
  ! u_*            file units                                          !
  !--------------------------------------------------------------------!

  type(InputData)                                :: inp

  type(Setup),       dimension(:),   allocatable :: stp
  double precision                               :: r_min, r_max
  type(TrnSet)                                   :: ts

  integer                                        :: nnb_max, nnb
  double precision,  dimension(:,:), allocatable :: nbcoo
  double precision,  dimension(:),   allocatable :: nbdist
  integer,           dimension(:),   allocatable :: nbtype

  double precision, dimension(:),     allocatable :: sfval
  double precision, dimension(:,:),   allocatable :: sfderiv_i
  double precision, dimension(:,:,:), allocatable :: sfderiv_j

  double precision                               :: E_coh
  integer                                        :: ifile
  double precision                               :: nFiles_inv

  character(len=PATHLEN)                         :: inFile
  character(len=PATHLEN)                         :: cooFile
  character(len=LINELEN)                         :: keyword

  integer                                        :: itype1
  integer                                        :: itype, iatom

  integer                                        :: iline
  character(len=1024)                            :: line

  integer                                        :: u_in, u_tng
  logical                                        :: do_debug = .false.
  integer                                        :: u_dbg, idbg

  integer :: i

  ! timing registers
  integer, parameter :: R_GEO = 1, R_NBL = 2, R_SF = 3

  !-------------------------- initialization --------------------------!

  call initialize(inFile)

  inp = read_InpGenerate(inFile)
  allocate(stp(inp%nTypes))
  call load_symmfunc_setups(inp, stp)

  ! call parse_input(inFile)

  if (inp%do_timing) then
     u_tng = io_unit()
     call tng_init(unit=u_tng, file='generate.time', registers=3)
     write(*,*) 'Timing info will be written to: ', 'generate.time'
     write(*,*)
  end if
  if (do_debug) then
     u_dbg = io_unit()
     open(u_dbg, file='generate.debug', status='replace', action='write')
  end if

  ! get interaction range and max. number of atoms within range
  call stp_get_range(inp%nTypes, stp, r_min, r_max)
  nnb_max = lcl_nmax_nbdist(r_min, r_max)
  allocate(nbcoo(3,nnb_max), nbdist(nnb_max), nbtype(nnb_max))

  ! initialize workspace for structural fingerprint basis:
  call stp_init(inp%nTypes, stp, nnb_max)
  if (inp%do_timing) call tng_timing('Structural fingerprint basis initialized.')

  ! allocate workspace for basis function evaluation:
  allocate(sfval(nsf_max), sfderiv_i(3,nsf_max), sfderiv_j(3,nsf_max,nnb_max))
  sfval(:) = 0.0d0
  sfderiv_i(:,:) = 0.0d0
  sfderiv_j(:,:,:) = 0.0d0

  call aeio_header('Generation of training set started')
  write(*,*)

  write(*,*) 'Number of atom types  : ', trim(io_adjustl(inp%nTypes))
  write(*,'(1x,"types                 : ")', advance='no')
  do itype = 1, inp%nTypes
     if (mod(itype,7) == 0) write(*,'(29x)')
     write(*,'(A5,1x)', advance='no') inp%typeName(itype)
  end do
  write(*,*)
  write(*,*) "Number of structures  : ", trim(io_adjustl(inp%nStrucs))
  write(*,*)

  !-------------- write basis function settings to stdout -------------!

  call aeio_header("Structural fingerprint basis set-up")
  write(*,*)

  do itype1 = 1, inp%nTypes
     call stp_print_info(stp(itype1))
  end do

  !----------- write training set header to the output file -----------!

  ts = new_TrnSet(inp%nTypes, inp%typeName, inp%atomicEnergy, &
                  inp%nStrucs, trim(inp%outFileName))

  if (inp%do_timing) call tng_timing('Training set file started.')

  !------------------ iterate over coordinates files ------------------!

  call aeio_header("Adding structures to the training set")
  write(*,*)

  u_in = io_unit()
  open(u_in, file=inFile, status='old', action='read')
  rewind(u_in)

  iline = 0
  do
     ! forward until the FILES keyword:
     call aeio_readline(u_in, iline, line)
     read(line,*) keyword
     if (trim(keyword) == 'FILES') then
        read(u_in,*)
        exit
     end if
  end do

  ! header for stdout
  write(*,'("#",A6,2x,A6,2x,A6,2x,A15,2x,A)') &
       'N', 'nAtoms', 'nTypes', 'E/atom', 'structure file (xsf)'

  nFiles_inv = 1.0d0/dble(inp%nStrucs)
  structures : do ifile = 1, inp%nStrucs

     if (inp%do_timing) call tng_timing('Structure: '// io_adjustl(ifile))

     call aeio_readline(u_in, iline, line)
     cooFile = trim(line)

     call geo_init(cooFile, 'xsf')
     if (inp%do_timing) call tng_timing3(register=R_GEO)
     if (.not. (hasForces .and. hasEnergy)) then
        write(0,*) ">>>", hasForces, hasEnergy
        write(0,*) "Error: incomplete output data in : ", trim(cooFile)
        call finalize()
        stop
     end if

     if (nTypes > inp%nTypes) then
        write(*,*) 'Skipping ', trim(adjustl(cooFile)), &
                   ': too many atomic species'
        call geo_final()
        cycle structures
     end if

     if (abs(cohesiveEnergy) /= 0.0d0) then
        E_coh = cohesiveEnergy
     else
        ! if only the total energy is available, we have to calculate
        ! the cohesive energy at this point
        E_coh = totalEnergy
        do iatom = 1, nAtoms
           itype1 = atomType(iatom)
           itype1 = geo_type_conv(itype1, nTypes, atomTypeName, &
                                  inp%nTypes, inp%typeName)
           E_coh = E_coh - inp%atomicEnergy(itype1)
        end do
     end if

     write(*,'(1x,I6,2x,I6,2x,I6,2x,ES15.8,2x,A)') &
          ifile, nAtoms, nTypes, E_coh/dble(nAtoms), &
          trim(adjustl(cooFile))

     call lcl_init(r_min, r_max, latticeVec, nAtoms, atomType, cooLatt, pbc)
     if (inp%do_timing) call tng_timing3(register=R_NBL)

     ! write structure info (atoms, types, energy) to training set file:
     call ts_write_structure_info(ts, cooFile, nAtoms, nTypes, E_coh)

     atoms : do iatom = 1, nAtoms

        ! determine the training atom type of atom `iatom' in global
        ! index terms
        itype1 = atomType(iatom)
        itype1 = geo_type_conv(itype1, nTypes, atomTypeName, &
                               inp%nTypes, inp%typeName)

        ! assert that atom type is included in the set-ups:
        if (itype1 == 0) then
           write(0,*) "Error: not a valid structure    : ", trim(cooFile)
           write(0,*) "       Additional species found."
           call finalize()
           stop
        end if

        ! write atom info (species, forces) to training set file:
        call ts_write_atom_info(ts, itype1, cooCart(iatom), forCart(iatom))

        ! get all atoms within cut-off:
        nnb = nnb_max
        call lcl_nbdist_cart(iatom, nnb, nbcoo, nbdist, r_cut=r_max, nbtype=nbtype)
        if (inp%do_timing) call tng_timing3(register=R_NBL)
        write(*,'(1x,I6,2x,A2,2x,I6)') &
             iatom, trim(atomTypeName(atomType(iatom))), nnb

        ! convert atom types to global index:
        do i = 1, nnb
           nbtype(i) = geo_type_conv(nbtype(i), nTypes, atomTypeName, &
                                    inp%nTypes, inp%typeName)
           if (nbtype(i) == 0) then
              write(0,*) "Error: atom type not found in setup."
              call finalize()
              stop
           end if
        end do

        ! evaluate the structural fingerprint basis function set-up:
        call stp_eval(itype1, cooCart(iatom), nnb, nbcoo, nbtype, &
                      stp(itype1), sfval=sfval)

        if (do_debug) then
           do idbg = 1, stp(itype1)%nsf
              write(u_dbg,'(1x,ES15.8,1x)', advance='no') sfval(idbg)
           end do
           write(u_dbg,*)
        end if

        if (inp%do_timing) call tng_timing3(register=R_SF)

        ! write basis function values and derivatives
        ! to the training set file:
        call ts_write_sf_info(ts, stp(itype1)%nsf, sfval(1:stp(itype1)%nsf))

     end do atoms

     if (inp%do_timing) then
        call tng_dump(R_GEO, 'time spent reading geometries (so far)')
        call tng_dump(R_NBL, 'time spent in the neighbor list (so far)')
        call tng_dump(R_SF,  'time spent evaluating structural fingerprints (so far)')
     end if

     call lcl_final()
     call geo_final()

  end do structures
  write(*,*)

  if (inp%do_timing) then
     call tng_timing('Loop over structures done.')
     call tng_dump(R_GEO, 'total time spent reading geometries')
     call tng_dump(R_NBL, 'total time spent in the neighbor list')
     call tng_dump(R_SF,  'total time spent evaluating structural fingerprints')
  end if

  !--------- save basis function setups with final statistics ---------!

  call ts_print_info(ts)

  !----------------------------- finalize -----------------------------!

  deallocate(nbcoo, nbdist, nbtype)
  close(u_in)

  call close_TrnSet(ts, stp=stp(1:inp%nTypes))
  call finalize()


contains !=============================================================!


  subroutine initialize(inFile)

    implicit none

    character(len=*), intent(out) :: inFile

    integer :: nargs
    logical :: fexists

    call aeio_header("generate.x - training set generation", char='=')
    write(*,*)

    call aeio_print_copyright('2015-2018', 'Nongnuch Artrith and Alexander Urban')

    nargs = command_argument_count()
    if (nargs < 1) then
       write(0,*) "Error: No input file provided."
       call print_usage()
       call finalize()
       stop
    end if

    call get_command_argument(1, value=inFile)
    inquire(file=trim(inFile), exist=fexists)
    if (.not. fexists) then
       write(0,*) "Error: File not found: ", trim(inFile)
       call print_usage()
       call finalize()
       stop
    end if

  end subroutine initialize

  !--------------------------------------------------------------------!

  subroutine finalize()

    implicit none

    integer :: itype

    if (allocated(sfval)) then
       deallocate(sfval, sfderiv_i, sfderiv_j)
    end if

    if (allocated(stp)) then
       call stp_final(inp%nTypes, stp)
       do itype = 1, inp%nTypes
          call del_Setup(stp(itype))
       end do
       deallocate(stp, inp%typeName, inp%atomicEnergy)
    end if

    if (ts%init) call close_TrnSet(ts)

    if (allocated(nbcoo)) deallocate(nbcoo, nbdist)

    if (inp%do_timing) call tng_final()
    if (do_debug)  close(u_dbg)

    call aeio_header(aeio_timestamp(), char=' ')
    call aeio_header("Training set generation done.", char='=')

  end subroutine finalize

  !--------------------------------------------------------------------!

  subroutine print_usage()

    implicit none

    write(*,*)
    write(*,*) "generate.x -- Generate training sets for use with `train.x'"
    write(*,'(1x,70("-"))')
    write(*,*) 'Usage: generate.x <input-file>'
    write(*,*)
    write(*,*) 'See the documentation or the source code for a description of the '
    write(*,*) 'input file format.'
    write(*,*)

  end subroutine print_usage

  !--------------------------------------------------------------------!

  subroutine load_symmfunc_setups(inp, stp)

    implicit none

    type(InputData),           intent(in)  :: inp
    type(Setup), dimension(:), intent(out) :: stp

    integer :: i

    do i = 1, inp%nTypes
       stp(i) = read_Setup_parameters(inp%setupFile(i), inp%typeName(:))
       if (.not. (trim(stp(i)%atomtype) == trim(inp%typeName(i)))) then
          write(0,*) "Error: Inconsistent atom type in setup:"
          write(0,*) "       type expected : ", trim(inp%typeName(i))
          write(0,*) "       type found    : ", trim(stp(i)%atomtype)
          call finalize()
          stop
       end if
    end do

  end subroutine load_symmfunc_setups

end program generate
