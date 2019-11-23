!-----------------------------------------------------------------------
!       predict.f90 - predict atomic energies of input structure
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
! 2011-11-17 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

program predict

  use aeio,      only: aeio_header,                    &
                       aeio_timestamp,                 &
                       aeio_print_copyright

  use aenet,     only: aenet_init,                     &
                       aenet_final,                    &
                       aenet_atomic_energy,            &
                       aenet_atomic_energy_and_forces, &
                       aenet_convert_atom_types,       &
                       aenet_free_atom_energy,         &
                       aenet_load_potential,           &
                       aenet_print_info,               &
                       aenet_Rc_min, aenet_Rc_max,     &
                       aenet_nnb_max

  use constants, only: PI

  use geometry,  only: geo_init,                       &
                       geo_final,                      &
                       pbc,                            &
                       latticeVec,                     &
                       recLattVec,                     &
                       geo_update_bounds,              &
                       origin,                         &
                       nAtoms,                         &
                       nTypes,                         &
                       atomType,                       &
                       atomTypeName,                   &
                       cooLatt

  use input,     only: InputData,                      &
                       read_InpPredict

  use io,        only: io_adjustl

  use lclist,    only: lcl_init,                       &
                       lcl_final,                      &
                       lcl_nmax_nbdist,                &
                       lcl_nbdist_cart

  use optimize,  only: opt_init,                       &
                       opt_final,                      &
                       opt_optimize_coords

  use parallel,  only: pp_init,                        &
                       pp_final,                       &
                       pp_bcast,                       &
                       pp_bcast_coo,                   &
                       pp_print_info,                  &
                       pp_bcast_InputData,             &
                       pp_bcast_latt,                  &
                       pp_sum,                         &
                       ppMaster, ppRank, ppSize

  implicit none

  !--------------------------------------------------------------------!
  ! A '*' in front of the variable name means that it is a broadcasted !
  ! variable and has the same value on each process.  A '+' means that !
  ! an array is allocated on all parallel processes, but does not      !
  ! necessarily have the same contents.                                !
  !                                                                    !
  !----------------------------- general ------------------------------!
  ! inp             structure with input data                          !
  ! inFile          name of the input file                             !
  !                                                                    !
  !---------------------------- structures ----------------------------!
  !*nFiles          number of input files/structures                   !
  ! cooFile         file name of structure file (atomic coordinates)   !
  !                                                                    !
  !------------------------------ output ------------------------------!
  ! Ecoh            cohesive energy of the current structure           !
  ! Etot            total energy                                       !
  !+forCart         cartesian atomic forces of the current structure   !
  !+atomicEnergy    energy of the individual atoms in the structure    !
  !--------------------------------------------------------------------!

  type(InputData)                               :: inp

  character(len=1024)                           :: inFile, cooFile, strucFile

  integer,          dimension(:),   allocatable :: atomType_orig

  integer                                       :: istruc, nStrucs

  double precision                              :: Ecoh, Etot, E0
  double precision, dimension(:,:), allocatable :: forCart
  double precision, dimension(:),   allocatable :: atomicEnergy

  double precision, dimension(3)                :: F_mav, F_max, F_avg
  double precision                              :: F_rms, F_rms_prev
  double precision, dimension(3)                :: dmax
  integer                                       :: imax

  integer                                       :: iter, stat
  logical                                       :: conv


  !-------------------------- initialization --------------------------!

  call initialize(inFile, strucFile, inp)

  if (ppMaster .and. (inp%verbosity > 0)) call aenet_print_info()

  ! number of structures to calculate
  if (len_trim(strucFile)>0) then
     nStrucs = 1
  else
     if (inp%nStrucs <= 0) then
        if (ppMaster) then
           write(0,*) "Error: no input structures specified in ", &
                      trim(inFile)
        end if
        call finalize()
        stop
     else
        nStrucs = inp%nStrucs
     end if
  end if

  !--------- loop over all the structures from the input file ---------!

  if (ppMaster .and. (inp%verbosity > 0)) then
     call aeio_header('Energy evaluation')
     write(*,*)
  end if

  do istruc = 1, nStrucs

     iter = 0

     ! only the master process reads the input structure:
     if (ppMaster) then
        if (len_trim(strucFile)>0) then
           cooFile = strucFile
        else
           cooFile = inp%strucFile(istruc)
        end if
        call geo_init(cooFile, 'xsf')
        ! convert atom type IDs to ANN potential IDs
        allocate(atomType_orig(nAtoms))
        atomType_orig = atomType
        call aenet_convert_atom_types(&
             atomTypeName, atomType_orig, atomType, stat)
     end if
     call pp_bcast(nAtoms)
     call pp_bcast(nTypes)
     call pp_bcast(pbc)
     call pp_bcast_latt(latticeVec)
     call pp_bcast_latt(recLattVec)

     ! allocate memory for other parallel processes and transfer data
     if (.not. ppMaster) then
        allocate(cooLatt(3,nAtoms),   &
                 atomType(nAtoms))
     end if
     if (inp%do_forces) then
        allocate(forCart(3,nAtoms))
     end if
     if (inp%print_atomic_energies) then
        allocate(atomicEnergy(nAtoms))
     end if
     call pp_bcast_coo(cooLatt, nAtoms)
     call pp_bcast(atomType, nAtoms)

     ! write out basic structure info
     if (ppMaster .and. (inp%verbosity > 0)) then
        call print_fileinfo(istruc, cooFile, latticeVec, nAtoms, nTypes)
     end if

     ! evaluate atomic energy and forces
     if (inp%do_forces .and. inp%print_atomic_energies) then
        call get_energy(latticeVec, nAtoms, cooLatt, atomType, pbc, &
                        Ecoh, Etot, forCart=forCart, atomicEnergy=atomicEnergy)
     else if (inp%do_forces) then
        call get_energy(latticeVec, nAtoms, cooLatt, atomType, pbc, &
                        Ecoh, Etot, forCart=forCart)
     else if (inp%print_atomic_energies) then
        call get_energy(latticeVec, nAtoms, cooLatt, atomType, pbc, &
                        Ecoh, Etot, atomicEnergy=atomicEnergy)
     else
        call get_energy(latticeVec, nAtoms, cooLatt, atomType, pbc, &
                        Ecoh, Etot)
     end if

     if (ppMaster) then
        if (inp%do_forces .and. inp%print_atomic_energies) then
           call print_coordinates(iter, latticeVec, nAtoms, nTypes, cooLatt, &
                                  atomType_orig, atomTypeName, origin, &
                                  forCart=forCart, atomicEnergy=atomicEnergy)
        else if (inp%do_forces) then
           call print_coordinates(iter, latticeVec, nAtoms, nTypes, cooLatt, &
                                  atomType_orig, atomTypeName, origin, &
                                  forCart=forCart)
        else if (inp%print_atomic_energies) then
           call print_coordinates(iter, latticeVec, nAtoms, nTypes, cooLatt, &
                                  atomType_orig, atomTypeName, origin, &
                                  atomicEnergy=atomicEnergy)
        else
           call print_coordinates(iter, latticeVec, nAtoms, nTypes, cooLatt, &
                                  atomType_orig, atomTypeName, origin)
        end if
     end if

     ! optimize coordinates, if requested:
     if (inp%do_relax .and. inp%do_forces) then
        E0 = Ecoh
        if (ppMaster) then
           call opt_init(inp%relax_method, 3*nAtoms, &
                         ftol=inp%relax_E_conv, gtol=inp%relax_F_conv)
           if (inp%verbosity > 0) then
              write(*,*) "Geometry optimization:"
              write(*,*)
              write(*,*) '       ', '   energy change  ', '     rms force    '
              write(*,*) '       ', '        (eV)      ', '     (eV/Ang)     '
              write(*,*) repeat('-',60)
           end if
        end if
        dmax = matmul(recLattVec, inp%relax_dmax)/(2.0d0*PI)
        conv = .false.
        relax : do while(iter < inp%relax_steps)
           iter = iter + 1
           if (ppMaster) then
              call calc_rms_force(forCart, F_mav, F_max, imax, F_avg, F_rms)
              if (inp%verbosity > 0) then
                 if (iter>1) then
                    write(*,'(1x,I5,2x,F16.8,2x,F16.8,2x,F16.8)') &
                         iter, Ecoh-E0, F_rms, F_rms - F_rms_prev
                 else
                    write(*,'(1x,I5,2x,F16.8,2x,F16.8)') iter, Ecoh-E0, F_rms
                 end if
              end if
              F_rms_prev = F_rms
              call opt_optimize_coords(Ecoh, nAtoms, cooLatt, forCart, &
                                       conv, dmax=(/0.1d0, 0.1d0, 0.1d0/))
              if (.not. pbc) then
                 call geo_update_bounds(cooLatt, latticeVec, recLattVec, origin)
              end if
           end if
           call pp_bcast(conv)
           if (conv) then
              if (ppMaster .and. (inp%verbosity > 0)) then
                 write(*,*) "   converged after ", &
                            trim(io_adjustl(iter)), " iterations."
              end if
              exit relax
           end if
           call pp_bcast_coo(cooLatt, nAtoms)
           call get_energy(latticeVec, nAtoms, cooLatt, atomType, pbc, &
                           Ecoh, Etot, forCart=forCart)
        end do relax
        if (ppMaster) then
           call opt_final()
           write(*,*)
           call print_coordinates(iter, latticeVec, nAtoms, nTypes, cooLatt, &
                                  atomType_orig, atomTypeName, origin, forCart)
        end if
     end if

     if (ppMaster) then
        if (inp%do_forces) then
           call print_energy(Ecoh, Etot, forCart)
        else
           call print_energy(Ecoh, Etot)
        end if
        call geo_final()
        deallocate(atomType_orig)
     else
        deallocate(cooLatt, atomType)
     end if
     if (inp%do_forces) then
        deallocate(forCart)
     end if
     if (inp%print_atomic_energies) then
        deallocate(atomicEnergy)
     end if

     if (ppMaster .and. istruc < nStrucs) then
        write(*,*) repeat('+', 70)
        write(*,*)
     end if

  end do

  !--------------------------- finalization ---------------------------!

  call finalize()

contains

  subroutine initialize(inFile, strucFile, inp)

    implicit none

    character(len=*), intent(out) :: inFile
    character(len=*), intent(out) :: strucFile
    type(InputData),  intent(out) :: inp


    logical :: fexists
    integer :: nargs
    integer :: stat
    integer :: itype

    call pp_init()

    if (ppMaster) then
       nargs = command_argument_count()
       if (nargs < 1) then
          write(0,*) "Error: No input file specified."
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

       ! read name of structure from command line, if present
       if (nargs > 1) then
          call get_command_argument(2, value=strucFile)
       else
          strucFile = ''
       end if

       ! read general input file
       inp = read_InpPredict(inFile)
    end if
    call pp_bcast(inFile)
    call pp_bcast(strucFile)
    call pp_bcast_InputData(inp)

    if (inp%verbosity > 0) call pp_print_info()

    ! initialize aenet
    call aenet_init(inp%typeName, stat)
    if (stat /= 0) then
       write(0,*) 'Error: aenet initialization failed'
       call finalize()
       stop
    end if

    ! load ANN potentials
    do itype = 1, inp%nTypes
       call aenet_load_potential(itype, inp%netFile(itype), stat)
       if (stat /= 0) then
       write(0,*) 'Error: could not load ANN potentials'
          call finalize()
          stop
       end if
    end do

    if (ppMaster .and. (inp%verbosity > 0)) then
       ! write header and copyright info
       call aeio_header("Atomic Energy Network Interpolation", char='=')
       call aeio_header(aeio_timestamp(), char=' ')
       write(*,*)
       call aeio_print_copyright('2015-2018', 'Nongnuch Artrith and Alexander Urban')
    end if

  end subroutine initialize

  !--------------------------------------------------------------------!

  subroutine finalize()

    implicit none

    integer :: stat

    if (allocated(atomType_orig)) deallocate(atomType_orig)

    if (ppMaster .and. (inp%verbosity > 0)) then
       call aeio_header(aeio_timestamp(), char=' ')
       call aeio_header("Atomic Energy Network done.", char='=')
    end if

    call aenet_final(stat)
    call pp_final()

 end subroutine finalize

  !--------------------------------------------------------------------!

  subroutine print_usage()

    implicit none

    write(*,*)
    write(*,*) "predict.x -- Predict/interpolate atomic energy."
    write(*,'(1x,70("-"))')
    write(*,*) 'Usage: predict.x <input-file> [<structure files>]'
    write(*,*)
    write(*,*) 'See the documentation or the source code for a description of the'
    write(*,*) 'input file format.  Structure files can either be listed in the'
    write(*,*) 'input file, or specified on the command line.'
    write(*,*)

  end subroutine print_usage

  !--------------------------------------------------------------------!
  !                           general output                           !
  !--------------------------------------------------------------------!

  subroutine print_fileinfo(istruc, file, latticeVec, nAtoms, nTypes)

    implicit none

    integer,                                         intent(in) :: istruc
    character(len=*),                                intent(in) :: file
    double precision, dimension(3,3),                intent(in) :: latticeVec
    integer,                                         intent(in) :: nAtoms
    integer,                                         intent(in) :: nTypes

    write(*,*) 'Structure number  : ', trim(io_adjustl(istruc))
    write(*,*) 'File name         : ', trim(adjustl(file))
    write(*,*) 'Number of atoms   : ', trim(io_adjustl(nAtoms))
    write(*,*) 'Number of species : ', trim(io_adjustl(nTypes))
    write(*,*)

    write(*,*) 'Lattice vectors (Angstrom):'
    write(*,*)
    write(*,'(3x,"a = ( ",3(2x,F15.8)," )")') latticeVec(1:3,1)
    write(*,'(3x,"b = ( ",3(2x,F15.8)," )")') latticeVec(1:3,2)
    write(*,'(3x,"c = ( ",3(2x,F15.8)," )")') latticeVec(1:3,3)
    write(*,*)

  end subroutine print_fileinfo

  !--------------------------------------------------------------------!

  subroutine print_coordinates(iter, latticeVec, nAtoms, nTypes, cooLatt, &
                               atomType, atomTypeName, origin, forCart, &
                               atomicEnergy)

    implicit none

    integer,                                         intent(in) :: iter
    double precision, dimension(3,3),                intent(in) :: latticeVec
    integer,                                         intent(in) :: nAtoms
    integer,                                         intent(in) :: nTypes
    double precision, dimension(3,nAtoms),           intent(in) :: cooLatt
    integer,          dimension(nAtoms),             intent(in) :: atomType
    character(len=*), dimension(nTypes),             intent(in) :: atomTypeName
    double precision, dimension(3),                  intent(in) :: origin
    double precision, dimension(3,nAtoms), optional, intent(in) :: forCart
    double precision, dimension(nAtoms),   optional, intent(in) :: atomicEnergy

    character(len=80)              :: header
    integer                        :: iat
    character(len=2)               :: symbol
    double precision, dimension(3) :: cooCart

    header = 'Cartesian atomic coordinates'
    if (iter == 0) then
       header = trim(header) // ' (input)'
    else
       header = trim(header) // ' (optimized)'
    end if
    if (present(forCart)) then
       header = trim(header) // ' and corresponding atomic forces'
    end if
    header = trim(header) // ':'

    write(*,*) trim(header)
    write(*,*)
    write(*,'(1x,2x,3(2x,A12))', advance='no') &
         '     x      ', '     y      ', '     z      '
    if (present(forCart)) then
       write(*,'(3(2x,A12))', advance='no') &
            '     Fx     ', '    Fy      ', '    Fz      '
    end if
    if (present(atomicEnergy)) then
       write(*,'(2x,A12)', advance='no') '   E_atom   '
    end if
    write(*,*)
    write(*,'(1x,2x,3(2x,A12))', advance='no') &
         '    (Ang)   ', '    (Ang)   ', '    (Ang)   '
    if (present(forCart)) then
       write(*,'(3(2x,A12))', advance='no') &
            '  (eV/Ang)  ', '  (eV/Ang)  ', '  (eV/Ang)  '
    end if
    if (present(atomicEnergy)) then
       write(*,'(2x,A12)', advance='no') '    (eV)    '
    end if
    write(*,*)
    write(*,'(1x,44("-"))', advance='no')
    if (present(forCart)) write(*,'(42("-"))', advance='no')
    if (present(atomicEnergy)) write(*,'(14("-"))', advance='no')
    write(*,*)
    do iat = 1, nAtoms
       symbol = atomTypeName(atomType(iat))
       cooCart(1:3) = matmul(latticeVec, cooLatt(1:3,iat)) + origin
       write(*,'(1x,A2,3(2x,F12.6))', advance='no') symbol, cooCart(1:3)
       if (present(forCart)) then
          write(*,'(3(2x,F12.6))', advance='no') forCart(1:3,iat)
       end if
       if (present(atomicEnergy)) then
          write(*,'(2x,F12.6)', advance='no') atomicEnergy(iat)
       end if
       write(*,*)
    end do
    write(*,*)

  end subroutine print_coordinates

  !--------------------------------------------------------------------!

  subroutine print_energy(Ecoh, Etot, forCart)

    implicit none

    double precision,                           intent(in) :: Ecoh, Etot
    double precision, dimension(:,:), optional, intent(in) :: forCart

    double precision, dimension(3) :: F_mav, F_max, F_avg
    double precision               :: F_rms
    integer                        :: imax

    write(*,'(1x,"Cohesive energy            :",2x,F20.8," eV")') Ecoh
    write(*,'(1x,"Total energy               :",2x,F20.8," eV")') Etot
    if (present(forCart)) then
       call calc_rms_force(forCart, F_mav, F_max, imax, F_avg, F_rms)
       write(*,'(1x,"Mean force (must be zero)  :",3(2x,F12.6))') F_avg(1:3)
       write(*,'(1x,"Mean absolute force        :",3(2x,F12.6))') F_mav(1:3)
       write(*,'(1x,"Maximum force              :",3(2x,F12.6))') F_max(1:3)
       write(*,'(1x,"RMS force                  :",2x,F12.6)')    F_rms
       write(*,*) 'The maximum force is acting on atom ', &
            trim(io_adjustl(imax)), '.'
       write(*,*) 'All forces are given in eV/Angstrom.'
    end if
    write(*,*)

  end subroutine print_energy

  !--------------------------------------------------------------------!
  !                 calculate atomic energy and forces                 !
  !--------------------------------------------------------------------!

  subroutine get_energy(latticeVec, nAtoms, cooLatt, atomType, &
                        pbc, Ecoh, Etot, forCart, atomicEnergy)

    implicit none

    double precision, dimension(3,3),                intent(in)  :: latticeVec
    integer,                                         intent(in)  :: nAtoms
    double precision, dimension(3,nAtoms),           intent(in)  :: cooLatt
    integer,          dimension(nAtoms),             intent(in)  :: atomType
    logical,                                         intent(in)  :: pbc
    double precision,                                intent(out) :: Ecoh
    double precision,                                intent(out) :: Etot
    double precision, dimension(3,nAtoms), optional, intent(out) :: forCart
    double precision, dimension(nAtoms),   optional, intent(out) :: atomicEnergy

#ifdef CHECK_FORCES
    double precision, dimension(3,nAtoms) :: forCart_num
    double precision                              :: E_i1, E_i2
    double precision :: d
    double precision, dimension(3,3) :: dd
    integer :: i, j
#endif

    logical                                       :: do_F, do_E_atom

    integer                                       :: nnb
    integer,          dimension(aenet_nnb_max)    :: nblist
    double precision, dimension(3,aenet_nnb_max)  :: nbcoo
    double precision, dimension(aenet_nnb_max)    :: nbdist
    integer,          dimension(aenet_nnb_max)    :: nbtype

    integer                                       :: type_i
    double precision, dimension(3)                :: coo_i
    double precision                              :: E_i

    integer                                       :: iatom, stat

    do_F = present(forCart)
    do_E_atom = present(atomicEnergy)
    if (do_F) forCart(1:3,1:nAtoms) = 0.0d0
    if (do_E_atom) atomicEnergy(1:nAtoms) = 0.0d0

    call lcl_init(aenet_Rc_min, aenet_Rc_max, latticeVec, nAtoms, &
                  atomType, cooLatt, pbc)

#ifdef CHECK_FORCES
    d = 0.01d0
    dd(:,1) = [d, 0.0d0, 0.0d0]
    dd(:,2) = [0.0d0, d, 0.0d0]
    dd(:,3) = [0.0d0, 0.0d0, d]
    forCart_num = 0.0d0
#endif

    Ecoh = 0.0d0
    Etot = 0.0d0
    atoms : do iatom = 1, nAtoms

       ! distribute atoms over processes:
       if (mod(iatom-1,ppSize) /= ppRank) cycle

       type_i = atomType(iatom)
       coo_i(1:3) = matmul(latticeVec, cooLatt(1:3,iatom))

       ! get all atoms of species type_i within the cut-off:
       nnb = aenet_nnb_max
       call lcl_nbdist_cart(iatom, nnb, nbcoo, nbdist, aenet_Rc_max, &
                            nblist=nblist, nbtype=nbtype)

       if (do_F) then
          call aenet_atomic_energy_and_forces( &
               coo_i, type_i, iatom, nnb, nbcoo, nbtype, nblist, &
               nAtoms, E_i, forCart, stat)
#ifdef CHECK_FORCES
          do i = 1, 3
             coo_i = coo_i - dd(:,i)
             call aenet_atomic_energy(coo_i, type_i, nnb, nbcoo, nbtype, &
                                   E_i1, stat)
             coo_i = coo_i + 2.0d0*dd(:,i)
             call aenet_atomic_energy(coo_i, type_i, nnb, nbcoo, nbtype, &
                                   E_i2, stat)
             coo_i = coo_i - dd(:,i)
             forCart_num(i,iatom) = forCart_num(i,iatom) - (E_i2 - E_i1)/(2.0d0*d)
          end do
          do j = 1, nnb
             do i = 1, 3
                nbcoo(:,j) = nbcoo(:,j) - dd(:,i)
                call aenet_atomic_energy(coo_i, type_i, nnb, nbcoo, nbtype, &
                                         E_i1, stat)
                nbcoo(:,j) = nbcoo(:,j) + 2.0d0*dd(:,i)
                call aenet_atomic_energy(coo_i, type_i, nnb, nbcoo, nbtype, &
                                         E_i2, stat)
                nbcoo(:,j) = nbcoo(:,j) - dd(:,i)
                forCart_num(i,nblist(j)) = forCart_num(i,nblist(j)) - (E_i2 - E_i1)/(2.0d0*d)
             end do
          end do
#endif
       else
          call aenet_atomic_energy(coo_i, type_i, nnb, nbcoo, nbtype, &
                                   E_i, stat)
       end if

       Etot = Etot + E_i
       Ecoh = Ecoh + E_i - aenet_free_atom_energy(type_i)
       if (do_E_atom) atomicEnergy(iatom) = E_i

    end do atoms

#ifdef CHECK_FORCES
    open(99, file='CHECK_FORCES.dat', status='replace', action='write')
    do iatom = 1, nAtoms
       write(99,'(9(1x,ES15.8))') &
            forCart(1:3,iatom), forCart_num(1:3,iatom), &
            forCart(1:3,iatom) - forCart_num(1:3,iatom)
    end do
    close(99)
#endif

    call lcl_final()

    call pp_sum(Ecoh)
    call pp_sum(Etot)

    ! gather results from all processes
    if (do_F) then
       if (ppSize>1) then
          do iatom = 1, nAtoms
             call pp_sum(forCart(1:3,iatom), 3)
          end do
       end if
    end if
    if (do_E_atom) call pp_sum(atomicEnergy(1:nAtoms), nAtoms)

  end subroutine get_energy

  !--------------------------------------------------------------------!
  !                       analyze atomic forces                        !
  !--------------------------------------------------------------------!

  subroutine calc_rms_force(forCart, F_mav, F_max, imax, F_avg, F_rms)

    implicit none

    double precision, dimension(:,:), optional, intent(in)  :: forCart
    double precision, dimension(3),             intent(out) :: F_mav
    double precision, dimension(3),             intent(out) :: F_max
    integer,                                    intent(out) :: imax
    double precision, dimension(3),             intent(out) :: F_avg
    double precision,                           intent(out) :: F_rms

    integer                        :: nAtoms
    double precision               :: F_abs2, F_abs2_max
    integer                        :: iat

    nAtoms = size(forCart(1,:))
    F_rms       = 0.0d0
    F_mav(1:3)  = 0.0d0
    F_max(1:3)  = 0.0d0
    F_avg(1:3)  = 0.0d0
    F_abs2      = 0.0d0
    F_abs2_max  = 0.0d0
    do iat = 1, nAtoms
       F_avg(1:3) = F_avg(1:3) + forCart(1:3,iat)
       F_mav(1:3) = F_mav(1:3) + abs(forCart(1:3,iat))
       F_abs2 = sum(forCart(1:3,iat)*forCart(1:3,iat))
       F_rms  = F_rms + F_abs2
       if (F_abs2 > F_abs2_max) then
          F_abs2_max  = F_abs2
          F_max(1:3) = forCart(1:3,iat)
          imax       = iat
       end if
    end do
    F_avg(1:3) = F_avg(1:3)/dble(nAtoms)
    F_mav(1:3) = F_mav(1:3)/dble(nAtoms)
    F_rms = sqrt(F_rms/dble(nAtoms))

  end subroutine calc_rms_force

end program predict
