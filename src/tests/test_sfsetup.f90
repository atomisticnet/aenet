!-----------------------------------------------------------------------
!       Unit tests for the structural fingerprint setup module
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
! 2014-09-29 Alexander Urban (AU) and Nongnuch Artrith (NA)
!-----------------------------------------------------------------------
program test_sfsetup

  use io,       only: io_unlink
  use unittest, only: tst_new, tst_check_passed, tst_equal
  use sfsetup,  only: Setup,                 &
                     read_Setup_parameters, &
                     save_Setup,            &
                     save_Setup_ASCII,      &
                     load_Setup,            &
                     load_Setup_ASCII,      &
                     del_Setup

  implicit none

  call test_IO()

contains

  subroutine test_IO()

    implicit none

    character(len=2), dimension(2), parameter :: global_types = ['H ', 'O ']

    type(Setup) :: stp1, stp2

    logical :: has_passed

    !--------------------------- chebyshev ----------------------------!

    call tst_new("Structural Fingerprint Test 1: save and restore (Chebyshev)")
    has_passed = .true.

    call write_setup_params_chebyshev('TEST_SFSETUP_01')

    stp1 = read_Setup_parameters('TEST_SFSETUP_01', global_types)
    call save_Setup(stp1, file='TEST_SETUP_BINARY_01')
    call save_Setup_ASCII(stp1, file='TEST_SETUP_ASCII_01')

    stp2 = load_Setup(global_types, file='TEST_SETUP_BINARY_01')
    has_passed = (has_passed .and. tst_equal(stp1%description, stp2%description))
    has_passed = (has_passed .and. tst_equal(stp1%atomtype, stp2%atomtype))
    has_passed = (has_passed .and. tst_equal(stp1%nenv, stp2%nenv))
    has_passed = (has_passed .and. tst_equal(stp1%envtypes, stp2%envtypes))
    has_passed = (has_passed .and. tst_equal(stp1%Rc_min, stp2%Rc_min, prec=1.0d-3))
    has_passed = (has_passed .and. tst_equal(stp1%Rc_max, stp2%Rc_max, prec=1.0d-3))
    has_passed = (has_passed .and. tst_equal(stp1%sftype, stp2%sftype))
    has_passed = (has_passed .and. tst_equal(stp1%nsf, stp2%nsf))
    has_passed = (has_passed .and. tst_equal(stp1%nsfparam, stp2%nsfparam))
    has_passed = (has_passed .and. tst_equal(stp1%sf, stp2%sf))
    has_passed = (has_passed .and. tst_equal(stp1%sfparam, stp2%sfparam, prec=1.0d-6))
    has_passed = (has_passed .and. tst_equal(stp1%sfenv, stp2%sfenv))
    if (.not. has_passed) write(*,*) 'Load from binary file differs.'
    call del_Setup(stp2)

    stp2 = load_Setup_ASCII(global_types, file='TEST_SETUP_ASCII_01')
    has_passed = (has_passed .and. tst_equal(stp1%description, stp2%description))
    has_passed = (has_passed .and. tst_equal(stp1%atomtype, stp2%atomtype))
    has_passed = (has_passed .and. tst_equal(stp1%nenv, stp2%nenv))
    has_passed = (has_passed .and. tst_equal(stp1%envtypes, stp2%envtypes))
    has_passed = (has_passed .and. tst_equal(stp1%Rc_min, stp2%Rc_min, prec=1.0d-3))
    has_passed = (has_passed .and. tst_equal(stp1%Rc_max, stp2%Rc_max, prec=1.0d-3))
    has_passed = (has_passed .and. tst_equal(stp1%sftype, stp2%sftype))
    has_passed = (has_passed .and. tst_equal(stp1%nsf, stp2%nsf))
    has_passed = (has_passed .and. tst_equal(stp1%nsfparam, stp2%nsfparam))
    has_passed = (has_passed .and. tst_equal(stp1%sf, stp2%sf))
    has_passed = (has_passed .and. tst_equal(stp1%sfparam, stp2%sfparam, prec=1.0d-6))
    has_passed = (has_passed .and. tst_equal(stp1%sfenv, stp2%sfenv))
    if (.not. has_passed) then
       write(*,*) 'Load from ASCII file differs. Compare "TEST_SETUP_ASCII2_01"'
       call save_Setup_ASCII(stp2, file="TEST_SETUP_ASCII2_01")
    end if
    call del_Setup(stp2)

    call del_Setup(stp1)
    call tst_check_passed(has_passed)

    if (has_passed) then
       call io_unlink('TEST_SFSETUP_01')
       call io_unlink('TEST_SETUP_BINARY_01')
       call io_unlink('TEST_SETUP_ASCII_01')
    end if

    !--------------------------- behler2011 ---------------------------!

    call tst_new("Structural Fingerprint Test 2: save and restore (Behler2011)")
    has_passed = .true.

    call write_setup_params_Behler2011('TEST_SFSETUP_02')

    stp1 = read_Setup_parameters('TEST_SFSETUP_02', global_types)
    call save_Setup(stp1, file='TEST_SETUP_BINARY_02')
    call save_Setup_ASCII(stp1, file='TEST_SETUP_ASCII_02')

    stp2 = load_Setup(global_types, file='TEST_SETUP_BINARY_02')
    has_passed = (has_passed .and. tst_equal(stp1%description, stp2%description))
    has_passed = (has_passed .and. tst_equal(stp1%atomtype, stp2%atomtype))
    has_passed = (has_passed .and. tst_equal(stp1%nenv, stp2%nenv))
    has_passed = (has_passed .and. tst_equal(stp1%envtypes, stp2%envtypes))
    has_passed = (has_passed .and. tst_equal(stp1%Rc_min, stp2%Rc_min, prec=1.0d-3))
    has_passed = (has_passed .and. tst_equal(stp1%Rc_max, stp2%Rc_max, prec=1.0d-3))
    has_passed = (has_passed .and. tst_equal(stp1%sftype, stp2%sftype))
    has_passed = (has_passed .and. tst_equal(stp1%nsf, stp2%nsf))
    has_passed = (has_passed .and. tst_equal(stp1%nsfparam, stp2%nsfparam))
    has_passed = (has_passed .and. tst_equal(stp1%sf, stp2%sf))
    has_passed = (has_passed .and. tst_equal(stp1%sfparam, stp2%sfparam, prec=1.0d-6))
    has_passed = (has_passed .and. tst_equal(stp1%sfenv, stp2%sfenv))
    if (.not. has_passed) write(*,*) 'Load from binary file differs.'
    call del_Setup(stp2)

    stp2 = load_Setup_ASCII(global_types, file='TEST_SETUP_ASCII_02')
    has_passed = (has_passed .and. tst_equal(stp1%description, stp2%description))
    has_passed = (has_passed .and. tst_equal(stp1%atomtype, stp2%atomtype))
    has_passed = (has_passed .and. tst_equal(stp1%nenv, stp2%nenv))
    has_passed = (has_passed .and. tst_equal(stp1%envtypes, stp2%envtypes))
    has_passed = (has_passed .and. tst_equal(stp1%Rc_min, stp2%Rc_min, prec=1.0d-3))
    has_passed = (has_passed .and. tst_equal(stp1%Rc_max, stp2%Rc_max, prec=1.0d-3))
    has_passed = (has_passed .and. tst_equal(stp1%sftype, stp2%sftype))
    has_passed = (has_passed .and. tst_equal(stp1%nsf, stp2%nsf))
    has_passed = (has_passed .and. tst_equal(stp1%nsfparam, stp2%nsfparam))
    has_passed = (has_passed .and. tst_equal(stp1%sf, stp2%sf))
    has_passed = (has_passed .and. tst_equal(stp1%sfparam, stp2%sfparam, prec=1.0d-6))
    has_passed = (has_passed .and. tst_equal(stp1%sfenv, stp2%sfenv))
    if (.not. has_passed) then
       write(*,*) 'Load from ASCII file differs. Compare "TEST_SETUP_ASCII2_02"'
       call save_Setup_ASCII(stp2, file="TEST_SETUP_ASCII2_02")
    end if
    call del_Setup(stp2)

    call del_Setup(stp1)
    call tst_check_passed(has_passed)

    if (has_passed) then
       call io_unlink('TEST_SFSETUP_02')
       call io_unlink('TEST_SETUP_BINARY_02')
       call io_unlink('TEST_SETUP_ASCII_02')
    end if

  end subroutine test_IO

  !--------------------------------------------------------------------!

  subroutine write_setup_params_chebyshev(file)

    implicit none

    character(len=*), intent(in) :: file
    integer, parameter :: u = 99

    open(u, file=trim(file), status='replace', action='write')

    write(u,'(A)') 'DESCR'
    write(u,'(A)') 'Just for testing.'
    write(u,'(A)') 'END DESCR'
    write(u,'(A)')
    write(u,'(A)') 'ATOM O'
    write(u,'(A)')
    write(u,'(A)') 'ENV  2'
    write(u,'(A)') 'O'
    write(u,'(A)') 'H'
    write(u,'(A)')
    write(u,'(A)') 'RMIN 0.75d0'
    write(u,'(A)')
    write(u,'(A)') 'BASIS type=chebyshev'
    write(u,'(A)') 'radial_Rc=6.5 radial_N=10 angular_Rc=6.5 angular_N=10'

    close(u)

  end subroutine write_setup_params_chebyshev

  !--------------------------------------------------------------------!

  subroutine write_setup_params_Behler2011(file)

    implicit none

    character(len=*), intent(in) :: file
    integer, parameter :: u = 99

    open(u, file=trim(file), status='replace', action='write')

    write(u,'(A)') 'DESCR'
    write(u,'(A)') 'Symmetry function setup for Oxygen'
    write(u,'(A)') '(just for testing)'
    write(u,'(A)') 'END DESCR'
    write(u,'(A)')
    write(u,'(A)') 'ATOM O'
    write(u,'(A)')
    write(u,'(A)') 'ENV  2'
    write(u,'(A)') 'O'
    write(u,'(A)') 'H'
    write(u,'(A)')
    write(u,'(A)') 'RMIN 0.75d0'
    write(u,'(A)')
    write(u,'(A)') 'SYMMFUNC type=Behler2011'
    write(u,'(A)') '52'
    write(u,'(A)') 'G=2 type2=O   eta=0.003214  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=H   eta=0.003214  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=O   eta=0.035711  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=H   eta=0.035711  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=O   eta=0.071421  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=H   eta=0.071421  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=O   eta=0.124987  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=H   eta=0.124987  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=O   eta=0.214264  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=H   eta=0.214264  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=O   eta=0.357106  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=H   eta=0.357106  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=O   eta=0.714213  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=H   eta=0.714213  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=O   eta=1.428426  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=2 type2=H   eta=1.428426  Rs=0.0000  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.000357 lambda= -1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.000357 lambda= -1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.000357 lambda= -1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.028569 lambda= -1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.028569 lambda= -1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.028569 lambda= -1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.089277 lambda= -1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.089277 lambda= -1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.089277 lambda= -1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.000357 lambda=  1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.000357 lambda=  1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.000357 lambda=  1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.028569 lambda=  1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.028569 lambda=  1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.028569 lambda=  1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.089277 lambda=  1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.089277 lambda=  1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.089277 lambda=  1.0  zeta= 1.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.000357 lambda= -1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.000357 lambda= -1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.000357 lambda= -1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.028569 lambda= -1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.028569 lambda= -1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.028569 lambda= -1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.089277 lambda= -1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.089277 lambda= -1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.089277 lambda= -1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.000357 lambda=  1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.000357 lambda=  1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.000357 lambda=  1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.028569 lambda=  1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.028569 lambda=  1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.028569 lambda=  1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=O    eta=0.089277 lambda=  1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=O  type3=H    eta=0.089277 lambda=  1.0  zeta= 2.0  Rc=6.5000'
    write(u,'(A)') 'G=4 type2=H  type3=H    eta=0.089277 lambda=  1.0  zeta= 2.0  Rc=6.5000'

    close(u)

  end subroutine write_setup_params_Behler2011

end program test_sfsetup
