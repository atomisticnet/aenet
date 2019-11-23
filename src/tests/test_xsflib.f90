!-----------------------------------------------------------------------
!       Unit tests for the XSF structure format parser module
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
! 2014-09-23 Alexander Urban (AU) and Nongnuch Artrith (NA)
!-----------------------------------------------------------------------
program xsflibtest

  use io,      only: io_unit, io_unlink
  use unittest, only: tst_new, tst_check_passed, tst_assert, tst_equal

  use xsflib, only: xsf_init,         &
                    xsf_final,        &
                    nAtoms, nTypes,   &
                    pbc,              &
                    cooCart, forCart, &
                    atomType,         &
                    atomTypeName,     &
                    hasForces

  implicit none

  call test_isolated()
  call test_periodic()

contains

  !------------------------ isolated structure ------------------------!

  subroutine test_isolated()

    implicit none

    logical :: has_passed

    call tst_new("XSFlib Test 1 - isolated structure")
    has_passed = .true.

    call write_test1('TEST1.xsf')
    call xsf_init('TEST1.xsf')

    call tst_assert(tst_equal(nAtoms, 21), 'Number of atoms')
    call tst_assert(tst_equal(nTypes, 3), 'Number of atom types')
    call tst_assert((.not. pbc), 'Not periodic')
    call tst_assert(hasForces, 'Has atomic forces')

    ! confirm some atomic coordinates
    call tst_assert(&
         tst_equal(cooCart(1:3, 3), &
         [2.56679921d0, 4.44582609d0, 6.00000283d0]), &
         "Coordinates of atom 3")
    call tst_assert(&
         tst_equal(cooCart(1:3, 19), &
         [5.90056206d0, 3.92211093d0, 10.85186772d0]), &
         "Coordinates of atom 19")
    call tst_assert(&
         tst_equal(cooCart(1:3, 20), &
         [5.13359842d0, 4.44582609d0, 10.09578481d0]), &
         "Coordinates of atom 20")

    ! confirm some atomic forces
    call tst_assert(&
         tst_equal(forCart(1:3, 3), &
         [0.43939200d0, -0.36265900d0, 0.73827500d0]), &
         "Forces of atom 3")
    call tst_assert(&
         tst_equal(forCart(1:3, 19), &
         [-0.00156829d0, 0.00139941d0, -0.00177746d0]), &
         "Forces of atom 19")
    call tst_assert(&
         tst_equal(forCart(1:3, 20), &
         [0.08246940d0, 0.10459800d0, 0.20628100d0]), &
         "Forces of atom 20")

    ! check if species have been assigned correctly
    call tst_assert(tst_equal(atomTypeName(atomType(3)), 'Cu'), &
                    'Type of atom 3')
    call tst_assert(tst_equal(atomTypeName(atomType(19)), 'O'), &
                    'Type of atom 3')
    call tst_assert(tst_equal(atomTypeName(atomType(20)), 'C'), &
                    'Type of atom 3')

    call tst_check_passed(has_passed)
    if (.not. has_passed) write(*,*) nAtoms, nTypes, pbc, hasForces

    call xsf_final()
    call io_unlink('TEST1.xsf')

  end subroutine test_isolated

  !------------------------ periodic structure ------------------------!

  subroutine test_periodic()

    implicit none

    logical :: has_passed

    call tst_new("XSFlib Test 2 - periodic structure")
    has_passed = .true.

    call write_test2('TEST2.xsf')
    call xsf_init('TEST2.xsf')

    call tst_assert(tst_equal(nAtoms, 24), 'Number of atoms')
    call tst_assert(tst_equal(nTypes, 2), 'Number of atom types')
    call tst_assert(pbc, 'Periodic')
    call tst_assert((.not. hasForces), 'No atomic forces')

    call tst_check_passed(has_passed)
    if (.not. has_passed) write(*,*) nAtoms, nTypes, pbc, hasForces

    call xsf_final()
    call io_unlink('TEST2.xsf')

  end subroutine test_periodic

  !--------------------------------------------------------------------!

  subroutine write_test1(fname)

    implicit none

    character(len=*), intent(in) :: fname
    integer :: u

    u = io_unit()

    open(u, file=trim(fname), status='replace', action='write')

    write(u,'(A)') '# total energy = -819543.67017695 eV'
    write(u,'(A)')
    write(u,'(A)') 'ATOMS'
    write(u,'(A)') 'Cu      0.00000000     0.00000000     6.00000283     0.90450300     0.49026800     0.41789800'
    write(u,'(A)') 'Cu      1.28339961     2.22291304     6.00000283     0.58819200    -0.39665700     0.30838700'
    write(u,'(A)') 'Cu      2.56679921     4.44582609     6.00000283     0.43939200    -0.36265900     0.73827500'
    write(u,'(A)') 'Cu      2.56679921     0.00000000     6.00000283    -0.04917570     0.74709300     0.34151800'
    write(u,'(A)') 'Cu      3.85019881     2.22291304     6.00000283    -0.26087400    -0.18366700    -0.04616380'
    write(u,'(A)') 'Cu      5.13359842     4.44582609     6.00000283    -0.03463600    -0.12564200     0.45479300'
    write(u,'(A)') 'Cu      5.13359841     0.00000000     6.00000283    -0.38312400     0.51523300     0.70528000'
    write(u,'(A)') 'Cu      6.41699802     2.22291304     6.00000283    -0.16177600     0.28501300     0.46754400'
    write(u,'(A)') 'Cu      7.70039763     4.44582609     6.00000283    -0.48853100    -0.27462200     0.74704900'
    write(u,'(A)') 'Cu      1.28339961     0.74097135     8.09578481     0.47004200     0.22508900    -0.75157500'
    write(u,'(A)') 'Cu      2.56679921     2.96388439     8.09578481     0.07553830    -0.42643900    -0.60240400'
    write(u,'(A)') 'Cu      3.85019882     5.18679743     8.09578481     0.02980850    -0.29619000    -0.98916000'
    write(u,'(A)') 'Cu      3.85019881     0.74097135     8.09578481    -0.08197910     0.22486800    -0.46016200'
    write(u,'(A)') 'Cu      5.13359842     2.96388439     8.09578481     0.36858900     0.05985530     0.15191600'
    write(u,'(A)') 'Cu      6.41699803     5.18679743     8.09578481     0.16791500    -0.82856700    -0.23678600'
    write(u,'(A)') 'Cu      6.41699802     0.74097135     8.09578481    -0.27475100     0.40082100    -0.60768200'
    write(u,'(A)') 'Cu      7.70039763     2.96388439     8.09578481    -0.48023800     0.34366800    -0.45382900'
    write(u,'(A)') 'Cu      8.98379723     5.18679743     8.09578481    -0.91356400    -0.50177400    -0.39263800'
    write(u,'(A)') 'O       5.90056206     3.92211093    10.85186772    -0.00156829     0.00139941    -0.00177746'
    write(u,'(A)') 'C       5.13359842     4.44582609    10.09578481     0.08246940     0.10459800     0.20628100'
    write(u,'(A)') 'O       4.10491663     5.15195854    10.08745039     0.00333637    -0.00192313     0.00050037'

    close(u)

  end subroutine write_test1

  !--------------------------------------------------------------------!

  subroutine write_test2(fname)

    implicit none

    character(len=*), intent(in) :: fname
    integer :: u

    u = io_unit()

    open(u, file=trim(fname), status='replace', action='write')

    write(u,'(A)') '# total energy = -6970801.09717 eV'
    write(u,'(A)')
    write(u,'(A)') 'CRYSTAL'
    write(u,'(A)')
    write(u,'(A)') 'PRIMVEC'
    write(u,'(A)') '      5.72106227     0.00000000     0.00000000'
    write(u,'(A)') '      2.86052948     4.64245205     0.00000000'
    write(u,'(A)') '     -0.00000354     2.56332638    27.03158552'
    write(u,'(A)')
    write(u,'(A)') 'PRIMCOORD'
    write(u,'(A)') '24 1'
    write(u,'(A)') 'Cu      5.72106694     4.53500342     0.06603627'
    write(u,'(A)') 'Au      4.29082051     2.17312786    -0.02831821'
    write(u,'(A)') 'Au      1.43026748     2.17314610    -0.02829993'
    write(u,'(A)') 'Cu      2.86054536     4.53498601     0.06610758'
    write(u,'(A)') 'Cu      5.72107830     1.67120952     2.24031171'
    write(u,'(A)') 'Au      7.15122123     3.99010491     2.29593743'
    write(u,'(A)') 'Au      4.29090136     3.99011832     2.29595968'
    write(u,'(A)') 'Cu      2.86052387     1.67131603     2.24054325'
    write(u,'(A)') 'Cu      2.86052955     3.51987257     4.52540414'
    write(u,'(A)') 'Au      1.42989027     1.19908417     4.52635549'
    write(u,'(A)') 'Au      4.29116879     1.19908661     4.52635704'
    write(u,'(A)') 'Cu      5.72106093     3.52111793     4.52699758'
    write(u,'(A)') 'Cu      5.72105984     5.28096712     6.77803627'
    write(u,'(A)') 'Au      4.29042055     2.96017872     6.77898762'
    write(u,'(A)') 'Au      1.43063680     2.96018116     6.77898916'
    write(u,'(A)') 'Cu      2.86052894     5.28221248     6.77962971'
    write(u,'(A)') 'Cu      5.72105085     2.48817946     9.06464966'
    write(u,'(A)') 'Au      7.15119625     4.81175635     9.00898939'
    write(u,'(A)') 'Au      4.29088783     4.81179881     9.00903595'
    write(u,'(A)') 'Cu      2.86050355     2.48828253     9.06487897'
    write(u,'(A)') 'Cu      2.86050096     4.26689936    11.23882313'
    write(u,'(A)') 'Au      1.43023999     1.98625486    11.33324544'
    write(u,'(A)') 'Au      4.29075303     1.98626284    11.33328944'
    write(u,'(A)') 'Cu      5.72103165     4.26688112    11.23889146'

    close(u)

  end subroutine write_test2

end program xsflibtest
