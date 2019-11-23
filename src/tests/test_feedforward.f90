!-----------------------------------------------------------------------
!             Unit tests for the feedforward module
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
! 2014-09-28 Alexander Urban (AU) and Nongnuch Artrith (NA)
!-----------------------------------------------------------------------
program test_feedforward

  use feedforward, only: Network,                &
                         new_Network,            &
                         save_Network,           &
                         save_Network_ASCII,     &
                         load_Network,           &
                         load_Network_ASCII,     &
                         del_Network,            &
                         ff_activate,            &
                         ff_get_nweights,        &
                         ff_print_info,          &
                         ff_eval,                &
                         ff_deriv,               &
                         ff_wderiv,              &
                         ff_random_init_weights, &
                         ff_change_activation

  use io,          only: io_unlink
  use unittest,    only: tst_new, tst_check_passed, tst_equal

  implicit none

  call test_IO()
  call test_eval()
  call test_activation()

contains

  subroutine test_IO()

    implicit none

    type(Network) :: nn1, nn2
    logical       :: has_passed

    call tst_new("Feedforward Test 1: save and restore")
    has_passed = .true.

    ! Set up new ANN with random weights
    nn1 = new_Network((/3, 5, 5, 5, 3/))
    call ff_random_init_weights(nn1)

    call save_Network(nn1, 'TEST_SAVE_BINARY')
    call save_Network_ASCII(nn1, 'TEST_SAVE_ASCII')

    ! Load ANN from binary file and compare with original one
    nn2 = load_Network('TEST_SAVE_BINARY')
    has_passed = (has_passed .and. all(nn2%f_a == nn1%f_a))
    has_passed = (has_passed .and. all(nn2%iw == nn1%iw))
    has_passed = (has_passed .and. all(nn2%iv == nn1%iv))
    has_passed = (has_passed .and. all(nn2%W == nn1%W))
    call del_network(nn2)

    ! Load ANN from ASCII file and compare to original one
    nn2 = load_Network_ASCII('TEST_SAVE_ASCII')
    has_passed = (has_passed .and. all(nn2%f_a == nn1%f_a))
    has_passed = (has_passed .and. all(nn2%iw == nn1%iw))
    has_passed = (has_passed .and. all(nn2%iv == nn1%iv))
    has_passed = (has_passed .and. all(nn2%W == nn1%W))
    call del_network(nn2)

    call del_Network(nn1)

    call io_unlink('TEST_SAVE_BINARY')
    call io_unlink('TEST_SAVE_ASCII')

    call tst_check_passed(has_passed)

  end subroutine test_IO

  !--------------------------------------------------------------------!

  subroutine test_eval()

    implicit none

    integer, parameter :: nx = 3
    integer, parameter :: ny = 3

    double precision, dimension(nx)               :: x
    double precision, dimension(ny)               :: y, y1, y2
    double precision, dimension(ny,nx)            :: dy_dx, dy_dx_num
    double precision, dimension(:,:), allocatable :: dy_dw, dy_dw_num
    integer                                       :: nw
    double precision                              :: d

    type(Network) :: nn
    logical       :: has_passed

    integer :: nvalues, nweights
    double precision, dimension(:), allocatable :: values
    double precision, dimension(:), allocatable :: derivs
    double precision, dimension(:), allocatable :: jacobian

    integer :: i, j

    call tst_new("Feedforward Test 2: evaluation and derivative")
    has_passed = .true.

    ! Set up new ANN with random weights
    nn = new_Network((/3, 5, 5, 5, 3/))
    call ff_random_init_weights(nn)

    ! allocate memory
    nvalues = nn%nvalues
    nweights = nn%Wsize
    nw = ff_get_nweights(nn)
    allocate(values(nvalues), derivs(nvalues), jacobian(nweights), &
             dy_dw(ny,nw), dy_dw_num(ny,nw))

    x = [1.0d0, 2.0d0, 3.0d0]
    call ff_eval(nn, nx, x, ny, values, derivs, y)
    call ff_deriv(nn, nx, ny, derivs, jacobian, dy_dx)
    call ff_wderiv(nn, nw, ny, values, derivs, jacobian, dy_dw)

    ! numerical derivative dy/dx
    d = 0.01d0
    do i = 1, nx
       x(i) = x(i) - d
       call ff_eval(nn, nx, x, ny, values, derivs, y1)
       x(i) = x(i) + 2.0d0*d
       call ff_eval(nn, nx, x, ny, values, derivs, y2)
       x(i) = x(i) - d
       dy_dx_num(1:3,i) = (y2 - y1)/(2.0d0*d)
    end do
    do j = 1, nx
       do i = 1, ny
          has_passed = tst_equal(dy_dx(i,j), dy_dx_num(i,j), prec=0.05d0)
       end do
    end do

    ! numerical derivative dy/dw
    d = 0.01d0
    do i = 1, nw
       nn%W(i) = nn%W(i) - d
       call ff_eval(nn, nx, x, ny, values, derivs, y1)
       nn%W(i) = nn%W(i) + 2.0d0*d
       call ff_eval(nn, nx, x, ny, values, derivs, y2)
       nn%W(i) = nn%W(i) - d
       dy_dw_num(1:3,i) = (y2 - y1)/(2.0d0*d)
    end do
    open(99, file='TEST_dy_dw.dat', status='replace', action='write')
    do j = 1, nw
       do i = 1, ny
          has_passed = tst_equal(dy_dw(i,j), dy_dw_num(i,j), prec=0.05d0)
          if (.not. has_passed) then
             write(*,*) dy_dw(i,j), dy_dw_num(i,j), dy_dw(i,j)-dy_dw_num(i,j)
             stop
          end if
       end do
       write(99,'(9(1x,ES24.17))') &
                   dy_dw(:,j), dy_dw_num(:,j), dy_dw(:,j)-dy_dw_num(:,j)
    end do
    close(99)

    deallocate(values, derivs, jacobian, dy_dw, dy_dw_num)

    call tst_check_passed(has_passed)
    if (has_passed) then
       call io_unlink('TEST_dy_dw.dat')
    else
       write(*,*) 'see file: TEST_dy_dw.dat'
    end if

  end subroutine test_eval

  !----------------------- activation functions -----------------------!

  subroutine test_activation()

    implicit none

    integer :: t
    double precision :: d
    double precision :: x0, x1
    double precision, dimension(:),   allocatable :: x
    double precision, dimension(:,:), allocatable :: y, dy
    integer :: i, N

    logical       :: has_passed

    call tst_new("Feedforward Test 3: activation functions")
    has_passed = .true.

    d  = 0.01d0
    x0 = -2.0d0
    x1 =  2.0d0
    N = ceiling((x1-x0)/d + 1)
    allocate(x(N), y(N,0:4), dy(N,0:4))
    x(1) = x0
    do i = 2, N
       x(i) = x(i-1) + d
    end do

    ! Function types
    !   0 : linear function f(x) = x
    !   1 : hyperbolic tangent, y in [-1:1]
    !   2 : sigmoid,            y in [ 0:1]
    !   3 : modified tanh,      y in [-1.7159:1.7159]  f(+/-1) = +/-1
    !   4 : tanh & linear twisting term
    ftype : do t = 0, 4
       call ff_activate(t, x, y(:,t), dy(:,t))
       ! numerical derivative
       assert : do i = 1, N-2, 2
          x0 = y(i,t)
          x1 = y(i+2,t)
          has_passed = (has_passed .and. &
                        tst_equal(dy(i+1,t), (x1-x0)/(2.0d0*d), prec=1.0d-3))
          if (.not. has_passed) then
             write(*,'(A,I2,1x)', advance='no') 'assertion failed for t = ', t
             exit assert
          end if
       end do assert
    end do ftype

    call tst_check_passed(has_passed)
    if (.not. has_passed) then
       ! write out function values and derivatives
       write(*,*) 'Check file: TEST_ACTIVATION'
       open(99, file='TEST_ACTIVATION', status='replace', action='write')
       do i = 1, N
          write(99,'(11(1x,ES15.8))') &
             x(i), y(i,0), dy(i,0), y(i,1), dy(i,1), y(i,2), dy(i,2), &
             y(i,3), dy(i,3), y(i,4), dy(i,4)
       end do
       close(99)
    end if

    deallocate(x, y, dy)

  end subroutine test_activation

end program test_feedforward
