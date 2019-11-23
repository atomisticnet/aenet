!-----------------------------------------------------------------------
!             sortlib.f90 - Library with sorting routines
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
! sort(A)        : sort double precision array A
! sort(A,n)      : sort first n elements in double precision array A
! argsort(A,idx) : return index idx that sorts dp array A
!-----------------------------------------------------------------------
! 2010-02-11 Nongnuch Artrith (NA), Alexander Urban (AU)
!-----------------------------------------------------------------------

module sortlib

  implicit none

  public  :: sort,           &
             sort_i,         &
             sort_d,         &
             argsort,        &
             argsort_d,      &
             argsort_i

  private :: qsort_d,        &
             qsort_i,        &
             qargsort_d,     &
             qargsort_i,     &
             partition_d,    &
             partition_i,    &
             argpartition_d, &
             argpartition_i

  interface sort
     module procedure sort_i, sort_d
  end interface
  interface argsort
     module procedure argsort_i, argsort_d
  end interface

contains

  !--------------------------------------------------------------------!
  !                          User interfaces                           !
  !--------------------------------------------------------------------!

  subroutine sort_d(A, n)

    implicit none

    integer, optional,              intent(in)    :: n
    double precision, dimension(:), intent(inout) :: A
    integer                                       :: m

    if (present(n)) then
       m = n
    else
       m = size(A)
    end if

    call qsort_d(A, m)

  end subroutine sort_d

  !--------------------------------------------------------------------!

  subroutine sort_i(A, n)

    implicit none

    integer, optional,     intent(in)    :: n
    integer, dimension(:), intent(inout) :: A
    integer                              :: m

    if (present(n)) then
       m = n
    else
       m = size(A)
    end if

    call qsort_i(A, m)

  end subroutine sort_i

  !--------------------------------------------------------------------!

  subroutine argsort_d(A, idx)

    implicit none

    double precision, dimension(:), intent(in)  :: A
    integer,          dimension(:), intent(out) :: idx
    integer                                     :: i, n

    n = size(idx)
    do i = 1, n
       idx(i) = i
    end do

    call qargsort_d(A, idx, n, n)

  end subroutine argsort_d

  !--------------------------------------------------------------------!

  subroutine argsort_i(A, idx)

    implicit none

    integer, dimension(:), intent(in)  :: A
    integer, dimension(:), intent(out) :: idx
    integer                            :: i, n

    n = size(idx)
    do i = 1, n
       idx(i) = i
    end do

    call qargsort_i(A, idx, n, n)

  end subroutine argsort_i

  !--------------------------------------------------------------------!
  !                      quicksort implementation                      !
  !--------------------------------------------------------------------!

  recursive subroutine qsort_d(A,n)

    implicit none

    integer,                        intent(in)    :: n
    double precision, dimension(n), intent(inout) :: A
    integer                                       :: iq

    if(n > 1) then
       call partition_d(A, n, iq)
       call qsort_d(A(:iq-1), iq-1)
       call qsort_d(A(iq:), n-iq+1)
    endif

  end subroutine qsort_d

  !--------------------------------------------------------------------!

  recursive subroutine qsort_i(A,n)

    implicit none

    integer,               intent(in)    :: n
    integer, dimension(n), intent(inout) :: A
    integer                              :: iq

    if(n > 1) then
       call partition_i(A, n, iq)
       call qsort_i(A(:iq-1), iq-1)
       call qsort_i(A(iq:), n-iq+1)
    endif

  end subroutine qsort_i

  !--------------------------------------------------------------------!

  subroutine partition_d(A, n, marker)

    implicit none

    integer,                        intent(in)    :: n
    double precision, dimension(n), intent(inout) :: A
    integer,                        intent(out)   :: marker
    integer                                       :: i, j
    double precision                              :: temp
    double precision                              :: x    ! pivot point

    x = A(1)
    i = 0
    j = n + 1

    do
       j = j-1
       do while (A(j) > x)
          j = j-1
       end do
       i = i+1
       do while (A(i) < x)
          i = i+1
       end do
       if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
       else if (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do

  end subroutine partition_d

  !--------------------------------------------------------------------!

  subroutine partition_i(A, n, marker)

    implicit none

    integer,               intent(in)    :: n
    integer, dimension(n), intent(inout) :: A
    integer,               intent(out)   :: marker
    integer                              :: i, j
    integer                              :: temp
    double precision                     :: x    ! pivot point

    x = A(1)
    i = 0
    j = n + 1

    do
       j = j-1
       do while (A(j) > x)
          j = j-1
       end do
       i = i+1
       do while (A(i) < x)
          i = i+1
       end do
       if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
       else if (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do

  end subroutine partition_i

  !--------------------------------------------------------------------!
  !                   Argument sort using quicksort                    !
  !--------------------------------------------------------------------!

  recursive subroutine qargsort_d(A, idx, n, m)

    implicit none

    integer,                        intent(in)    :: n
    integer,                        intent(in)    :: m
    double precision, dimension(n), intent(in)    :: A
    integer,          dimension(m), intent(inout) :: idx
    integer                                       :: iq

    if(m > 1) then
       call argpartition_d(A, idx, n, m, iq)
       call qargsort_d(A, idx(:iq-1), n, iq-1)
       call qargsort_d(A, idx(iq:),   n, m-iq+1)
    endif

  end subroutine qargsort_d

  !--------------------------------------------------------------------!

  recursive subroutine qargsort_i(A, idx, n, m)

    implicit none

    integer,               intent(in)    :: n
    integer,               intent(in)    :: m
    integer, dimension(n), intent(in)    :: A
    integer, dimension(m), intent(inout) :: idx
    integer                              :: iq

    if(m > 1) then
       call argpartition_i(A, idx, n, m, iq)
       call qargsort_i(A, idx(:iq-1), n, iq-1)
       call qargsort_i(A, idx(iq:),   n, m-iq+1)
    endif

  end subroutine qargsort_i

  !--------------------------------------------------------------------!

  subroutine argpartition_d(A, idx, n, m, marker)

    implicit none

    integer,                        intent(in)    :: n
    integer,                        intent(in)    :: m
    double precision, dimension(n), intent(in)    :: A
    integer,          dimension(m), intent(inout) :: idx
    integer,                        intent(out)   :: marker
    integer                                       :: i, j
    integer                                       :: temp
    double precision                              :: x    ! pivot point

    x = A(idx(1))
    i = 0
    j = m + 1

    do
       j = j-1
       do while (A(idx(j)) > x)
          j = j-1
       end do
       i = i+1
       do while (A(idx(i)) < x)
          i = i+1
       end do
       if (i < j) then
          ! exchange A(idx(i)) and A(idx(j))
          temp   = idx(i)
          idx(i) = idx(j)
          idx(j) = temp
       else if (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do

  end subroutine argpartition_d

  !--------------------------------------------------------------------!

  subroutine argpartition_i(A, idx, n, m, marker)

    implicit none

    integer,               intent(in)    :: n
    integer,               intent(in)    :: m
    integer, dimension(n), intent(in)    :: A
    integer, dimension(m), intent(inout) :: idx
    integer,               intent(out)   :: marker
    integer                              :: i, j
    integer                              :: temp
    double precision                     :: x    ! pivot point

    x = A(idx(1))
    i = 0
    j = m + 1

    do
       j = j-1
       do while (A(idx(j)) > x)
          j = j-1
       end do
       i = i+1
       do while (A(idx(i)) < x)
          i = i+1
       end do
       if (i < j) then
          ! exchange A(idx(i)) and A(idx(j))
          temp   = idx(i)
          idx(i) = idx(j)
          idx(j) = temp
       else if (i == j) then
          marker = i+1
          return
       else
          marker = i
          return
       endif
    end do

  end subroutine argpartition_i

end module sortlib
