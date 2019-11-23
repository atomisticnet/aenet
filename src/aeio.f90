!-----------------------------------------------------------------------
!                 aeio.f90 - Atomic Energy I/O module
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
! 2011-10-20 Alexander Urban (AU), Nongnuch Artrith (NA)
!-----------------------------------------------------------------------

module aeio

  use io, only: io_lower,    &
                io_readval,  &
                io_adjustl,  &
                io_center

  implicit none

  public :: aeio_readline,             &
            aeio_header,               &
            aeio_timestamp,            &
            aeio_print_copyright,      &
            aeio_assert_file_exists,   &
            aeio_assert_file_notexists

  private :: aeio_readline_c1, aeio_readline_cn, &
             aeio_readline_i1, aeio_readline_in, &
             aeio_readline_d1, aeio_readline_dn, &
             aeio_readline_l1, aeio_readline_ln

  !---------------------------- constants -----------------------------!
  ! LINELEN    max. length of a line read from an input file           !
  ! PATHLEN    max. number of characters available for a file path     !
  ! TYPELEN    max. number of characters for atmoic species names      !
  ! STDIN      file unit of standard in                                !
  ! STDOUT     file unit of standard out                               !
  ! STDERR     file unit of standard error                             !
  !--------------------------------------------------------------------!

  integer, parameter, public :: LINELEN = 1024
  integer, parameter, public :: PATHLEN = 1024
  integer, parameter, public :: TYPELEN = 2
  integer, parameter, public :: STDIN   = 5
  integer, parameter, public :: STDOUT  = 6
  integer, parameter, public :: STDERR  = 0

  !--------------------------------------------------------------------!
  !   aeio_readline() - read next line with contents from input file   !
  !                                                                    !
  ! A line is skipped if                                               !
  !                                                                    !
  ! - it only contains blanks                                          !
  ! - the first non-blank character is `!', `#', or `%'                !
  !                                                                    !
  ! usage: call aeio_readline(unit, iline, line[, n][, stat])          !
  !        implementations available for 'line' as character, integer, !
  !        double precision, and logical (also array of length 'n')    !
  !--------------------------------------------------------------------!

  interface aeio_readline
     module procedure aeio_readline_c1, aeio_readline_cn, &
                      aeio_readline_i1, aeio_readline_in, &
                      aeio_readline_d1, aeio_readline_dn, &
                      aeio_readline_l1, aeio_readline_ln
  end interface aeio_readline

contains

  !--------------------------------------------------------------------!
  !           write a centered header (for formatted output)           !
  !--------------------------------------------------------------------!

  subroutine aeio_header(str, char, unit)

    character(len=*),    intent(in) :: str
    character, optional, intent(in) :: char
    integer,   optional, intent(in) :: unit

    character :: c

    if (present(char)) then
       c = char
    else
       c = '-'
    end if

    if (present(unit)) then
       write(unit,*) repeat(c,70)
       write(unit,*) io_center(trim(str),70)
       write(unit,*) repeat(c,70)
    else
       write(*,*) repeat(c,70)
       write(*,*) io_center(trim(str),70)
       write(*,*) repeat(c,70)
    end if

  end subroutine aeio_header

  !--------------------------------------------------------------------!
  !              return a formatted date and time string               !
  !--------------------------------------------------------------------!

  function aeio_timestamp() result(date)

    implicit none

    character(len=30)     :: date
    integer, dimension(8) :: v

    call date_and_time(values=v)

    write(date, '(I4.4,"-",I2.2,"-",I2.2,2x,I2.2,":",I2.2,":",I2.2)') &
         v(1:3), v(5:7)

  end function aeio_timestamp

  !--------------------------------------------------------------------!
  !                        auxiliary procedures                        !
  !--------------------------------------------------------------------!

  subroutine aeio_assert_file_exists(file)
    implicit none
    character(len=*), intent(in) :: file
    logical :: fexists
    inquire(file=trim(adjustl(file)), exist=fexists)
    if (.not. fexists) then
       write(0,*) "Error: file not found: ", trim(adjustl(file))
       stop
    end if
  end subroutine aeio_assert_file_exists

  !--------------------------------------------------------------------!

  subroutine aeio_assert_file_notexists(file)
    implicit none
    character(len=*), intent(in) :: file
    logical :: fexists
    inquire(file=trim(adjustl(file)), exist=fexists)
    if (fexists) then
       write(0,*) "Error: file already exists: ", trim(adjustl(file))
       stop
    end if
  end subroutine aeio_assert_file_notexists

  !--------------------------------------------------------------------!
  !                        print copyright info                        !
  !--------------------------------------------------------------------!

  subroutine aeio_print_copyright(year, authors)

    implicit none

    character(len=*), intent(in) :: year
    character(len=*), intent(in) :: authors

    write(*,*) 'Copyright (C) ', trim(adjustl(year)), ' ', trim(adjustl(authors))
    write(*,*)
    write(*,*) "This program is distributed in the hope that it will be useful, but"
    write(*,*) "WITHOUT ANY WARRANTY; without even the implied warranty of"
    write(*,*) "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
    write(*,*) "Mozilla Public License, v. 2.0, for more details."
    write(*,*)

  end subroutine aeio_print_copyright

  !--------------------------------------------------------------------!
  !     Implementation of aeio_readline() for different data types     !
  !--------------------------------------------------------------------!

  subroutine aeio_readline_c1(u_in, iline, line, stat)

    implicit none

    integer,           intent(in)    :: u_in
    integer,           intent(inout) :: iline
    character(len=*),  intent(out)   :: line
    integer, optional, intent(out)   :: stat

    integer :: stat2

    stat2 = 0
    do
       read(u_in, '(A)', iostat=stat2) line
       if (stat2 == 0) then
          iline = iline + 1
          line  = trim(adjustl(line))
          if (line(1:1) == '!')    cycle
          if (line(1:1) == '#')    cycle
          if (line(1:1) == '%')    cycle
          if (len_trim(line) == 0) cycle
       end if
       exit
    end do
    if (present(stat)) stat = stat2

  end subroutine aeio_readline_c1

  !--------------------------------------------------------------------!

  subroutine aeio_readline_cn(u_in, iline, line, n, stat)

    implicit none

    integer,                        intent(in)    :: u_in
    integer,                        intent(inout) :: iline
    integer,                        intent(in)    :: n
    character(len=*), dimension(n), intent(out)   :: line
    integer, optional,              intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call aeio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line(1:n)
    end if
    if (present(stat)) stat = stat2

  end subroutine aeio_readline_cn

  !--------------------------------------------------------------------!

  subroutine aeio_readline_i1(u_in, iline, line, stat)

    implicit none

    integer,           intent(in)    :: u_in
    integer,           intent(inout) :: iline
    integer,           intent(out)   :: line
    integer, optional, intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call aeio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line
    end if
    if (present(stat)) stat = stat2

  end subroutine aeio_readline_i1

  !--------------------------------------------------------------------!

  subroutine aeio_readline_in(u_in, iline, line, n, stat)

    implicit none

    integer,               intent(in)    :: u_in
    integer,               intent(inout) :: iline
    integer,               intent(in)    :: n
    integer, dimension(n), intent(out)   :: line
    integer, optional,     intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call aeio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line(1:n)
    end if
    if (present(stat)) stat = stat2

  end subroutine aeio_readline_in

  !--------------------------------------------------------------------!

  subroutine aeio_readline_d1(u_in, iline, line, stat)

    implicit none

    integer,           intent(in)    :: u_in
    integer,           intent(inout) :: iline
    double precision,  intent(out)   :: line
    integer, optional, intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call aeio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line
    end if
    if (present(stat)) stat = stat2

  end subroutine aeio_readline_d1

  !--------------------------------------------------------------------!

  subroutine aeio_readline_dn(u_in, iline, line, n, stat)

    implicit none

    integer,                        intent(in)    :: u_in
    integer,                        intent(inout) :: iline
    integer,                        intent(in)    :: n
    double precision, dimension(n), intent(out)   :: line
    integer, optional,              intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call aeio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line(1:n)
    end if
    if (present(stat)) stat = stat2

  end subroutine aeio_readline_dn

  !--------------------------------------------------------------------!

  subroutine aeio_readline_l1(u_in, iline, line, stat)

    implicit none

    integer,           intent(in)    :: u_in
    integer,           intent(inout) :: iline
    logical,           intent(out)   :: line
    integer, optional, intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call aeio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line
    end if
    if (present(stat)) stat = stat2

  end subroutine aeio_readline_l1

  !--------------------------------------------------------------------!

  subroutine aeio_readline_ln(u_in, iline, line, n, stat)

    implicit none

    integer,                        intent(in)    :: u_in
    integer,                        intent(inout) :: iline
    integer,                        intent(in)    :: n
    logical, dimension(n),          intent(out)   :: line
    integer, optional,              intent(out)   :: stat

    character(len=LINELEN) :: line2
    integer                :: stat2

    call aeio_readline_c1(u_in, iline, line2, stat2)
    if (stat2==0) then
       read(line2, *) line(1:n)
    end if
    if (present(stat)) stat = stat2

  end subroutine aeio_readline_ln


end module aeio
