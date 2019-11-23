!-----------------------------------------------------------------------
!    A tool to convert binary training set files to an ASCII format
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
! 2019-09-25 Alexander Urban (AU) and Nongnuch Artrith (NA)
!-----------------------------------------------------------------------
program trnset_info

  use trainset, only: open_TrnSet,       &
                      close_TrnSet,      &
                      ts_print_info,     &
                      TrnSet

  implicit none

  character(len=1024) :: infile
  type(TrnSet)        :: ts

  call initialize(infile)
  ts = open_TrnSet(infile, raw=raw)
  call ts_print_info(ts)
  call close_TrnSet(ts)

contains

  subroutine initialize(infile)

    implicit none

    character(len=*), intent(out) :: infile
    integer :: iarg, nargs
    character(len=100) :: arg

    nargs = command_argument_count()
    if (nargs < 1) then
       write(0,*) "Error: No input file provided."
       call print_usage()
       call finalize()
       stop
    end if

    infile = ' '

    iarg = 1
    do while(iarg <= nargs)
       call get_command_argument(iarg, value=arg)
       select case(trim(arg))
       case('--help', '-h')
          call print_usage()
          call finalize()
          stop
       case default
          if (len_trim(infile) == 0) then
             infile = trim(arg)
          else
             write(0,*) 'Error: Unknown argument: ', trim(arg)
             call finalize()
             stop
          end if
       end select
       iarg = iarg + 1
    end do

    if ((len(infile) == 0) .or. (len(outfile) == 0))then
       write(0,*) 'Error: No input file specified.'
       call finalize()
       stop
    end if

  end subroutine initialize

  !--------------------------------------------------------------------!

  subroutine finalize()

    implicit none

  end subroutine finalize

  !--------------------------------------------------------------------!

  subroutine print_usage()

    implicit none

    write(*,*) "Print information about training set file."
    write(*,*) "----------------------------------------------------------------------"
    write(*,*) "Usage: trnset_info.x <trnset_file>"
    write(*,*)
    write(*,*) "  <trnset_file>   Training set file generated with 'generate.x'."
    write(*,*)

  end subroutine print_usage

end program trnset_info
