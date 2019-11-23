!-----------------------------------------------------------------------
!        tsanalyze.f90 -- a tool to analyze training set files
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
! 2014-01-08 Nongnuch Artrith (NA), Alexander Urban (AU)
!-----------------------------------------------------------------------

program tsanalyze

  use arglib,   only: args_init,   &
                      args_final,  &
                      args_switch

  use trainset, only: TrnSet,      &
                      open_TrnSet, &


  implicit none

  logical             :: do_save_sf
  character(len=1024) :: filename_sf

  call initialize()

contains

  subroutine initialize()

    implicit none

    integer :: stat
    logical :: is_help

    call args_init('--help:-h:--in:-i:--save-sf', status=stat)
    call args_switch('--help:-h', value=is_help)
    if (is_help .or. (stat/=0) .or. (nargs==0)) then
       call print_usage()
       stop
    end if

    call args_switch('--save-sf', value=filename_sf, default='')
    do_save_sf = (len_trim(filename_sf) > 0)

  end subroutine initialize

  !--------------------------------------------------------------------!

  subroutine finalize()

    implicit none

  end subroutine finalize

  !--------------------------------------------------------------------!

  subroutine print_usage()

    implicit none

  end subroutine print_usage

end program tsanalyze
