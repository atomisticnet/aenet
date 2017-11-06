!-----------------------------------------------------------------------
!        tsanalyze.f90 -- a tool to analyze training set files
!-----------------------------------------------------------------------
!+ This file is part of the AENET package.
!+
!+ Copyright (C) 2012-2017 Nongnuch Artrith and Alexander Urban
!+
!+ This program is free software: you can redistribute it and/or modify
!+ it under the terms of the GNU General Public License as published by
!+ the Free Software Foundation, either version 3 of the License, or
!+ (at your option) any later version.
!+
!+ This program is distributed in the hope that it will be useful, but
!+ WITHOUT ANY WARRANTY; without even the implied warranty of
!+ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!+ General Public License for more details.
!+
!+ You should have received a copy of the GNU General Public License
!+ along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
