!-----------------------------------------------------------------------
!                timing.f90 - simple timing procedures
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
! 2011-02-07 Nongnuch Artrith (NA), Alexander Urban (AU)
!-----------------------------------------------------------------------

module timing

  implicit none
  save

  public :: tng_init,    &
            tng_final,   &
            tng_timing,  &
            tng_timing2, &
            tng_timing3, &
            tng_dump

  double precision,            public  :: tng_t_ini
  double precision,            public  :: tng_t_prev
  double precision,            public  :: tng_t_lastlog
  double precision,            public  :: tng_t_now
  double precision,            public  :: tng_t

  !--------------------------------------------------------------------!

  character(len=*), parameter, private :: TNG_DEFAULT_FILE = 'timing.log'
  integer,          parameter, private :: TNG_DEFAULT_UNIT = 55

  !--------------------------------------------------------------------!

  integer,                                     private :: tng_nregisters = 0
  double precision, dimension(:), allocatable, private :: tng_register

  logical,                                     private :: tng_isinit = .false.

  character(len=50),                           private :: tng_file
  integer,                                     private :: u_tng = 6

contains

  subroutine tng_init(file, unit, registers)

    implicit none

    character(len=*), optional, intent(in) :: file
    integer,          optional, intent(in) :: unit
    integer,          optional, intent(in) :: registers

    if (present(file)) then
       tng_file = trim(file)
    else
       tng_file = TNG_DEFAULT_FILE
    end if

    if (present(unit)) then
       u_tng = unit
    else
       u_tng = TNG_DEFAULT_UNIT
    end if

    if (present(registers)) then
       tng_nregisters = registers
       allocate(tng_register(registers))
       tng_register(:) = 0.0d0
    end if

    open(u_tng, file=trim(tng_file), status='replace', action='write')
    call tng_timing()

  end subroutine tng_init

  !--------------------------------------------------------------------!

  subroutine tng_final()

    implicit none

    call tng_timing('Timing finished.')
    close(u_tng)

    if (allocated(tng_register)) then
       deallocate(tng_register)
       tng_nregisters = 0
    end if

    tng_isinit = .false.

  end subroutine tng_final

  !--------------------------------------------------------------------!
  !           Timing with respect to the initial time t_ini            !
  !--------------------------------------------------------------------!

  subroutine tng_timing(msg, silent)

    implicit none

    character(len=*), optional, intent(in) :: msg
    logical,          optional, intent(in) :: silent
    integer                                :: cnt, cnt_rate
    logical :: be_silent

    be_silent = .false.
    if (present(silent)) then
       if (silent) then
          be_silent = .true.
       else
          be_silent = .false.
       end if
    end if

    call system_clock(count=cnt, count_rate=cnt_rate)
    tng_t_now = dble(cnt)/dble(cnt_rate)

    if (.not. tng_isinit) then
       tng_t_ini     = tng_t_now
       tng_t_prev    = tng_t_now
       tng_t_lastlog = tng_t_now
       tng_isinit    = .true.
       tng_t         = 0.0d0
       return
    endif

    tng_t         = tng_t + tng_t_now - tng_t_prev
    tng_t_prev    = tng_t_now
    tng_t_lastlog = tng_t_now

    if (.not. be_silent) then
       if (present(msg)) then
          write(u_tng, '(1x,F10.2," s",5x,A)') tng_t, msg
       else
          write(u_tng, '(1x,F10.2," s")') tng_t
       end if
    end if

  end subroutine tng_timing

  !--------------------------------------------------------------------!
  !  Measure just the time passed since the last call to tng_timing()  !
  !--------------------------------------------------------------------!

  subroutine tng_timing2(msg)

    implicit none

    character(len=*), optional, intent(in) :: msg

    integer                                :: cnt, cnt_rate

    call system_clock(count=cnt, count_rate=cnt_rate)
    tng_t_now = dble(cnt)/dble(cnt_rate)

    if (.not. (tng_isinit)) return

    if (present(msg)) then
       write(u_tng, '(5x,F10.2," s",1x,A)') tng_t_now - tng_t_prev, msg
    else
       write(u_tng, '(5x,F10.2," s")') tng_t_now - tng_t_prev
    end if

  end subroutine tng_timing2

  !--------------------------------------------------------------------!
  ! log time elapsed since last call to tng_timing() or tng_timing3()  !
  !--------------------------------------------------------------------!

  subroutine tng_timing3(register)

    implicit none

    integer, intent(in) :: register

    integer          :: cnt, cnt_rate

    if (register > tng_nregisters) return
    if (.not. tng_isinit) return

    call system_clock(count=cnt, count_rate=cnt_rate)
    tng_t_now = dble(cnt)/dble(cnt_rate)

    tng_register(register) = tng_register(register) &
                           + (tng_t_now - tng_t_lastlog)

    tng_t_lastlog = tng_t_now

  end subroutine tng_timing3

  !--------------------------------------------------------------------!
  !        dump current value of timing register to output file        !
  !--------------------------------------------------------------------!

  subroutine tng_dump(register, msg)

    implicit none

    integer, intent(in) :: register
    character(len=*), optional, intent(in) :: msg

    if (register > tng_nregisters) return
    if (.not. tng_isinit) return

    if (present(msg)) then
       write(u_tng, '(1x,F10.2," s",5x,A)') tng_register(register), msg
    else
       write(u_tng, '(1x,F10.2," s",5x,"accumulated in register",1x,I3)') &
            tng_register(register), register
    end if

  end subroutine tng_dump

end module timing
