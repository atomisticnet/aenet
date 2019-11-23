!-----------------------------------------------------------------------
!          arglib.f90 - Easy access to command line arguments
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
! 2011-02-18 Alexander Urban (AU), Nongnuch Artrith (NA)
! 2011-03-19 -- new interface: args_env()
!-----------------------------------------------------------------------

module arglib

  implicit none

  public  :: args_init,     &
             args_final,    &
             args_switch,   &
             args_env,      &
             args_get,      &
             in_list

  private :: args_switch_c1, &
             args_switch_cn, &
             args_switch_l1, &
             args_switch_i1, &
             args_switch_in, &
             args_switch_d1, &
             args_switch_dn, &
             args_env_c,     &
             args_env_i,     &
             args_env_d,     &
             args_env_l


  integer, public :: nargs

  integer, parameter, private :: lenarg = 100

  !--------------------------------------------------------------------!
  !                      args_switch() interface                       !
  !                                                                    !
  ! Get the value(s) of a single command line switch.                  !
  !                                                                    !
  ! Examples:                                                          !
  !                                                                    !
  !    character(len=50) :: input                                      !
  !    call args_switch('--input', value=infile, default='INP')        !
  !    --> infile == 'INP' if the '--input' switch was not found       !
  !        otherwise the value from the command line will be returned  !
  !                                                                    !
  !    logical :: ishelp                                               !
  !    call args_switch('--help:-h', value=ishelp)                     !
  !    --> ishelp is .true., if either the '--help' or the '-h' switch !
  !        was found                                                   !
  !                                                                    !
  !    integer, dimension(3) :: vector                                 !
  !    call args_switch('-n', value=vector, n=3, default=(/ 0, 0, 0 /))!
  !                                                                    !
  ! Available implementations:                                         !
  !                                                                    !
  !    args_switch_c1 : single character string                        !
  !    args_switch_cn : n character strings                            !
  !    args_switch_l1 : logical (see above)                            !
  !    args_switch_i1 : single integer                                 !
  !    args_switch_in : n integers                                     !
  !    args_switch_d1 : singlle doublle precisions                     !
  !    args_switch_dn : n double precisions                            !
  !                                                                    !
  !--------------------------------------------------------------------!

  interface args_switch
     module procedure args_switch_c1, args_switch_cn, args_switch_l1, &
                      args_switch_i1, args_switch_in, args_switch_d1, &
                      args_switch_dn
  end interface

  !--------------------------------------------------------------------!
  !                        args_env() interface                        !
  !                                                                    !
  ! Get the value of an environment variable.                          !
  !                                                                    !
  ! Usage:                                                             !
  !                                                                    !
  !   call args_env(name=n, value=v, stat=s)                           !
  !   n : name of the environment variable (character string)          !
  !   v : value of the variable, if present (character, integer,       !
  !       double precision, or logical)                                !
  !   s : is the return status (optional);                             !
  !       the value of s is 0 if the env. variable was present         !
  !                                                                    !
  ! Available implementations:                                         !
  !                                                                    !
  !    args_env_c   : character string                                 !
  !    args_env_i   : integer value                                    !
  !    args_env_d   : double precision value                           !
  !    args_env_l   : logical value (.true. if the variable is set)    !
  !                                                                    !
  !--------------------------------------------------------------------!

  interface args_env
     module procedure args_env_c, args_env_i, args_env_d, args_env_l
  end interface

contains

  subroutine args_init(check, status)

    implicit none

    character(len=*),  optional :: check
    integer,           optional :: status

    character(len=*), parameter :: figures = '0:1:2:3:4:5:6:7:8:9:.'

    integer                     :: istat, iarg
    character(len=100)          :: arg

    nargs = command_argument_count()

    istat = 0

    ! check if all occuring switches are recognized:
    if (present(check)) then
       chk : do iarg = 1, nargs
          call get_command_argument(iarg, value=arg)
          if ((arg(1:1)=='-') .and. (.not. in_list(figures,arg(2:2)))) then
             if (.not. in_list(check, arg)) then
                write(0,*) 'Warning: Unrecognized command line switch: ', trim(arg)
                istat = iarg
                exit chk
             end if
          end if
       end do chk
    end if

    if (present(status)) then
       status = istat
    end if

  end subroutine args_init

  !--------------------------------------------------------------------!

  subroutine args_final()

    implicit none

    continue
    return

  end subroutine args_final

  !--------------------------------------------------------------------!
  !          alias for the `get_command_aregument' subroutine          !
  !--------------------------------------------------------------------!

  subroutine args_get(iarg, value)

    implicit none

    integer,          intent(in)  :: iarg
    character(len=*), intent(out) :: value

    integer :: i

    if (iarg < 0) then
       i = nargs - iarg + 1
    else
       i = iargs
    end if

    if (i <= nargs) then
       call get_command_argument(iarg, value=value)
    end if

  end subroutine args_get



  !============================= PRIVATE ==============================!



  !--------------------------------------------------------------------!
  !                   args_switch_? implementations                    !
  !--------------------------------------------------------------------!

  subroutine args_switch_c1(switch, value, pos, default)

    implicit none

    character(len=*),           intent(in)    :: switch
    character(len=*), optional, intent(inout) :: value
    integer,          optional, intent(out)   :: pos
    character(len=*), optional, intent(in)    :: default

    integer               :: iarg
    character(len=lenarg) :: arg

    if (present(value) .and. present(default)) value = trim(default)
    if (present(pos)) pos = 0

    do iarg = 1, nargs
       call get_command_argument(iarg, value=arg)
       if (in_list(switch, arg)) then
          if (present(pos)) pos = iarg
          if (present(value)) then
             call get_command_argument(iarg+1, value=value)
          end if
       end if
    end do

  end subroutine args_switch_c1

  !--------------------------------------------------------------------!

  subroutine args_switch_cn(switch, value, n, default)

    implicit none

    integer,                                  intent(in)    :: n
    character(len=*),                         intent(in)    :: switch
    character(len=*), dimension(n),           intent(inout) :: value
    character(len=*), dimension(n), optional, intent(in)    :: default

    integer               :: iarg, i
    character(len=lenarg) :: arg

    if (present(default)) value(1:n) = default(1:n)

    do iarg = 1, nargs
       call get_command_argument(iarg, value=arg)
       if (in_list(switch, arg)) then
          do i = 1, n
             call get_command_argument(iarg+i, value=value(i))
          end do
       end if
    end do

  end subroutine args_switch_cn

  !--------------------------------------------------------------------!

  subroutine args_switch_l1(switch, value)

    implicit none

    character(len=*),           intent(in)    :: switch
    logical,                    intent(inout) :: value

    integer                                   :: ipos

    value = .false.

    call args_switch_c1(switch, pos=ipos)
    if (ipos > 0) then
       value = .true.
    end if

  end subroutine args_switch_l1

  !--------------------------------------------------------------------!

  subroutine args_switch_i1(switch, value, default)

    implicit none

    character(len=*),           intent(in)    :: switch
    integer,                    intent(inout) :: value
    integer,          optional, intent(in)    :: default

    integer                                   :: ipos
    character(len=lenarg)                     :: arg

    if (present(default)) value = default
    call args_switch_c1(switch, value=arg, pos=ipos)
    if (ipos > 0) then
       read(arg, *) value
    end if

  end subroutine args_switch_i1

  !--------------------------------------------------------------------!

  subroutine args_switch_in(switch, value, n, default)

    implicit none

    integer,                         intent(in)    :: n
    character(len=*),                intent(in)    :: switch
    integer, dimension(n),           intent(inout) :: value
    integer, dimension(n), optional, intent(in)    :: default

    integer                                        :: ipos, i
    character(len=lenarg), dimension(n)            :: arg

    if (present(default)) value(1:n) = default(1:n)
    call args_switch_c1(switch, pos=ipos)
    if (ipos > 0) then
       call args_switch_cn(switch, value=arg, n=n)
       do i = 1, n
          read(arg(i), *) value(i)
       end do
    end if

  end subroutine args_switch_in

  !--------------------------------------------------------------------!

  subroutine args_switch_d1(switch, value, default)

    implicit none

    character(len=*),           intent(in)    :: switch
    double precision,           intent(inout) :: value
    double precision, optional, intent(in)    :: default

    integer                                   :: ipos
    character(len=lenarg)                     :: arg

    if (present(default)) value = default
    call args_switch_c1(switch, value=arg, pos=ipos)
    if (ipos > 0) then
       read(arg, *) value
    end if

  end subroutine args_switch_d1

  !--------------------------------------------------------------------!

  subroutine args_switch_dn(switch, value, n, default)

    implicit none

    integer,                                  intent(in)    :: n
    character(len=*),                         intent(in)    :: switch
    double precision, dimension(n),           intent(inout) :: value
    double precision, dimension(n), optional, intent(in)    :: default

    integer                                                 :: ipos, i
    character(len=lenarg), dimension(n)                     :: arg

    if (present(default)) value(1:n) = default(1:n)
    call args_switch_c1(switch, pos=ipos)
    if (ipos > 0) then
       call args_switch_cn(switch, value=arg, n=n)
       do i = 1, n
          read(arg(i), *) value(i)
       end do
    end if

  end subroutine args_switch_dn

  !--------------------------------------------------------------------!
  !            implementations of the args_env_? interface             !
  !--------------------------------------------------------------------!

  subroutine args_env_c(name, value, stat)

    implicit none

    character(len=*),           intent(in)  :: name
    character(len=*),           intent(out) :: value
    integer,          optional, intent(out) :: stat

    integer :: status

    status = 0
    call get_environment_variable(name=trim(name), value=value, status=status)
    if (present(stat)) stat = status

  end subroutine args_env_c

  !--------------------------------------------------------------------!

  subroutine args_env_i(name, value, stat)

    implicit none

    character(len=*),           intent(in)  :: name
    integer,                    intent(out) :: value
    integer,          optional, intent(out) :: stat

    character(len=100) :: str
    integer            :: status

    call args_env_c(name=name, value=str, stat=status)
    if (status == 0) then
       read(str, *) value
    end if
    if (present(stat)) stat = status

  end subroutine args_env_i

  !--------------------------------------------------------------------!

  subroutine args_env_d(name, value, stat)

    implicit none

    character(len=*),           intent(in)  :: name
    double precision,           intent(out) :: value
    integer,          optional, intent(out) :: stat

    character(len=100) :: str
    integer            :: status

    call args_env_c(name=name, value=str, stat=status)
    if (status == 0) then
       read(str, *) value
    end if
    if (present(stat)) stat = status

  end subroutine args_env_d

  !--------------------------------------------------------------------!

  subroutine args_env_l(name, value)

    implicit none

    character(len=*),           intent(in)  :: name
    logical,                    intent(out) :: value

    character(len=100) :: str
    integer            :: status

    call args_env_c(name=name, value=str, stat=status)
    if (status == 0) then
       value = .true.
    else
       value = .false.
    end if

  end subroutine args_env_l

  !--------------------------------------------------------------------!
  ! Search input string `string' for a shorter string `search'.        !
  ! The input string contains records separated by the character       !
  ! `delim' (default = ':').                                           !
  !--------------------------------------------------------------------!

  function in_list(string, search, delim) result(found)

    implicit none

    character(len=*),    intent(in) :: string, search
    character, optional, intent(in) :: delim
    logical                         :: found

    character :: c
    integer   :: i1, i2
    integer   :: slen

    if (present(delim)) then
       c = delim
    else
       c = ':'
    end if

    slen = len_trim(string)

    found = .false.
    i1    = 1
    i2    = scan(string, c)

    do while(i2 >= i1)
       if (trim(adjustl(string(i1:i2-1))) == trim(adjustl(search))) then
          found = .true.
          exit
       end if
       i1 = i2 + 1
       i2 = i1 + scan(string(i1:slen), c) - 1
    end do

    if (trim(adjustl(string(i1:slen))) == trim(adjustl(search))) then
       found = .true.
    end if

  end function in_list

end module arglib
