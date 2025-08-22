!gfortran predprey.f90 -o predprey -framework Accelerate # mac
!gfortran predprey.f90 -o predprey -LC:\rtools45\x86_64-w64-mingw32.static.posix\lib
program predprey
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: nmax = 15001, un = 51
  real(dp), parameter :: dt = 0.001d0

  real(dp), dimension(6) :: a
  real(dp), dimension(nmax) :: x, y
  integer :: n

  x(1) = 1
  y(1) = 1
  a = (/4.0_dp, -2.0_dp, -4.0_dp, -6.0_dp, 2.0_dp, 4.0_dp/)

  do n = 1, nmax-1
    x(n+1) = x(n) + dt * (x(n) * (a(1) + a(2) * x(n) + a(3) * y(n)))
    y(n+1) = y(n) + dt * (y(n) * (a(4) + a(5) * y(n) + a(6) * x(n)))
  end do

  open(unit = un, file = "xy.dat", access = "stream", &
    form = "unformatted",  action = "write", status = "replace")
  write(unit = un) nmax
  write(unit = un) dt
  write(unit = un) x
  write(unit = un) y
  close(un)

end program predprey