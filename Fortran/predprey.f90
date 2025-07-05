program predprey
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: nmax = 15001, un = 51
  real(dp), parameter :: dt = 0.01d0

  real(dp), dimension(6) :: a
  real(dp), dimension(nmax) :: x, y
  integer :: n

  x(1) = 1
  y(1) = 1

  do n = 1, nmax-1
    x(n+1) = x(n) + dt * (x(n) * (a(1) + a(2) * x(n) + a(3) * y(n)))
    y(n+1) = y(n) + dt * (y(n) * (a(4) + a(5) * y(n) + a(6) * x(n)))
  end do

  open(unit = un, file = "xy.dat", access = "stream", &
    form = "unformatted",  action = "write", status = "replace")
  write(unit = un) x
  write(unit = un) y
  close(un)

end program predprey