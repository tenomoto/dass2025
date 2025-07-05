subroutine adjoint(dt, a, x, y, xo, yo, tobs, aa, ax_out, ay_out)
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  real(dp), intent(in) :: dt
  real(dp), dimension(6), intent(in) :: a
  real(dp), dimension(:), intent(in) :: x, y
  real(dp), dimension(:), intent(in) :: xo, yo
  integer, intent(in) :: tobs

  real(dp), dimension(6), intent(out) :: aa
  real(dp), intent(out) :: ax_out, 
  real(dp), intent(out) :: ay_out

  integer :: nmax, n
  real(dp), dimension(size(x)) :: ax, ay

  nmax = size(x)
  aa(:) = 0.0_dp
  ax(:) = 0.0_dp
  ay(:) = 0.0_dp
  do nmax-1, 1, -1
    aa(6) = aa(6) + dt * x(n) * y(n) * ay(n+1)
    aa(5) = aa(5) + dt * y(n) * y(n) * ay(n+1)
    aa(4) = aa(4) + dt * y(n) * ay(n+1)
    ax(n) = ax(n) + dt * a(6) * y(n) * ay(n+1)
    ay(n) = ay(n) + dt * a(5) * y(n) * ay(n+1)
    ay(n) = ay(n) + (1 + dt * (a(4) + a(5) * y(n) + a(6) * x(n))) * ay(n+1)
    aa(3) = aa(3) + dt * y(n) * x(n) * ax(n+1)
    aa(2) = aa(2) + dt * x(n) * x(n) * ax(n+1)
    aa(1) = aa(1) + dt * x(n) * ax(n+1)
    ay(n) = ay(n) + dt * a(3) * x(n) * ax(n+1)
    ax(n) = ax(n) + dt * a(2) * x(n) * ax(n+1)
    ax(n) = ax(n) + (1 + dt * (a(1) + a(2) * x(n) + a(3) * y(n))) * ax(n+1)
    if (any(n == tobs)) then 
      ax(n) = ax(n) + (x(n) - xo(n))
      ay(n) = ay(n) + (y(n) - yo(n))
    end if
  end do

end subroutine adjoint
