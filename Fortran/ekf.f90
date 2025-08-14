program ekf
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: &
    seed = 514, nmax = 501, iobs1 = 10, dobs = 10
  real(dp), parameter :: dt = 0.01, sr = 0.01_dp, sb = 0.1_dp, sq = 0.1_dp

  real(dp), dimension(:), parameter :: at = [4, -2, -4, -6, 2, 4]
  real(dp), dimension(:), parameter :: a = [1, 0, 0, -1, 0, 0]
  real(dp), parameter, dimension(2, 2) :: rmat = [[sr**2, 0], [0, sr**2]]
  real(dp), parameter, dimension(2, 2) :: hmat = [[1, 0], [0, 1]]
  real(dp), parameter, dimension(2, 2) :: qmat = [[sq**2, 0], [0, sq**2]]

  integer :: ntobs, i, n, t, nf
  integer, dimension(nmax + ntobs) :: t_hist
  real(dp), dimension(8) :: w1, xa
  real(dp), dimension(nmax) :: x, y
  real(dp), dimension(:), allocatable :: tobs
  real(dp), dimension(2, 2) :: pamat = [[sb**2, 0], [0, sb**2]]
  real(dp), dimension(2, 2) :: mmat = [[1, 0], [0, 1]], mmat1
  real(dp), dimension(nmax + ntobs) :: x_hist, y_hist

!set.seed(seed)
  w1(1:2) = [1, 1]
  w1(3:8) = at
  call predict_state(w1, dt, nmax, x, y)

  ntobs = (nmax - iobs1 + 1) / dobs + 1
  allocate(tobs(ntobs))
  tobs = [i, i = iobs1, nmax, dobs]
  yo = wt(, tobs) !+ matrix(rnorm(2 * ntobs, 0, sr), 2, ntobs)

  xa(1:2) = [2, 2]
  xa(3:8) = a(:)

  n = 0
  do t = 1, ntobs
    nf = tobs(t) - n + 1
    t_hist = c(t_hist, n:tobs(t), tobs(t))
    call predict_state(xa, dt, nf, x, y)
    mmat = diag(2)
    do i = 1, nf - 1
      call calc_jacobian(c(xf(, i), a), dt, mmat1)
      mmat = matmul(mmat1(1:2, 1:2), mmat)
    end do
    pfmat = matmul(mmat, matmul(pamat, transpose(mmat))) + qmat
    kmat = matmul(matmul(pfmat, transpose(hmat)), solve(hmat %*% pfmat %*% t(hmat) + rmat))
    xa(1:2) = w + matmul(kmat, (yo - w))
    xa(3:8) = a(:)
    n = nrow(pfmat)
    ikhmat = diag(n) - matmul(kmat, hmat)
    pamat = matmul(ikhmat, matmul(pfmat, t(ikhmat))) + matmul(kmat, matmul(rmat, transpose(kmat)))
    x_hist = c(x_hist, xf(1, ), xa(1))
    y_hist = c(y_hist, xf(2, ), xa(2))
    n = tobs(t)
  end do

  deallocate(tobs)

contains

  subroutine predict_state(w, dt, nmax, x, y)
    real(kind=dp), dimension(:), intent(in) :: w
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(nmax), intent(inout) :: x, y

    integer :: n
    real(kind=dp), dimension(6) :: a

    x(1) = w(1)
    y(1) = w(2)
    a(:) = w(3:8)
    do n = 1, nmax-1
      x(n + 1) = x(n) + dt * (x(n) * (a(1) + a(2) * x(n) + a(3) * y(n)))
      y(n + 1) = y(n) + dt * (y(n) * (a(4) + a(5) * y(n) + a(6) * x(n)))
    end do

  end subroutine predict_state

  subrotuine calc_jacobian(w, dt)
    real(kind=dp), dimension(:), intent(in) :: w
    real(kind=dp), intent(in) :: dt
    real(kind=dp), dimension(8, 8) :: jmat, mmat

    integer :: i
    real(kind=dp)  :: x, y
    real(kind=dp), dimension(6) :: a

    x = w(1)
    y = w(2)
    a = w(3:8)

    jmat(:, :) = 0.0_dp
    jmat(1, 1) = a(1) + 2 * a(2) * x + a(3) * y
    jmat(1, 2) = a(3) * x
    jmat(1, 3:5) = [x, x**2, x * y]
    jmat(2, 1) = a(6) * y
    jmat(2, 2) = a(4) + 2 * a(5) * y + a(6) * x
    jmat(2, 6:8) = [y, y**2, x * y]
    do i = 3, 8
      jmat(i, i) = 1.0_dp
    end do
    mmat = dt * jmat
    do i = 1, 8
      mmat(i, i) = 1.0_dp + mmat(i, i)
    end do

  end subroutine predict_state

end program ekf

