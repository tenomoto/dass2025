!gfortran ekf.f90 -o ekf -framework Accelerate # mac
!gfortran ekf.f90 -o ekf -LC:\rtools45\x86_64-w64-mingw32.static.posix\lib -llapack -lblas
program ekf
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: &
    seed = 514, nmax = 501, iobs1 = 10, dobs = 10, un = 51
  real(dp), parameter :: &
    dt = 0.001_dp, sr = 0.05_dp, sb = 0.1_dp, sq = 0.03_dp
  real(dp), dimension(6), parameter :: &
    at = [4.0_dp, -2.0_dp, -4.0_dp, -6.0_dp, 2.0_dp, 4.0_dp]

  integer :: ntobs, i, n, t, nf, m
  real(dp), dimension(6) :: a = at
  real(dp), dimension(8) :: xa
  real(dp), dimension(nmax) :: xt, yt
  real(dp), dimension(2, 2) :: &
    ihess, pamat, pfmat, rmat, qmat, hmat, kmat, ikhmat
  real(dp), dimension(8, 8) :: mmat
  integer, dimension(:), allocatable :: tobs, t_hist
  real(dp), dimension(:), allocatable :: x, y, x_hist, y_hist, s
  real(dp), dimension(:, :), allocatable :: p_hist
  real(dp), dimension(:, :), allocatable :: yo

  hmat = diag(2)
  rmat = sr ** 2 * diag(2)
  qmat = sq ** 2 * diag(2)
  pamat = sb ** 2 * diag(2)

  s = set_seed(seed)

  xa(1:2) = [1, 1]
  xa(3:8) = a
  call predict_state(xa, dt, nmax, xt, yt)

  ntobs = (nmax - iobs1) / dobs + 1
  allocate(tobs(ntobs), yo(2, ntobs), t_hist(nmax+ntobs), &
    x_hist(nmax+ntobs), y_hist(nmax+ntobs))
  allocate(p_hist(nmax+ntobs, 3))
  tobs(:) = [(i, i = iobs1, nmax, dobs)]
  yo(1, :) = [(xt(tobs(i)), i = 1, ntobs)] + rnorm(ntobs, 0.0_dp, sr)
  yo(2, :) = [(yt(tobs(i)), i = 1, ntobs)] + rnorm(ntobs, 0.0_dp, sr)

  xa(1:2) = [2, 2]
  xa(3:8) = a(:)

  n = 0
  m = 1
  do t = 1, ntobs
    nf = tobs(t) - n + 1
    t_hist((n + m):(tobs(t) + m)) = [(i, i = n, tobs(t))]
    allocate(x(nf), y(nf))
    p_hist(n+m,:) = [pamat(1,1), pamat(2,2), pamat(1,2)]
    call predict_state(xa, dt, nf, x, y)
    pfmat = pamat
    do i = 1, nf - 1
      call calc_jacobian([x(i), y(i), a(:)], dt, mmat)
      pfmat = matmul(mmat(1:2, 1:2), matmul(pfmat, transpose(mmat(1:2, 1:2)))) + qmat
      p_hist(n+m+i,:) = [pfmat(1,1), pfmat(2,2), pfmat(1,2)]
    end do
    call inv(matmul(hmat, matmul(pfmat, transpose(hmat))) + rmat, ihess)
    kmat = matmul(pfmat, matmul(transpose(hmat), ihess))
    xa(1:2) = [x(nf), y(nf)] + matmul(kmat, (yo(1:2, t) - [x(nf), y(nf)]))
    xa(3:8) = a(:)
    ikhmat = diag(2) - matmul(kmat, hmat)
    pamat = matmul(ikhmat, matmul(pfmat, transpose(ikhmat))) + &
            matmul(kmat, matmul(rmat, transpose(kmat)))
    x_hist((n + m):(tobs(t) + m)) = x(:)
    y_hist((n + m):(tobs(t) + m)) = y(:)
    n = tobs(t)
    m = m + 1
    deallocate(x, y)
    print '(a,i3,2(a,es9.3))', 'Assimilation ',tobs(t), &
    & ': error x = ',abs(xa(1) - xt(tobs(t))),', error y = ',abs(xa(2) - yt(tobs(t)))
  end do
  t_hist(nmax + ntobs) = nmax
  x_hist(nmax + ntobs) = xa(1)
  y_hist(nmax + ntobs) = xa(2)
  p_hist(nmax + ntobs,:) = [pamat(1,1), pamat(2,2), pamat(1,2)]

  print *, sr**2
  open(unit = un, file = "out_ekf.dat", access = "stream", &
    form = "unformatted", status = "replace", action = "write")
  write(unit = un) size(xt), size(t_hist), xt, yt, t_hist, x_hist, y_hist &
  &, p_hist(:,1), p_hist(:,2), p_hist(:,3), sr
  close(unit = un)

  deallocate(tobs, yo, t_hist, x_hist, y_hist, p_hist)

contains

  function diag(n)
    integer, intent(in) :: n
    real(dp), dimension(n, n) :: diag
    
    integer :: i
    
    diag(:, :) = 0.0_dp
    do i = 1, n
      diag(i, i) = 1.0_dp
    end do
    
  end function diag
  
  subroutine inv(amat, iamat)
    real(dp), dimension(:, :), intent(in) :: amat
    real(dp), dimension(:, :), intent(inout) :: iamat
    
    character(len=1), parameter :: uplo = "u"
    integer, parameter :: nb = 64 ! system dependent
    integer :: n, lda, lwork, info
    integer, dimension(:), allocatable :: ipiv
    real(dp), dimension(:), allocatable :: work
    
    n = size(amat, 1)
    lda = n
    lwork = n * nb
    iamat(:, :) = amat(:, :)
    allocate(ipiv(n), work(lwork))
    call dsytrf(uplo, n, iamat, lda, ipiv, work, lwork, info)
    call dsytri(uplo, n, iamat, lda, ipiv, work, info)
    deallocate(ipiv, work)

  end subroutine inv
  
  subroutine predict_state(w, dt, nmax, x, y)
    real(dp), dimension(:), intent(in) :: w
    real(dp), intent(in) :: dt
    integer, intent(in) :: nmax
    real(dp), dimension(nmax), intent(inout) :: x, y

    integer :: n
    real(dp), dimension(6) :: a

    x(1) = w(1)
    y(1) = w(2)
    a(:) = w(3:8)
    do n = 1, nmax-1
      x(n + 1) = x(n) + dt * (x(n) * (a(1) + a(2) * x(n) + a(3) * y(n)))
      y(n + 1) = y(n) + dt * (y(n) * (a(4) + a(5) * y(n) + a(6) * x(n)))
    end do

  end subroutine predict_state

  subroutine calc_jacobian(w, dt, mmat)
    real(dp), dimension(:), intent(in) :: w
    real(dp), intent(in) :: dt
    real(dp), dimension(:, :), intent(inout) ::  mmat

    integer :: i, n
    real(dp)  :: x, y
    real(dp), dimension(6) :: a
    real(dp), dimension(8, 8) :: jmat

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
    jmat(3:8, 3:8) = diag(6)
    mmat = diag(8) + dt * jmat

  end subroutine calc_jacobian

  function set_seed(seed) result(s)
    integer :: seed
    integer :: n
    integer, allocatable :: s(:)

    call random_seed(size=n)
!    print *, "seed size=", n
    allocate(s(n))
    s(:) = seed
    call random_seed(put=s)
    call random_seed(get=s)
!    print *, "seed=", s

  end function set_seed

  function runif(n) result(u)
    integer, intent(in) :: n
    real(dp), allocatable :: u(:)

    allocate(u(n))
    call random_number(u)
    u = 1.0_dp - u

  end function runif

  function rnorm(n, mean, sd) result(x)
    integer, intent(in) :: n 
    real(dp), intent(in), optional :: mean, sd
    real(dp), allocatable :: x(:)

    real(dp), allocatable :: u1(:), u2(:)
    real(dp) :: tau

    tau = 2.0_dp * acos(-1.0_dp)
    allocate(u1(n), u2(n), x(n))
    u1 = runif(n)
    u2 = runif(n)
    x = sqrt(-2.0_dp * log(u1)) * cos(tau * u2)
    if (present(sd)) x = sd * x
    if (present(mean)) x = mean + x
    deallocate(u1, u2)

  end function rnorm

end program ekf

