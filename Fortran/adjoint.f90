! gfortran adjoint.f90 -o adjoint
program adjoint_main
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  logical, parameter :: debug = .true.
  integer, parameter :: &
    nmax = 501, iobs1 = 2, dobs = 2, maxit = 500, logit=20, un = 51
  real(dp), parameter :: dt = 0.001_dp, alpha = 1.0d-3, gtol = 1.0d-5
  real(dp), dimension(6), parameter :: &
    at = [4.0_dp, -2.0_dp, -4.0_dp, -6.0_dp, 2.0_dp, 4.0_dp]
  real(dp), parameter :: x1 = 1.0_dp, y1 = 1.0_dp
  
  integer :: i, n, ntobs, res, iter
  integer, dimension(:), allocatable :: tobs
  real(dp), dimension(nmax) :: x, y, xt, yt, xo, yo, ax, ay
  real(dp), dimension(6) :: a, aa
  real(dp), dimension(8) :: par
  real(dp) :: gnorm
  real(dp), dimension(2) :: g, d
  real(dp), dimension(:), allocatable :: hcost, hgnorm
  real(dp), dimension(:, :), allocatable :: hpar

  ntobs = (nmax - iobs1) / dobs + 1
  allocate(tobs(ntobs), hcost(maxit+1), hgnorm(maxit+1), hpar(2, maxit+1))
  tobs(:) = [(i, i = iobs1, nmax, dobs)]

  a = at
  x(1) = x1
  y(1) = y1
  call forward
  xt(:) = x
  yt(:) = y
  xo(:) = xt(:)
  yo(:) = yt(:)
  
  ! initialize control variables
  par(1:6) = a(:)
  par(7:8) = [1.5_dp, 1.5_dp]
  print *, par

  print *, "start optimization"
  res = 0
  do iter = 1, maxit
    ! evaluate cost function and its gradient
    hcost(iter) = fn(par)
    g = gr(par)
    gnorm = sqrt(sum(g**2))
    hgnorm(iter) = gnorm
    hpar(:, iter) = par(7:8)
    if (gnorm < gtol) then
      print *, "Convergence, Iteration: ", iter
      print *, "x,y = ", par(7:8)
      res = 1
      exit
    end if
    if (mod(iter,logit).eq.0) &
    &  print '(a,i3,2(a,es9.3))', "Iteration:", iter, " Cost:", hcost(iter), " Gradient Norm:", gnorm
    ! update control variables with steepest descent
    d = -alpha * g 
    par(7:8) = par(7:8) + d
  end do
  if (res.ne.1) then
    print *, "reaching maximum iterations"
    print *, "x,y = ", par(7:8)
  end if

  open(unit = un, file = "out_adjoint.dat", access = "stream", &
    form = "unformatted", status = "replace", action = "write")
  write(unit = un) iter, hcost(1:iter), hgnorm(1:iter), hpar(:, 1:iter)
  close(unit = un)

  deallocate(tobs, hcost, hgnorm, hpar)

contains

  subroutine forward

    do n = 1, nmax-1
      x(n+1) = x(n) + dt * (x(n) * (a(1) + a(2) * x(n) + a(3) * y(n)))
      y(n+1) = y(n) + dt * (y(n) * (a(4) + a(5) * y(n) + a(6) * x(n)))
    end do

  end subroutine forward

  subroutine adjoint

    aa(:) = 0.0_dp
    ax(:) = 0.0_dp
    ay(:) = 0.0_dp
    do n = nmax-1, 1, -1
      if (any(n == tobs)) then 
        ax(n) = ax(n) + (x(n) - xo(n))
        ay(n) = ay(n) + (y(n) - yo(n))
      end if
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
    end do

  end subroutine adjoint

  function calc_cost() result(cost)
    real(dp) :: cost

    cost = 0.5_dp * sum((x(tobs) - xo(tobs))**2 + (y(tobs) - yo(tobs))**2)

  end function calc_cost

  function fn(par) result(cost)
    real(dp), dimension(8), intent(in) :: par
    real(dp) :: cost

    x(1) = par(7)
    y(1) = par(8)
    a(:) = par(1:6)
    call forward
    cost = calc_cost()

  end function fn

  function gr(par) result(grad)
    real(dp), dimension(8), intent(in) :: par
    real(dp), dimension(2) :: grad

    x(1) = par(7)
    y(1) = par(8)
    a(:) = par(1:6)
    call forward
    call adjoint
    grad = [ax(1), ay(1)]

  end function gr

  function diag(n)
    integer, intent(in) :: n
    real(dp), dimension(n, n) :: diag
    
    integer :: i
    
    diag(:, :) = 0.0_dp
    do i = 1, n
      diag(i, i) = 1.0_dp
    end do

  end function diag

  function outer(a, b) result(c)
    real(dp), dimension(:), intent(in) :: a, b
    real(dp), dimension(:, :), allocatable :: c

    integer :: i, j

    allocate(c(size(a), size(b)))
    do j = 1, size(b)
      do i = 1, size(a)
        c(i, j) = a(i) * b(j)
      end do
    end do

  end function outer

end program adjoint_main
