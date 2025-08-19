! gfortran adjoint_params.f90 -o adjoint_params
program adjoint_params_main
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none

  logical, parameter :: debug = .true.
  integer, parameter :: &
    nmax = 501, iobs1 = 2, dobs = 2, maxit = 100, un = 51
  real(dp), parameter :: dt = 0.001_dp, &
    ctol = 1.0d-10, gtol = 1.0d-5, stol = 1.0e-7
  real(dp), dimension(6), parameter :: &
    at = [4.0_dp, -2.0_dp, -4.0_dp, -6.0_dp, 2.0_dp, 4.0_dp]
  real(dp), parameter :: x1 = 1.0_dp, y1 = 1.0_dp
  
  integer :: i, n, ntobs, res, iter
  integer, dimension(:), allocatable :: tobs
  real(dp), dimension(nmax) :: x, y, xt, yt, xo, yo, ax, ay
  real(dp), dimension(6) :: a, aa
  real(dp), dimension(8) :: par
  real(dp), dimension(:), allocatable :: hcost, hgnorm
  real(dp), dimension(:, :), allocatable :: hpar

  ntobs = (nmax - iobs1) / dobs + 1
  allocate(tobs(ntobs), hcost(maxit+1), hgnorm(maxit+1), hpar(8, maxit+1))
  tobs(:) = [(i, i = iobs1, nmax, dobs)]

  a = at
  x(1) = x1
  y(1) = y1
  call forward
  xt(:) = x
  yt(:) = y
  xo(:) = xt(:)
  yo(:) = yt(:)
  
  a = [1.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp]
  par(1:6) = a(:)
  par(7:8) = [2.0_dp, 2.0_dp]

  res = bfgs(par)
  print *, "Convergence: ", res, " Iteration: ", iter
  print *, "x,y = ", par(7:8)
  print *, "xparms = ", par(1:3)
  print *, "yparms = ", par(4:6)

  open(unit = un, file = "out_adjoint_params.dat", access = "stream", &
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
    real(dp), dimension(8) :: grad

    x(1) = par(7)
    y(1) = par(8)
    a(:) = par(1:6)
    call forward
    call adjoint
    grad = [aa(1), aa(2), aa(3), aa(4), aa(5), aa(6), ax(1), ay(1)]

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
    
  function bfgs(par) result(convergence)
    real(dp), dimension(8), intent(inout) :: par
 
    integer :: n, convergence
    real(dp) :: alpha, ys, rho, gnorm
    real(dp), dimension(size(par)) :: g, new_par, s, y, new_gr
    real(dp), dimension(:, :), allocatable :: H

    n = size(par)
    allocate(H(n, n))
    H = diag(n)
    g = gr(par)
    alpha = min(1.0_dp, 1.0_dp / sum(abs(g)))
    new_par = par - alpha * g
    s = -alpha * g
    y = gr(new_par) - g
    ys = sum(y * s)
    H = ys / sum(y * y) * diag(n)
    if (debug) print *, "Initial Hessian:", H
  
    do iter = 1, maxit
      hcost(iter) = fn(par)
      gnorm = sqrt(sum(g**2))
      hgnorm(iter) = gnorm
      hpar(:, iter) = par(:)
      if (debug) print *, &
        "Iteration:", iter, "Cost:", hcost(iter), "Gradient Norm:", gnorm
      if (debug) print *, par(1:3)
      if (debug) print *, par(4:6)
      if (debug) print *, par(7:8)
      if (gnorm < gtol) then
        convergence = 0
        return
      end if
      if (maxval(abs(s)) < stol) then
        convergence = 1
        return
      end if

      new_par = par - matmul(H, g)
      s = new_par - par
      new_gr = gr(new_par)
      y = new_gr - g
      ys = sum(y * s)
      if (ys > ctol) then
        rho = 1 / ys
        H = matmul((diag(n) - rho * outer(s, y)), matmul(H, (diag(n) - rho * outer(y, s)))) &
            + rho * outer(s, s)
!        if (debug) print *, "Updated Hessian:", H
      else
        H = ys / sum(y * y) * diag(n)
!        if (debug) print *, "Reset Hessian:", H
      end if
      if (iter == maxit) convergence = 2
    
      par(:) = new_par(:)
      g = new_gr
    end do

    deallocate(H)
  
  end function bfgs

end program adjoint_params_main
