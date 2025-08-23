library(torch)

calc_tendency <- function(a, x, y) {
  dxdt <- x * (a[1] + a[2] * x + a[3] * y)
  dydt <- y * (a[4] + a[5] * y + a[6] * x)
  torch_stack(c(dxdt, dydt), dim = -1)
}

calc.cost <- function(wf, wo) {
  0.5 * (torch_sum((wf - wo)^2))
}

dev <- "cpu"
dtype <- torch_float()
nmax <- 501
dt <- 0.001
loss_old <- Inf
stol <- 1e-8

at <- torch_tensor(c(4, -2, -4, -6, 2, 4),
                   device = dev, dtype = dtype)
wt <- torch_zeros(nmax, 2, device = dev, dtype = dtype)
wt[1, 1:2] <- c(1, 1)
for (n in 1:(nmax-1)) {
  dwdt <- calc_tendency(at, wt[n, 1], wt[n, 2])
  wt[n+1, ] <- wt[n, ] + dt * dwdt
}
xt <- as.numeric(wt[, 1])
yt <- as.numeric(wt[, 2])

tobs <- seq(2, nmax, by = 2)
wo <- wt[tobs]

a <- c(1, 0, 0, -1, 0, 0)
x1 <- 2
y1 <- 2
par <- torch_tensor(c(a, x1, y1), requires_grad = TRUE,
                    device = dev, dtype = dtype)
print(par$dtype)
print(par$device)
hist <- list(cost = numeric(0), gnorm = numeric(0), par = vector(length=0))

optimizer <- optim_lbfgs(par, max_iter = 20, tolerance_change = 1e-9, tolerance_grad = 1e-5)

maxit <- 10
for (i in 1:maxit) {
  closure <- function() {
    hist$par <<- rbind(hist$par, as.numeric(par))
    optimizer$zero_grad()
    w_hist <- list()
    w_curr <- par[7:8]$view(c(1, 2))
    w_hist[[1]] <- w_curr
    for (n in 1:(nmax-1)) {
      dwdt <- calc_tendency(par[1:6], w_curr[1, 1], w_curr[1, 2])
      w_next <- w_curr + dt * dwdt
      w_hist[[n + 1]] <- w_next
      w_curr <- w_next
    }
    w <- torch_cat(w_hist, dim = 1)

    cost <- calc.cost(w[tobs, ], wo)
    cost$backward()
    hist$cost <<- c(hist$cost, as.numeric(cost))
    grad <- par$grad
    gnorm <- as.numeric(grad$dot(grad)$sqrt())
    hist$gnorm <<- c(hist$gnorm, gnorm)
    cost
  }
  
  loss <- optimizer$step(closure)
  
  x <- as.numeric(par[7])
  y <- as.numeric(par[8])
  cat("Iteration:", i, "loss=", as.numeric(loss), "x=", x, "y=", y, "\n")
  if (abs(loss_old - as.numeric(loss)) < stol) {
    cat("Converged at iteration", i, "\n")
    break 
  }
  loss_old <- as.numeric(loss)
}

alg <- "lbfgs"
cat("Optimized Parameters:", as.numeric(par), "\n")
save(alg, hist, file = "out_adjoint_t4R.RData")