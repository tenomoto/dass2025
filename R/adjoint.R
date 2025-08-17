forward <- function(dt, a, x1, y1, nmax) {
  x <- rep(0, nmax)
  y <- rep(0, nmax)
  x[1] <- x1
  y[1] <- y1
  for (n in 1:(nmax-1)) {
    x[n+1] <- x[n] + dt * (x[n] * (a[1] + a[2] * x[n] + a[3] * y[n]))
    y[n+1] <- y[n] + dt * (y[n] * (a[4] + a[5] * y[n] + a[6] * x[n]))
  }
  list(x = x, y = y)
}

adjoint <- function(dt, a, x, y, xo, yo, tobs) {
  nmax <- length(x)
  aa <- rep(0, length(a))
  ax <- rep(0, nmax)
  ay <- rep(0, nmax)
  for (n in (nmax-1):1) {
    aa[6] <- aa[6] + dt * x[n] * y[n] * ay[n+1]
    aa[5] <- aa[5] + dt * y[n] * y[n] * ay[n+1]
    aa[4] <- aa[4] + dt * y[n] * ay[n+1]
    ax[n] <- ax[n] + dt * a[6] * y[n] * ay[n+1]
    ay[n] <- ay[n] + dt * a[5] * y[n] * ay[n+1]
    ay[n] <- ay[n] + (1 + dt * (a[4] + a[5] * y[n] + a[6] * x[n])) * ay[n+1]
    aa[3] <- aa[3] + dt * y[n] * x[n] * ax[n+1]
    aa[2] <- aa[2] + dt * x[n] * x[n] * ax[n+1]
    aa[1] <- aa[1] + dt * x[n] * ax[n+1]
    ay[n] <- ay[n] + dt * a[3] * x[n] * ax[n+1]
    ax[n] <- ax[n] + dt * a[2] * x[n] * ax[n+1]
    ax[n] <- ax[n] + (1 + dt * (a[1] + a[2] * x[n] + a[3] * y[n])) * ax[n+1]
    if (n %in% tobs) {
      ax[n] <- ax[n] + (x[n] - xo[n])
      ay[n] <- ay[n] + (y[n] - yo[n])
    }
  }
  c(aa, ax[1], ay[1])
}

fn <- function(par, dt, nmax, xo, yo, tobs) {
  xyf <- forward(dt, par[1:6], par[7], par[8], nmax)
  cost <- calc.cost(xyf$x[tobs], xyf$y[tobs], xo[tobs], yo[tobs])
  cost
}

gr <- function(par, dt, nmax, xo, yo, tobs){
  xyf <- forward(dt, par[1:6], par[7], par[8], nmax)
  hist$par <<- rbind(hist$par, par)
  cost <- calc.cost(xyf$x[tobs], xyf$y[tobs], xo[tobs], yo[tobs])
  hist$cost <<- c(hist$cost, cost)
  grad <- adjoint(dt, par[1:6], xyf$x, xyf$y, xo, yo, tobs)
  hist$gnorm <<- c(hist$gnorm, sqrt(sum(grad^2)))
  grad
}

calc.cost <- function(xf, yf, xo, yo) {
  0.5 * (sum((xf - xo)^2 + (yf - yo)^2)) 
}

nmax <- 501
dt <- 0.001
at <- c(4, -2, -4, -6, 2, 4)
x1 <- 1
y1 <- 1

tobs <- seq(2, nmax, by = 2)
forward.result <- forward(dt, at, x1, y1, nmax)
xt <- forward.result$x
yt <- forward.result$y
xo <- xt
yo <- yt

a <- c(1, 0, 0, -1, 0, 0)
x1 <- 2
y1 <- 2
par <- c(a, x1, y1)
cntl <- list(maxit = 100, reltol = 1e-5)

hist <- list(cost = numeric(0), gnorm = numeric(0), par = vector(length=0))

alg <- "BFGS"
res <- optim(par, fn, gr, method = alg, control = cntl, dt, nmax, xo, yo, tobs)

save(alg, hist, file = "out_adjoint.RData")
