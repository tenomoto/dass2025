library(optimx)

forward <- function(dt, x1, y1, z1, nmax) {
  x <- rep(0, nmax)
  y <- rep(0, nmax)
  z <- rep(0, nmax)
  x[1] <- x1
  y[1] <- y1
  z[1] <- z1
  for (n in 1:(nmax - 1)) {
    x[n + 1] <- x[n] + dt * s * (-x[n] + y[n])
    y[n + 1] <- y[n] + dt * (-x[n] * z[n] + r * x[n] - y[n])
    z[n + 1] <- z[n] + dt * (x[n] * y[n] - b * z[n])
  }
  list(x = x, y = y, z = z)
}

adjoint <- function(dt, w, wo, tobs) {
  nmax <- length(w$x)
  ax <- rep(0, nmax)
  ay <- rep(0, nmax)
  az <- rep(0, nmax)
  ar <- 0
  as <- 0
  ab <- 0
  x <- w$x
  y <- w$y
  z <- w$z
  for (n in (nmax - 1):1) {
    if (n %in% tobs) {
      ax[n] <- ax[n] + (x[n] - wo$x[n])
      ay[n] <- ay[n] + (y[n] - wo$y[n])
      az[n] <- az[n] + (z[n] - wo$z[n])
    }
    az[n] <- az[n] + (1 - dt * b) * az[n + 1]
    ay[n] <- ay[n] + dt * x[n] * az[n + 1]
    ax[n] <- ax[n] + dt * y[n] * az[n + 1]
    az[n] <- az[n] - dt * x[n] * ay[n + 1]
    ay[n] <- ay[n] + (1 - dt) * ay[n + 1]
    ax[n] <- ax[n] + dt * (-z[n] + r) * ay[n + 1]
    ay[n] <- ay[n] + s * dt * ax[n + 1]
    ax[n] <- ax[n] + (1 - s * dt) * ax[n + 1]
  }
  c(ax[1], ay[1], az[1])
}

fn <- function(par, dt, nmax, wo, tobs) {
  w <- forward(dt, par[1], par[2], par[3], nmax)
  cost <- calc.cost(w, wo, tobs)
  cost
}

gr <- function(par, dt, nmax, wo, tobs){
  w <- forward(dt, par[1], par[2], par[3], nmax)
  hist$par <<- rbind(hist$par, par)
  cost <- calc.cost(w, wo, tobs)
  hist$cost <<- c(hist$cost, cost)
  grad <- adjoint(dt, w, wo, tobs)
  hist$gnorm <<- c(hist$gnorm, sqrt(sum(grad^2)))
  grad
}

calc.cost <- function(wf, wo, tobs) {
  0.5 * (sum((wf$x[tobs] - wo$x[tobs])^2 + 
               (wf$y[tobs] - wo$y[tobs])^2 + (wf$z[tobs] - wo$z[tobs])^2)) 
}

nmax <- 200
r <- 32.0
s <- 10.0
b <- 8 / 3
dt <- 0.01
x1 <- 1
y1 <- 3
z1 <- 5

dobs <- 60
iobs <- dobs

tobs <- seq(iobs, nmax, by = dobs)
wt <- forward(dt, x1, y1, z1, nmax)
wo <- wt

x1 <- 1.1
y1 <- 3.3
z1 <- 5.5
par <- c(x1, y1, z1)
cntl <- list(maxit = 100, reltol = 1e-8)

hist <- list(cost = numeric(0), gnorm = numeric(0), par = vector(length=0))

#alpha <- 5e-4
#for (i in 1:100) {
#  grad <- gr(par, dt, nmax, wo, tobs)
#  cost <- fn(par, dt, nmax, wo, tobs)
#  par <- par - alpha * grad
#}

#alg <- "BFGS"
#res <- optim(par, fn, gr, method = alg, control = cntl, dt, nmax, wo, tobs)

alg <- "nvm"
res <- optimr(par, fn, gr, method = alg,
              dt = dt, nmax = nmax, wo = wo, tobs = tobs)