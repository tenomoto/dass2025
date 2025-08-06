forward <- function(w, dt, nmax) {
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

tlm <- function(dw, x, y, dt, a, da = rep(0, length(a))) {
  dx <- dw[1]
  dy <- dw[2]
  nmax <- length(x)
  for (n in 1:(nmax-1)) {
    dx <- dx + dt * (a[1] + 2 * a[2] * x[n] + a[3] * y[n]) * dx + a[3] * x[n] * dy +
      a[1] * x[n] * da[1] + a[2] * x[n]^2 * da[2] + a[3] * x[n] * y[n] * da[3]
    dy <- dy + dt * a[6] * y[n] * dx +  (a[4] + 2 * a[5] * y[n] + a[6] * x[n]) * dy +
      a[4] * y[n] * da[4] + a[5] * y[n]^2 * da[5] + a[6] * x[n] * y[n] * da[6]
  }
  c(dx, dy)
}

forecast_state <- function(xa, mfunc, ...) {
  mfunc(xa, ...)
}

forecast_ecov <- function(pamat, tfunc, sq, ...) {
  n <- nrow(pamat)
  mpmat <- matrix(0, n, n)
  pfmat <- matrix(0, n, n)
  for (i in 1:n) {
    mpmat[, i] <- tfunc(pamat[, i], ...)
  }
  for (i in 1:n) {
    pfmat[, i] <- tfunc(mpmat[i, ], ...)
  }
  pfmat + matrix(rnorm(n * n, 0, sq), n, n)
}

calc_kgain <- function(pfmat, hmat, rmat) {
  pfhtmat <- pfmat %*% t(hmat)
  pfhtmat %*% solve(hmat %*% pfhtmat + rmat)
}

analyze_state <- function(xf, kmat, yo, hfunc) {
  xf + kmat %*% (yo - hfunc(xf))
}

analyze_ecov <- function(pfmat, kmat, hmat, rmat) {
  n <- nrow(pfmat)
  ikhmat <- -kmat %*% hmat
  diag(ikhmat) <- diag(ikhmat) + 1
  ikhmat %*% pfmat %*% t(ikhmat) + kmat %*% rmat %*% t(kmat)
}

nmax <- 501
dt <- 0.001
at <- c(4, -2, -4, -6, 2, 4)
x1 <- 1
y1 <- 1
tobs <- seq(2, nmax, by = 2)
ntobs <- length(tobs)
forward.result <- forward(dt, at, x1, y1, nmax)
yo <- rbind(forward.result$x, forward.result$y)

xf <- c(2, 2)
#a <- c(1, 0, 0, -1, 0, 0)
a <- at
x_hist <- numeric(0)
y_hist <- numeric(0)
t_hist <- numeric(0)

hfunc <- function(x) x
hmat <- diag(1, 2, 2)
sb <- 0.1
sr <- 0.1
pamat <- diag(sb, 2, 2)
rmat <- diag(sr, 2, 2)
sq <- 0.1

n <- 1
for (t in 1:ntobs) {
  nmax <- tobs[t] - n + 1
  n <- nmax
  t_hist <- c(t_hist, n:tobs[t], tobs[t])
  xa <- c(xf[1], xf[2], a)
  xf <- forecast_state(xa, forward, dt = dt, nmax = nmax)
  pfmat <- forecast_ecov(pamat, tlm, sq, x = xf$x, y = xf$y, dt = dt, a = a)
  kmat <- calc_kgain(pfmat, hmat, rmat)
  xa <- analyze_state(xf, kmat, yo[,t], hfunc)
  pamat <- analyze_ecov(pfmat, kmat, hmat, rmat)
  xa_hist <- c(xa_hist, xf$x, xa[1])
  ya_hist <- c(ya_hist, yf$y, xa[2])
}