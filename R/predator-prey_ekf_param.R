forward <- function(w, dt, nmax) {
  x <- rep(0, nmax)
  y <- rep(0, nmax)
  x[1] <- w[1]
  y[1] <- w[2]
  a <- w[3:8]
  for (n in 1:(nmax-1)) {
    x[n+1] <- x[n] + dt * (x[n] * (a[1] + a[2] * x[n] + a[3] * y[n]))
    y[n+1] <- y[n] + dt * (y[n] * (a[4] + a[5] * y[n] + a[6] * x[n]))
  }
  rbind(x, y)
}

tlm <- function(dw, x, y, dt, a, da = rep(0, length(a))) {
  dx <- dw[1]
  dy <- dw[2]
  nmax <- length(x)
  for (n in 1:(nmax-1)) {
    dx <- dx + dt * (
      (a[1] + 2 * a[2] * x[n] + a[3] * y[n]) * dx + a[3] * x[n] * dy +
      x[n] * da[1] + x[n]^2 * da[2] + x[n] * y[n] * da[3])
    dy <- dy + dt * (
      a[6] * y[n] * dx + (a[4] + 2 * a[5] * y[n] + a[6] * x[n]) * dy +
      y[n] * da[4] + y[n]^2 * da[5] + x[n] * y[n] * da[6])
  }
  c(dx, dy)
}

forecast_state <- function(xa, mfunc, ...) {
  mfunc(xa, ...)
}

forecast_ecov <- function(pamat, tfunc, sq, nx, ...) {
  n <- nrow(pamat)
  mpmat <- matrix(0, n, n)
  pfmat <- matrix(0, n, n)
  for (i in 1:n) {
    mpmat[, i] <- tfunc(pamat[1:nx, i], ..., da = pamat[(nx+1):n, i])
  }
  for (i in 1:n) {
    pfmat[, i] <- tfunc(mpmat[i, 1:nx], ..., da = mpmat[i, (nx+1):n])
  }
  q <- rnorm(nx, 0, sq)
  qmat <- outer(q, q)
  q <- rnorm(n - nx, 0, sq)
  qmat <- rbind(cbind(qmat, matrix(0, nx, n - nx)),
                cbind(matrix(0, n - nx, nx), outer(q, q)))
  pfmat + qmat
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

seed <- 514
set.seed(seed)

nmax <- 501
dt <- 0.001
at <- c(4, -2, -4, -6, 2, 4)
w1 <- c(1, 1, at)
iobs1 <- 10
dobs <- 10
tobs <- seq(iobs1, nmax, by = dobs)
ntobs <- length(tobs)
wt <- forward(w1, dt, nmax)
yo <- wt[, tobs]

a <- c(1, 0, 0, -1, 0, 0)
#a <- at
xa <- c(2, 2, a)
x_hist <- numeric(0)
y_hist <- numeric(0)
t_hist <- numeric(0)
xf_hist <- rep(0, ntobs)
yf_hist <- rep(0, ntobs)
xa_hist <- rep(0, ntobs)
ya_hist <- rep(0, ntobs)
a_hist <- matrix(0, 0, 6)
a_hist <- rbind(a_hist, a)

hfunc <- function(x) x
hmat <- cbind(diag(1, 2, 2), matrix(0, 2, 6))
sb <- 0.1
sr <- 0.1
sa <- 0.1
pamat <- rbind(cbind(diag(sb^2, 2, 2), matrix(0, 2, 6)),
               cbind(matrix(0, 6, 2), diag(sa^2, 6, 6)))
rmat <- diag(sr^2, 2, 2)
sq <- 0.2
n <- 1
#ntobs <- 1
for (t in 1:ntobs) {
  nf <- tobs[t] - n + 1
  t_hist <- c(t_hist, n:tobs[t], tobs[t])
  xf <- forecast_state(xa, forward, dt = dt, nmax = nf)
  pfmat <- forecast_ecov(pamat, tlm, sq, nx = 2,
                         x = xf[1,], y = xf[2,], dt = dt, a = xa[3:8])
  kmat <- calc_kgain(pfmat, hmat, rmat)
  xa <- analyze_state(xf[, nf], kmat, yo[, t], hfunc)
  xa[3:8] <- a
  pamat <- analyze_ecov(pfmat, kmat, hmat, rmat)
  x_hist <- c(x_hist, xf[1, ], xa[1])
  y_hist <- c(y_hist, xf[2, ], xa[2])
  a_hist <- rbind(a_hist, xa[3:8])
  xf_hist[t] <- xf[1, nf]
  yf_hist[t] <- xf[2, nf]
  xa_hist[t] <- xa[1]
  ya_hist[t] <- xa[2]
  n <- tobs[t]
}

plot(t_hist, x_hist, type = "l", lwd = 1, xlab = "t", ylab = "x, y",
     ylim = c(0, 2))
lines(t_hist, y_hist, lwd = 1, col = "red")
lines(1:nmax, wt[1,], lwd = 3, lty = 3)
lines(1:nmax, wt[2,], lwd = 3, lty = 3, col = "red")
legend("topright", c("x", "y", "xt", "yt"),
       lwd = c(1, 1, 2, 2), lty = c(1, 1, 3, 3),
       col = c("black", "red", "black", "red"))

tobs1 <- c(1, tobs[1:ntobs])
plot(tobs1, a_hist[, 1], type = "l", lwd = 1, xlab = "t", ylab = "x parameters",
     ylim = c(-4, 4))
lines(tobs1, a_hist[, 2], lwd = 1, col = "red")
lines(tobs1, a_hist[, 3], lwd = 1, col = "blue")
abline(h = at[1], lwd = 1, lty = 3)
abline(h = at[2], lwd = 1, lty = 3, col = "red")
abline(h = at[3], lwd = 1, lty = 3, col = "blue")
legend("topright", c("a1", "a2", "a3"),
       lwd = c(1, 1, 1), lty = c(1, 1, 1),
       col = c("black", "red", "blue"))

plot(tobs1, a_hist[, 4], type = "l", lwd = 1, xlab = "t", ylab = "y parameters",
     ylim = c(-6, 4))
lines(tobs1, a_hist[, 5], lwd = 1, col = "red")
lines(tobs1, a_hist[, 6], lwd = 1, col = "blue")
abline(h = at[4], lwd = 1, lty = 3)
abline(h = at[5], lwd = 1, lty = 3, col = "red")
abline(h = at[6], lwd = 1, lty = 3, col = "blue")
legend("topright", c("a4", "a5", "a6"),
       lwd = c(1, 1, 1), lty = c(1, 1, 1),
       col = c("black", "red", "blue"))