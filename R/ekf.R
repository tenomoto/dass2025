predict_state <- function(w, dt, nmax) {
  x <- rep(0, nmax)
  y <- rep(0, nmax)
  x[1] <- w[1]
  y[1] <- w[2]
  a <- w[3:8]
  for (n in 1:(nmax-1)) {
    x[n + 1] <- x[n] + dt * (x[n] * (a[1] + a[2] * x[n] + a[3] * y[n]))
    y[n + 1] <- y[n] + dt * (y[n] * (a[4] + a[5] * y[n] + a[6] * x[n]))
  }
  rbind(x, y)
}

calc_jacobian <- function(w, dt) {
  x <- w[1]
  y <- w[2]
  a <- w[3:8]
  mmat <- matrix(0, 8, 8)
  mmat[1, 1] <- a[1] + 2 * a[2] * x + a[3] * y
  mmat[1, 2] <- a[3] * x
  mmat[1, 3:5] <- c(x, x^2, x * y)
  mmat[2, 1] <- a[6] * y
  mmat[2, 2] <- a[4] + 2 * a[5] * y + a[6] * x
  mmat[2, 6:8] <- c(y, y^2, x * y)
  diag(8) + dt * mmat
}

seed <- 514
set.seed(seed)
nmax <- 500
dt <- 0.001
at <- c(4, -2, -4, -6, 2, 4)
w1 <- c(1, 1, at)
wt <- predict_state(w1, dt, nmax)

sr <- 5e-2
iobs1 <- 10
dobs <- 10
tobs <- seq(iobs1, nmax, by = dobs)
ntobs <- length(tobs)
yo <- wt[, tobs] + matrix(rnorm(2 * ntobs, 0, sr), 2, ntobs)
rmat <- diag(sr^2, 2, 2)

a <- c(1, 0, 0, -1, 0, 0)
xa <- c(2, 2, a)

hfunc <- function(x) x
hmat <- diag(1, 2, 2)
sb <- 0.1
sq <- 0.03
pamat <- diag(sb^2, 2, 2)
qmat <- diag(sq^2, 2, 2)

x_hist <- numeric(0)
y_hist <- numeric(0)
t_hist <- numeric(0)
p_hist <- matrix(nrow = 0, ncol = 3)

n <- 0
for (t in 1:ntobs) {
  nf <- tobs[t] - n + 1
  t_hist <- c(t_hist, n:tobs[t])
  p_hist <- rbind(p_hist, c(pamat[1, 1], pamat[2, 2], pamat[1, 2]))
  xf <- predict_state(xa, dt = dt, nmax = nf)
  pfmat <- pamat
  for (i in 1:(nf - 1)) {
    mmat <- calc_jacobian(c(xf[, i], a), dt)[1:2, 1:2]
    pfmat <- mmat %*% pfmat %*% t(mmat) + qmat
    p_hist <- rbind(p_hist, c(pfmat[1, 1], pfmat[2, 2], pfmat[1, 2]))
  }
  kmat <- pfmat %*% t(hmat) %*% solve(hmat %*% pfmat %*% t(hmat) + rmat)
  xa <- xf[1:2, nf] + kmat %*% (yo[, t] - hfunc(xf[1:2, nf]))
  xa <- c(xa, a)
  ikhmat <- diag(nrow(pfmat)) - kmat %*% hmat
  pamat <- ikhmat %*% pfmat %*% t(ikhmat) + kmat %*% rmat %*% t(kmat)
  x_hist <- c(x_hist, xf[1, ])
  y_hist <- c(y_hist, xf[2, ])
  
  n <- tobs[t]
}
t_hist <- c(t_hist, n)
x_hist <- c(x_hist, xa[1])
y_hist <- c(y_hist, xa[2])
p_hist <- rbind(p_hist, c(pamat[1, 1], pamat[2, 2], pamat[1, 2]))

plot(t_hist, x_hist, type = "l", lwd = 1, xlab = "t", ylab = "x, y",
     ylim = c(0, 2))
lines(t_hist, y_hist, lwd = 1, col = "red")
lines(1:nmax, wt[1,], lwd = 3, lty = 3)
lines(1:nmax, wt[2,], lwd = 3, lty = 3, col = "red")
legend("topright", c("x", "y", "xt", "yt"),
       lwd = c(1, 1, 2, 2), lty = c(1, 1, 3, 3),
       col = c("black", "red", "black", "red"))

par(mar = c(2, 5, 2, 2))
plot(t_hist, p_hist[, 1], type = "l", lwd = 2, col = "blue",
     xlab = "t", ylab = "P", ylim = c(-0.02, 0.02),
     main = "EKF error covariance")
lines(t_hist, p_hist[, 2], lwd = 2, col = "red")
lines(t_hist, p_hist[, 3], lwd = 2, lty = 3, col = "purple")
lines(c(1,nmax), c(0.05*0.05,0.05*0.05), lwd=1, lty=5)
legend("bottomleft", c("Pxx", "Pyy", "Pxy=Pyx", "R"), cex = 1,
       lwd = c(2, 2, 5, 1), lty = c(1, 1, 3, 5),
       col = c("blue", "red", "purple", "black"))