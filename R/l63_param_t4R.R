library(torch)

calc_tendency <- function(a, x, y, z) {
  r <- a[1]
  s <- a[2]
  b <- a[3]
  dxdt <- -s * x + s * y
  dydt <- -x * z + r * x - y
  dzdt <- x * y - b * z
  torch_stack(c(dxdt, dydt, dzdt), dim = -1)
}

calc.cost <- function(wf, wo) {
  0.5 * (torch_sum((wf - wo)^2))
}

dev <- "cpu"
dtype <- torch_double()
nmax <- 200
dt <- 0.01
loss_old <- Inf
stol <- 1e-8

at <- c(32, 10, 8/3)
w1t <- c(1, 3, 5) 
wt <- torch_zeros(nmax, 3, device = dev, dtype = dtype)
wt[1, 1:3] <- w1t
for (n in 1:(nmax-1)) {
  dwdt <- calc_tendency(torch_tensor(at, device = dev, dtype = dtype), 
                        wt[n, 1], wt[n, 2], wt[n, 3])
  wt[n+1, ] <- wt[n, ] + dt * dwdt
}

tobs <- seq(60, nmax, by = 60)
wo <- wt[tobs]

a <- c(30, 11, 2)
w1 <- c(1.1, 3.3, 5.5)
par <- torch_tensor(c(a, w1), requires_grad = TRUE,
                    device = dev, dtype = dtype)
print(par$dtype)
print(par$device)
hist <- list(cost = numeric(0), gnorm = numeric(0), par = vector(length=0))

optimizer <- optim_lbfgs(par, max_iter = 20, tolerance_change = 1e-9, tolerance_grad = 1e-5, line_search_fn = "strong_wolfe")

maxit <- 10
for (i in 1:maxit) {
  closure <- function() {
    hist$par <<- rbind(hist$par, as.numeric(par))
    optimizer$zero_grad()
    w_hist <- list()
    w_curr <- par[4:6]$view(c(1, 3))
    w_hist[[1]] <- w_curr
    for (n in 1:(nmax-1)) {
      dwdt <- calc_tendency(par[1:3],
                            w_curr[1, 1], w_curr[1, 2], w_curr[1, 3])
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
  
  x <- as.numeric(par[4])
  y <- as.numeric(par[5])
  z <- as.numeric(par[6])
  cat("Iteration:", i, "loss=", as.numeric(loss), "x, y, z=", x, y, z, "\n")
  if (abs(loss_old - as.numeric(loss)) < stol) {
    cat("Converged at iteration", i, "\n")
    break 
  }
  loss_old <- as.numeric(loss)
}

alg <- "lbfgs"
cat("Optimized Parameters:", as.numeric(par), "\n")
#save(alg, hist, file = "out_l63_param_t4R.RData")

hgnorm <- as.numeric(hist$gnorm)
hcost <- as.numeric(hist$cost)
hpar <- as.matrix(hist$par)

#png("cost_l63_param_t4R.png", 900, 450)
plot(log10(hcost), type = "l", lwd = 2,
     main = paste("cost", alg), xlab = "Iteration", ylab = "log10(J)",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
#dev.off()

#png("gnorm_l63_param_t4R.png", 900, 450)
plot(log10(hgnorm), type = "l", lwd = 2,
     main = paste("gnorm", alg), xlab = "Iteration", ylab = "log10|g|",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
#dev.off()

#png("init_l63_param_t4R.png", 900, 450)
plot(hpar[, 4], ylim = c(0, 8), type = "l", lwd = 2, xlab = "Iteration", ylab = "Initial conditions",
     main = paste("init", alg), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
lines(hpar[, 5], lwd = 2, col = "red")
lines(hpar[, 6], lwd = 2, col = "blue")
abline(h = w1t, lty = 2, col = c("black", "red", "blue"))
legend("topright", legend = c("X", "Y", "Z"),
       col = c("black", "red", "blue"), lwd = 2, cex = 1)
#dev.off()

#png("param_l63_param_t4R.png", 900, 450)
plot(hpar[, 1] - at[1], ylim = c(-2, 2), type = "l", lwd = 2,
     main = paste("param", alg), xlab = "Iteration", ylab = "parameter error",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
lines(hpar[, 2] - at[2], lwd = 2, col = "red")
lines(hpar[, 3] - at[3], lwd = 2, col = "blue")
abline(h = at, lty = 2, col = c("black", "blue", "red"))
legend("topright", legend = c(expression("r" - r[t]), expression(sigma-sigma[t]), expression(beta-beta[t])),
       col = c("black", "red", "blue"), lwd = 2, cex = 1)
#dev.off()