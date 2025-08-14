bfgs <- function(par, fn, gr, ..., control){
  control <- modifyList(list(
    maxit = 100, debug = FALSE,
    ctol = 1e-10, gtol = 1e-7, stol = 1e-7), control)
 if (control$debug)  print(control)
  
  n <- length(par)
  H <- diag(n)
  g <- gr(par, ...)
  alpha <- min(1, 1. / sum(abs(g)))
  new_par <- par - alpha * g
  s <- -alpha * g
  y <- gr(new_par, ...) - g
  ys <- sum(y * s)
  H <- ys / sum(y * y) * diag(n)
  if (control$debug) cat("Initial Hessian:", H, "\n")
  
  cost <- rep(0, control$maxit + 1)
  gnorm.hist <- rep(0, control$maxit + 1)
  par.hist <- matrix(0, nrow = control$maxit + 1, ncol = n)

  for (iter in seq_len(control$maxit)) {
    cost[iter] <- fn(par, ...)
    
    gnorm <- sqrt(sum(g^2))
    gnorm.hist[iter] <- gnorm
    par.hist[iter, ] <- par

    if (control$debug) cat("Iteration:", iter, "Cost:", cost[iter], "Gradient Norm:", gnorm, "\n")   
    if (gnorm < control$gtol) {
      convergence <- 0
      break
    }
    if (max(abs(s)) < control$stol) {
      convergence <- 1
      break
    }

    new_par <- par - drop(H %*% g)
    s <- new_par - par
    new_gr <- gr(new_par, ...)
    y <- new_gr - g
    ys <- sum(y * s)
    if (ys > control$ctol) {
      rho <- 1 / ys
      H <- (diag(n) - rho * outer(s, y)) %*% H %*% (diag(n) - rho * outer(y, s)) + rho * outer(s, s)
      if (control$debug)  cat("Updated Hessian:", H, "\n")
    } else {
      H <- ys / sum(y * y) * diag(n)
      if (control$debug)  cat("Reset Hessian:", H, "\n")
    }
    if (iter == control$maxit) {
      convergence <- 2
    }
    
    par <- new_par
    g <- new_gr
  }
  
  list(par = par.hist[1:iter, ], convergence = convergence,
       cost = cost[1:iter], gnorm = gnorm.hist[1:iter])
}