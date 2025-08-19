alg <- "Steepest descent"
con <- file("out_adjoint.dat", "rb")
  iter <- readBin(con, "integer", 1)
  cost <- readBin(con, "double", iter)
  gnorm <- readBin(con, "double", iter)
  par <- matrix(readBin(con, "double", 2 * iter), nrow = 2)
close(con)

png("cost.png", 900, 450)
plot(log10(cost), type = "l", lwd = 2,
     main = paste("cost", alg), xlab = "Iteration", ylab = "J",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

png("gnorm.png", 900, 450)
plot(log10(gnorm), type = "l", lwd = 2,
     main = paste("gnorm", alg), xlab = "Iteration", ylab = "log10|g|",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

png("init.png", 900, 450)
plot(par[1, ], ylim = c(0, 2), type = "l", lwd = 2, xlab = "Iteration", ylab = "Initial conditions",
     main = paste("init", alg), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
lines(par[2, ], lwd = 2, col = "red")
legend("topright", legend = c("x1", "y1"),
       col = c("black", "red"), lwd = 2, cex = 1.5)
dev.off()
