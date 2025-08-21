load("out_adjoint.RData")

png("cost.png", 900, 450)
plot(log10(hist$cost), type = "l", lwd = 2,
     main = paste("cost", alg), xlab = "Iteration", ylab = "log10(J)",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

png("gnorm.png", 900, 450)
plot(log10(hist$gnorm), type = "l", lwd = 2,
     main = paste("gnorm", alg), xlab = "Iteration", ylab = "log10|g|",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

png("init.png", 900, 450)
plot(hist$par[, 7], ylim = c(0, 2), type = "l", lwd = 2, xlab = "Iteration", ylab = "Initial conditions",
     main = paste("init", alg), cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
lines(hist$par[, 8], lwd = 2, col = "red")
legend("topright", legend = c("x1", "y1"),
       col = c("black", "red"), lwd = 2, cex = 1.5)
dev.off()

png("xparam.png", 900, 450)
plot(hist$par[, 1], ylim = c(-10, 10), type = "l", lwd = 2,
     main = paste("xparam", alg), xlab = "Iteration", ylab = "x parameters",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
lines(hist$par[, 2], lwd = 2, col = "red")
lines(hist$par[, 3], lwd = 2, col = "blue")
legend("topleft", legend = c("a1", "a2", "a3"),
       col = c("black", "red", "blue"), lwd = 2, cex = 1)
dev.off()

png("yparam.png", 900, 450)
plot(hist$par[, 4], ylim = c(-10, 10), type = "l", lwd = 2, 
     main = paste("yparam", alg), xlab = "Iteration", ylab = "y parameters",
     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
lines(hist$par[, 5], lwd = 2, col = "red")
lines(hist$par[, 6], lwd = 2, col = "blue")
legend("topleft", legend = c("a4", "a5", "a6"),
       col = c("black", "red", "blue"), lwd = 2, cex = 1)
dev.off()
