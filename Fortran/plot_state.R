con <- file("xy.dat", "rb")
  nmax <- readBin(con, "integer", 1)
  dt <- readBin(con, "double", 1)
  x <- readBin(con, "double", nmax)
  y <- readBin(con, "double", nmax)
close(con)

png("state.png", 900, 450)
par(mar = c(5, 5, 2, 2))
t <- c(1:nmax)*dt
plot(t, x, type = "l", lwd = 2, col = "red", xlab = "time(days)", ylab = "Density",
    ylim = c(0, 2), 
    cex.main = 2, cex.lab = 2, cex.axis = 2)
lines(t, y, lwd = 2, col = "blue")
legend("topright", c("Prey(x)", "Predator(y)"), cex = 2, 
    lwd = c(2, 2), lty = c(1, 1), 
    col = c("red", "blue"))
dev.off()

png("phase.png", 600, 600)
par(mar = c(5, 5, 2, 2))
plot(x, y, type = "l", lwd = 2, xlab = "Prey(x)", ylab = "Predator(y)",
    xlim = c(0, 2), ylim = c(0, 2), 
    cex.main = 2, cex.lab = 2, cex.axis = 2)
lines(c(0, 2), c(1, 0), lwd = 1, lty = 3)
lines(c(0.5, 1.5), c(2, 0), lwd = 1, lty = 3)
dev.off()