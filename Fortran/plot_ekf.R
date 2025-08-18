con <- file("out_ekf.dat", "rb")
  nmax <- readBin(con, "integer", 1)
  nhist <- readBin(con, "integer", 1)
  xt <- readBin(con, "double", nmax)
  yt <- readBin(con, "double", nmax)
  t_hist <- readBin(con, "integer", nhist)
  x_hist <- readBin(con, "double", nhist)
  y_hist <- readBin(con, "double", nhist)
  p11_hist <- readBin(con, "double", nhist)
  p22_hist <- readBin(con, "double", nhist)
  p12_hist <- readBin(con, "double", nhist)
  p21_hist <- readBin(con, "double", nhist)
close(con)


png("state.png", 900, 450)
plot(t_hist, x_hist, type = "l", lwd = 2, xlab = "t", ylab = "x, y",
     ylim = c(0, 2), main = "EKF", cex.main = 2, cex.lab = 2, cex.axis =2)
lines(t_hist, y_hist, lwd = 2, col = "red")
lines(1:nmax, xt, lwd = 5, lty = 3)
lines(1:nmax, yt, lwd = 5, lty = 3, col = "red")
legend("topright", c("x", "y", "xt", "yt"), cex = 2, 
       lwd = c(2, 2, 5, 5), lty = c(1, 1, 3, 3),
       col = c("black", "red", "black", "red"))
dev.off()

png("errcov.png", 900, 450)
plot(t_hist, p11_hist, type = "l", lwd = 2, col = "blue", xlab = "t", ylab = "P",
     ylim = c(-0.02, 0.02), main = "EKF error covariance", cex.main = 2, cex.lab = 2, cex.axis =2)
lines(t_hist, p22_hist, lwd = 2, col = "red")
lines(t_hist, p12_hist, lwd = 5, lty = 3, col = "blue")
#lines(t_hist, p21_hist, lwd = 5, lty = 3, col = "red")
lines(c(1,nmax), c(0.05*0.05,0.05*0.05), lwd=1, lty=5)
legend("bottomleft", c("Pxx", "Pyy", "Pxy=Pyx", "R"), cex = 2,
       lwd = c(2, 2, 5, 1), lty = c(1, 1, 3, 5),
       col = c("blue", "red", "blue", "black"))
dev.off()
