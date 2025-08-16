con <- file("out_ekf.dat", "rb")
nmax <- readBin(con, "integer", 1)
nhist <- readBin(con, "integer", 1)
xt <- readBin(con, "double", nmax)
yt <- readBin(con, "double", nmax)
t_hist <- readBin(con, "integer", nhist)
x_hist <- readBin(con, "double", nhist)
y_hist <- readBin(con, "double", nhist)
close(con)

plot(t_hist, x_hist, type = "l", lwd = 1, xlab = "t", ylab = "x, y",
     ylim = c(0, 2))
lines(t_hist, y_hist, lwd = 1, col = "red")
lines(1:nmax, xt, lwd = 3, lty = 3)
lines(1:nmax, yt, lwd = 3, lty = 3, col = "red")
legend("topright", c("x", "y", "xt", "yt"),
       lwd = c(1, 1, 2, 2), lty = c(1, 1, 3, 3),
       col = c("black", "red", "black", "red"))
