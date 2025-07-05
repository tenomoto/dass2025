def adjoint(dt, a, x, y, xo, tobs) :
    nmax = x.size
    aa = np.zeros(6)
    ax = np.zeros(nmax)
    ay = np.zeros(nmax)
    for n in reversed(range(nmax-1)):
        aa[5] = aa[5] + dt * x[n] * y[n] * ay[n+1]
        aa[4] = aa[4] + dt * y[n] * y[n] * ay[n+1]
        aa[3] = aa[3] + dt * y[n] * ay[n+1]
        ax[n] <- ax[n] + dt * a[5] * y[n] * ay[n+1]
        ay[n] <- ay[n] + dt * a[4] * y[n] * ay[n+1]
        ay[n] <- ay[n] + (1 + dt * (a[3] + a[4] * y[n] + a[5] * x[n])) * ay[n+1]
        aa[2] = aa[2] + dt * y[n] * x[n] * ax[n+1]
        aa[1] = aa[1] + dt * x[n] * x[n] * ax[n+1]
        aa[0] = aa[0] + dt * x[n] * ax[n+1]
      　ay[n] = ay[n] + dt * a[2] * x[n] * ax[n+1]
    　　ax[n] = ax[n] + dt * a[1] * x[n] * ax[n+1]
    　　ax[n] = ax[n] + (1 + dt * (a[0] + a[1] * x[n] + a[2] * y[n])) * ax[n+1]
        if np.isin(n, tobs):
            ax[n] = ax[n] + (x[n] - xo[n])
            ay[n] = ay[n] + (y[n] - yo[n])
    return aa, ax[1], ay[1]
