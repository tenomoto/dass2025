import numpy as np
#from scipy.optimize import minimize
from bfgs import bfgs

def forward(dt, a, x1, y1, nmax):
    x = np.zeros(nmax)
    y = np.zeros(nmax)
    x[0] = x1
    y[0] = y1
    for n in range(nmax - 1):
        x[n + 1] = x[n] + dt * x[n] * (a[0] + a[1] * x[n] + a[2] * y[n])
        y[n + 1] = y[n] + dt * y[n] * (a[3] + a[4] * y[n] + a[5] * x[n])
    return x, y
  
def adjoint(dt, a, x, y, xo, yo, tobs) :
    nmax = x.size
    aa = np.zeros(6)
    ax = np.zeros(nmax)
    ay = np.zeros(nmax)
    for n in reversed(range(nmax - 1)):
        if n in tobs:
            ax[n] = ax[n] + (x[n] - xo[n])
            ay[n] = ay[n] + (y[n] - yo[n])
        aa[5] = aa[5] + dt * x[n] * y[n] * ay[n+1]
        aa[4] = aa[4] + dt * y[n] * y[n] * ay[n+1]
        aa[3] = aa[3] + dt * y[n] * ay[n+1]
        ax[n] = ax[n] + dt * a[5] * y[n] * ay[n+1]
        ay[n] = ay[n] + dt * a[4] * y[n] * ay[n+1]
        ay[n] = ay[n] + (1 + dt * (a[3] + a[4] * y[n] + a[5] * x[n])) * ay[n+1]
        aa[2] = aa[2] + dt * y[n] * x[n] * ax[n+1]
        aa[1] = aa[1] + dt * x[n] * x[n] * ax[n+1]
        aa[0] = aa[0] + dt * x[n] * ax[n+1]
        ay[n] = ay[n] + dt * a[2] * x[n] * ax[n+1]
        ax[n] = ax[n] + dt * a[1] * x[n] * ax[n+1]
        ax[n] = ax[n] + (1 + dt * (a[0] + a[1] * x[n] + a[2] * y[n])) * ax[n+1]
    return np.concatenate([aa, [ax[0], ay[0]]])

def fn(par, dt, nmax, xo, yo, tobs):
    x, y = forward(dt, par[0:6], par[6], par[7], nmax)
    return calc_cost(x[tobs], y[tobs], xo[tobs], yo[tobs])

def gr(par, dt, nmax, xo, yo, tobs):
  x, y = forward(dt, par[0:6], par[6], par[7], nmax)
  hist["par"].append(par)
  cost = calc_cost(x[tobs], y[tobs], xo[tobs], yo[tobs])
  hist["cost"].append(cost)
  grad = adjoint(dt, par[0:6], x, y, xo, yo, tobs)
  hist["gnorm"].append(np.sqrt(np.dot(grad, grad)))
  return grad

def calc_cost(xf, yf, xo, yo):
    return 0.5 * (np.sum((xf - xo)**2 + (yf - yo)**2)) 

if __name__ == "__main__":
    nmax = 501
    dt = 0.001
    at = [4, -2, -4, -6, 2, 4]
    x1, y1 = 1, 1
      
    tobs = np.arange(1, nmax, 2)
    xt, yt = forward(dt, at, x1, y1, nmax)
    xo, yo = xt, yt
     
    a = [1, 0, 0, -1, 0, 0]
    x1, y1 = 2, 2
    par = a + [x1, y1]
      
    cntl = {"maxiter": 100}
    hist = {"cost":[], "gnorm":[], "par":[]}
    args = (dt, nmax, xo, yo, tobs)
    alg = "BFGS"
    res = bfgs(par, fn, gr, args, control = cntl)
#    alg = "L-BFGS-B" # BFGS diverges
#    res = minimize(fn, par, args, method = alg, jac = gr, options = cntl)

    np.savez("out_adjoint.npz", alg = alg,
             cost = hist["cost"], gnorm = hist["gnorm"],
             par = hist["par"])

