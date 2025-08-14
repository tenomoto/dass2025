import numpy as np
import matplotlib.pyplot as plt
import sys

def predict_state(w, dt, nmax):
  x = np.zeros(nmax)
  y = np.zeros(nmax)
  x[0] = w[0]
  y[0] = w[1]
  a = w[2:8]
  for n in range(nmax-1):
    x[n + 1] = x[n] + dt * (x[n] * (a[0] + a[1] * x[n] + a[2] * y[n]))
    y[n + 1] = y[n] + dt * (y[n] * (a[3] + a[4] * y[n] + a[5] * x[n]))
  return np.vstack([x, y])

def calc_jacobian(w, dt):
  x = w[0]
  y = w[1]
  a = w[2:8]
  mmat = np.identity(8)
  mmat[0, 0] = a[0] + 2 * a[1] * x + a[2] * y
  mmat[0, 1] = a[2] * x
  mmat[0, 2:5] = np.array([x, x**2, x * y])
  mmat[1, 0] = a[5] * y
  mmat[1, 1] = a[3] + 2 * a[4] * y + a[5] * x
  mmat[1, 5:8] = np.array([y, y**2, x * y])
  return np.identity(8) + dt * mmat

def predict_covariance(mmat, pamat, qmat):
  return mmat @ pamat @ mmat.transpose() + qmat

def calc_gain(pfmat, hmat, rmat):
  return pfmat @ hmat.transpose() @ np.linalg.inv(hmat @ pfmat @ hmat.transpose() + rmat)

def analyze_state(w, kmat, yo, hfunc):
  return w + kmat @ (yo - hfunc(w))

def analyze_covariance(pfmat, kmat, hmat, rmat):
  n = pfmat.shape[0]
  ikhmat = np.identity(n) - kmat @ hmat
  return ikhmat @ pfmat @ ikhmat.transpose() + kmat @ rmat @ kmat.transpose()

seed = 514
rng = np.random.default_rng(seed)

nmax = 501
dt = 0.001
at = [4, -2, -4, -6, 2, 4]
w1 = [1, 1] + at
wt = predict_state(w1, dt, nmax)

sr = 1e-2
iobs1 = 10
dobs = 10
tobs = np.arange(iobs1 - 1, nmax, dobs)
ntobs = tobs.size
yo = wt[:, tobs] + rng.normal(0, sr, 2 * ntobs).reshape(2, ntobs)
rmat = np.diag(np.repeat(sr**2, 2))

a = [1, 0, 0, -1, 0, 0]
xa = np.array([2, 2] + a)

def hfunc(x): return x
hmat = np.identity(2)
sb = 0.1
sq = 0.1
pamat = np.diag(np.repeat(sb**2, 2))
qmat = np.diag(np.repeat(sq**2, 2))

x_hist = []
y_hist = []
t_hist = [] 

n = 0
for t in range(ntobs):
  nf = tobs[t] - n + 1
  t_hist += list(range(n, tobs[t] + 1)) + [tobs[t]]
  xf = predict_state(xa, dt, nf)
  mmat = np.identity(2)

  for i in range(nf - 1):
      mmat = calc_jacobian(np.concatenate([xf[:, i], a]), dt)[0:2, 0:2] @ mmat
  pfmat = predict_covariance(mmat, pamat, qmat)
  kmat = calc_gain(pfmat, hmat, rmat)
  xa = analyze_state(xf[0:2, nf - 1], kmat, yo[:, t], hfunc)
  xa = np.concatenate([xa, a])
  pamat = analyze_covariance(pfmat, kmat, hmat, rmat)
  x_hist = x_hist + xf[0, ].tolist() + [xa[0]]
  y_hist = y_hist + xf[1, ].tolist() + [xa[1]]
  n = tobs[t]

fig, ax = plt.subplots()
ax.plot(t_hist, x_hist, color = "black", label = "x")
ax.plot(t_hist, y_hist, color = "red", label = "y")
ax.plot(range(nmax), wt[0, :], color = "black", linewidth = 3, linestyle = "--", label = "xt")
ax.plot(range(nmax), wt[1, :], color = "red", linewidth = 3, linestyle = "--", label = "xt")
ax.set_xlabel("t")
ax.set_ylabel("x,y")
ax.set_ylim([0, 2])
ax.legend(loc = "upper right")
plt.show()
