import numpy as np
import matplotlib.pyplot as plt

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
    return np.array([ax[0], ay[0]])

def fn(par, dt, nmax, xo, yo, tobs):
    x, y = forward(dt, par[0:6], par[6], par[7], nmax)
    return calc_cost(x[tobs], y[tobs], xo[tobs], yo[tobs])

def gr(par, dt, nmax, xo, yo, tobs):
    x, y = forward(dt, par[0:6], par[6], par[7], nmax)
    hist["par"].append(par.copy())
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
     
    #a = [1, 0, 0, -1, 0, 0] #inexact parameters
    a = at.copy()
    x1, y1 = 1.5, 1.5
    par = np.concatenate([a, [x1, y1]])
      
    alg = "steepest descent"
    maxiter = 200
    alpha = 1.0e-3
    gtol = 1.0e-5
    hist = {"cost":[], "gnorm":[], "par":[]}

    nplot = 1
    maxplot = 10
    tplot_state = []
    dpi = 144
    cmap = plt.get_cmap('tab10')
    fig, ax = plt.subplots(figsize=(6,6))
    ax.plot(xt,yt,c='black',lw=3.0,label='nature')
    ax.plot(xt[0],yt[0],c='black',marker='*',ms=10)
    
    niter = 0
    while (niter < maxiter):
        # evaluate cost function & calculate gradient
        g = gr(par, dt, nmax, xo, yo, tobs)
        print("{} iterations, cost = {:.4e}, gradient norm = {:.4e}"\
            .format(niter,hist["cost"][-1],hist["gnorm"][-1]))
        print("x, y = {}, {}"\
            .format(hist["par"][-1][6],hist["par"][-1][7]))
        if niter%nplot==0 and len(tplot_state) < maxplot:
            x, y = forward(dt, par[0:6], par[6], par[7], nmax)
            ax.plot(x,y,ls='--',c=cmap(len(tplot_state)),label=f'iter{niter}')
            ax.plot(x[0],y[0],c=cmap(len(tplot_state)),marker='*')
            tplot_state.append(niter)
        # check convergence
        if hist["gnorm"][-1] < gtol:
            print("Convergence: {} iterations".format(niter))
            print("Final cost = {:.4e}, gradient norm = {:.4e}"\
                .format(hist["cost"][-1],hist["gnorm"][-1]))
            print("x, y = {}, {}"\
                .format(hist["par"][-1][6],hist["par"][-1][7]))
            break
        # update control variables
        par[6:8] = par[6:8] - alpha*g
        niter += 1
    print(tplot_state)
    x, y = forward(dt, par[0:6], par[6], par[7], nmax)
    ax.plot(x,y,c='r',ls='--',lw=3.0,label='final')
    ax.plot(x[0],y[0],c='r',marker='*')
    ax.plot([0.0,2.0],[1.0,0.0],ls=':',c='k',zorder=0)
    ax.plot([0.5,1.5],[2.0,0.0],ls=':',c='k',zorder=0)
    ax.set_xlim(0.0,2.0)
    ax.set_ylim(0.0,2.0)
    ax.grid()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.legend()
    fig.tight_layout()
    fig.savefig("phase.png", dpi=dpi)
    plt.close(fig)
    
    fig_size = (9, 4.5)

    fig, ax = plt.subplots(figsize=fig_size)
    ax.plot(np.log10(hist["cost"]), lw=2)
    ax.set_title(f"cost {alg}", fontsize=15)
    ax.set_xlabel("Iteration", fontsize=15)
    ax.set_ylabel(r"log$_{10}$ J", fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    fig.tight_layout()
    fig.savefig("cost.png", dpi=dpi)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=fig_size)
    ax.plot(np.log10(hist["gnorm"]), lw=2, color = "black")
    ax.set_title(f"gnorm {alg}", fontsize=15)
    ax.set_xlabel("Iteration", fontsize=15)
    ax.set_ylabel(r"log$_{10}$|g|", fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    fig.tight_layout()
    fig.savefig("gnorm.png", dpi=dpi)
    plt.close(fig)

    par = np.array(hist["par"])
    fig, ax = plt.subplots(figsize=fig_size)
    ax.plot(par[:, 6], label=r"x$_0$", lw=2, color="black")
    ax.plot(par[:, 7], label=r"y$_0$", lw=2, color="red")
    ax.set_title(f"init {alg}", fontsize=15)
    ax.set_xlabel("Iteration", fontsize=15)
    ax.set_ylabel("Initial conditions", fontsize=15)
    ax.set_ylim(0, 2)
    ax.legend(loc="upper right", fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    fig.tight_layout()
    fig.savefig("init.png", dpi=dpi)
    plt.close(fig)
