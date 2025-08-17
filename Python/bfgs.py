import numpy as np

def bfgs(par, fn, gr, args = (), **kwargs):
    control = {
        'maxit': 100,
        'debug': False,
        'ctol': 1e-10,
        'gtol': 1e-7,
        'stol': 1e-7
    }
    control.update(kwargs.pop('control', {}))

    if control['debug']:
        print(control)

    n = len(par)
    H = np.diag(np.ones(n))
    g = gr(par, *args)
    alpha = min(1.0, 1.0 / np.sum(np.abs(g)))
    new_par = par - alpha * g
    s = -alpha * g
    y = gr(new_par, *args) - g
    ys = np.dot(y, s)
    H = (ys / np.dot(y, y)) * np.diag(np.ones(n))

    if control['debug']:
        print("Initial Hessian:", H)

    cost_hist = np.zeros(control['maxit'] + 1)
    gnorm_hist = np.zeros(control['maxit'] + 1)
    par_hist = np.zeros((control['maxit'] + 1, n))

    convergence = -1

    for iter_ in range(control['maxit']):
        cost_hist[iter_] = fn(par, *args)
        
        gnorm = np.linalg.norm(g)
        gnorm_hist[iter_] = gnorm
        par_hist[iter_, :] = par

        if control['debug']:
            print(f"Iteration: {iter_ + 1}, Cost: {cost_hist[iter_]}, Gradient Norm: {gnorm}")

        if gnorm < control['gtol']:
            convergence = 0
            break
        if np.max(np.abs(s)) < control['stol']:
            convergence = 1
            break

        new_par = par - H @ g
        s = new_par - par
        new_gr = gr(new_par, *args)
        y = new_gr - g
        ys = np.dot(y, s)

        if ys > control['ctol']:
            rho = 1.0 / ys
            term1 = np.eye(n) - rho * np.outer(s, y)
            term2 = np.eye(n) - rho * np.outer(y, s)
            term3 = rho * np.outer(s, s)
            H = term1 @ H @ term2 + term3
            if control['debug']:
                print("Updated Hessian:", H)
        else:
            H = (ys / np.dot(y, y)) * np.diag(np.ones(n))
            if control['debug']:
                print("Reset Hessian:", H)

        if iter_ == control['maxit'] - 1:
            convergence = 2
        
        par = new_par
        g = new_gr

    return {
        'par': par_hist[:iter_ + 1],
        'convergence': convergence,
        'cost': cost_hist[:iter_ + 1],
        'gnorm': gnorm_hist[:iter_ + 1]
    }

if __name__ == '__main__':
    def rosenbrock_fn(x):
        return (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2

    def rosenbrock_gr(x):
        dfdx = -2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0]**2)
        dfdy = 200 * (x[1] - x[0]**2)
        return np.array([dfdx, dfdy])

    initial_par = np.array([-1.2, 1.0])
    result = bfgs(initial_par, fn=rosenbrock_fn, gr=rosenbrock_gr)
    print("Optimization Result:", result)
