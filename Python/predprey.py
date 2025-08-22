import numpy as np 
import matplotlib.pyplot as plt 

nmax = 15001
dt = 0.001
x = np.zeros(nmax)
y = np.zeros(nmax)
a = np.array([4.0,-2.0,-4.0,-6.0,2.0,4.0])

x[0] = 1.0
y[0] = 1.0
for n in range(nmax-1):
    x[n+1] = x[n] + dt*(x[n]*(a[0] + a[1]*x[n] + a[2]*y[n]))
    y[n+1] = y[n] + dt*(y[n]*(a[3] + a[4]*y[n] + a[5]*x[n]))

plt.rcParams['font.size'] = 14
dpi = 144
fig, ax = plt.subplots(figsize=(8,4))
ax.plot(np.arange(nmax)*dt, x, c="blue", label="Prey(x)")
ax.plot(np.arange(nmax)*dt, y, c="red", label="Predator(y)")
ax.legend(loc="upper right")
ax.set_ylim(0.0,2.0)
ax.set_ylabel('Density')
ax.set_xlabel('time(days)')
fig.tight_layout()
fig.savefig("state.png", dpi=dpi)
plt.show()
plt.close()

fig, ax = plt.subplots(figsize=(6,6))
ax.plot(x, y)
ax.plot([0.0,2.0], [1.0,0.0], ls='dashed',c='black')
ax.plot([0.5,1.5], [2.0,0.0], ls='dashed',c='black')
ax.set_xlim(0.0,2.0)
ax.set_ylim(0.0,2.0)
ax.set_xlabel('Prey(x)')
ax.set_ylabel('Predator(y)')
fig.tight_layout()
fig.savefig("phase.png", dpi=dpi)
plt.show()
plt.close()