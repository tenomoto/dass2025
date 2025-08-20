import numpy as np
import matplotlib.pyplot as plt

with np.load("out_adjoint.npz") as data:
    alg = data["alg"]
    cost = data["cost"]
    gnorm = data["gnorm"]
    par = data["par"]

fig_size = (9, 4.5)
dpi = 144

fig, ax = plt.subplots(figsize=fig_size)
ax.plot(np.log10(cost), lw=2)
ax.set_title(f"cost {alg}", fontsize=15)
ax.set_xlabel("Iteration", fontsize=15)
ax.set_ylabel("log10|g|", fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
fig.tight_layout()
fig.savefig("cost.png", dpi=dpi)
plt.close(fig)

fig, ax = plt.subplots(figsize=fig_size)
ax.plot(np.log10(gnorm), lw=2, color = "black")
ax.set_title(f"gnorm {alg}", fontsize=15)
ax.set_xlabel("Iteration", fontsize=15)
ax.set_ylabel("log10|g|", fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
fig.tight_layout()
fig.savefig("gnorm.png", dpi=dpi)
plt.close(fig)

fig, ax = plt.subplots(figsize=fig_size)
ax.plot(par[:, 6], label="x1", lw=2, color="black")
ax.plot(par[:, 7], label="y1", lw=2, color="red")
ax.set_title(f"init {alg}", fontsize=15)
ax.set_xlabel("Iteration", fontsize=15)
ax.set_ylabel("Initial conditions", fontsize=15)
ax.set_ylim(0, 2)
ax.legend(loc="upper right", fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
fig.tight_layout()
fig.savefig("init.png", dpi=dpi)
plt.close(fig)

fig, ax = plt.subplots(figsize=fig_size)
ax.plot(par[:, 0], label="a1", lw=2, color="black")
ax.plot(par[:, 1], label="a2", lw=2, color="red")
ax.plot(par[:, 2], label="a3", lw=2, color="blue")
ax.set_ylim(-10, 10)
ax.set_title(f"xparam {alg}", fontsize=15)
ax.set_xlabel("Iteration", fontsize=15)
ax.set_ylabel("x parameters", fontsize=15)
ax.legend(loc="upper left", fontsize=10)
ax.tick_params(axis='both', which='major', labelsize=15)
fig.tight_layout()
fig.savefig("xparam.png", dpi=dpi)
plt.close(fig)

fig, ax = plt.subplots(figsize=fig_size)
ax.plot(par[:, 3], label="a4", lw=2, color="black")
ax.plot(par[:, 4], label="a5", lw=2, color="red")
ax.plot(par[:, 5], label="a6", lw=2, color="blue")
ax.set_ylim(-10, 10)
ax.set_title(f"yparam {alg}", fontsize=15)
ax.set_xlabel("Iteration", fontsize=15)
ax.set_ylabel("y parameters", fontsize=15)
ax.legend(loc="upper left", fontsize=10)
ax.tick_params(axis='both', which='major', labelsize=15)
fig.tight_layout()
fig.savefig("yparam.png", dpi=dpi)
plt.close(fig)

