import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

cmap = plt.get_cmap('Set2')
#colors = ["#7570B3","#D95F02"]
colors = [cmap(i) for i in np.linspace(0, 1, 8)][1:3]
hex_colors = [mcolors.to_hex(c) for c in colors]
print(hex_colors)
colors = ['#8da0cb', '#fc8d62']

markers = ["o","s"]
def half_life(P0, a, b):
    p0 = P0/100
    tau = np.log((2*a-b*p0)/(a-b*p0))/a
    return tau

mu=0.6
eta=0#0.01
kappa=5e-5
D=100
alpha=1.15
gamma=1
c=1
dt=0.1
time = np.arange(0, 504, dt)
para=[1,mu,eta,kappa,D,alpha,gamma,c]
fig, ax = plt.subplots(1, 1, figsize=(2.05, 1.9))
popt = [0.6626225, 0.66151889]
pp = np.linspace(1,100,300)
halflife_theo=half_life(pp,*popt)
ax.plot(pp,halflife_theo,c='k',zorder=-10,lw=1)

for i,p0 in enumerate([50,100]):
    ax.scatter(p0,half_life(p0,*popt),color=colors[i],s=80, linewidth=1, edgecolors='black', marker=markers[i])

ax.set_xlim([10, 108])
ax.set_xticks([50, 100])
ax.set_ylim([0,11])
ax.set_yticks([])

ax.set_xlabel("P$_0$%")
ax.set_ylabel(r"$\tau_{1/2}$")
fig.subplots_adjust(left=0.25,right=0.95,bottom=0.23,top=0.95)
fig.savefig("./figures/pulsed_halflives.png",dpi=300)
fig.savefig("./figures/pulsed_halflives.svg")
plt.show()