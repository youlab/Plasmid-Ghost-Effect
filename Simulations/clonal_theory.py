import numpy as np
from equations import func
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]

def half_life(P0, a, b):
    p0 = P0/100
    tau = np.log((2*a-b*p0)/(a-b*p0))/a
    return tau

mu=0.6
eta=0
kappa=5e-5
D=0.2
alpha=1.1
gamma=1
c=1
A=0
dt=0.1
time = np.arange(0,30*24,dt)
t0 = time/24
para=[1,mu,eta,kappa,D,alpha,gamma,c,A]
fig, ax = plt.subplots(1, 1, figsize=(2.05, 1.9))
fig2, ax2 = plt.subplots(1,1,figsize=(2.05, 1.9))
halflife_sim=[]
p0_list=np.array([100,99,90,50,20])
for i,p0 in enumerate(p0_list):
    init = np.array([100 - p0, p0])/100
    results = odeint(func, init, time, args=(para,))
    p=results[:,1]/np.sum(results,axis=1)*100
    ax.plot(t0,p,c=colors[i],lw=2)

    t_new = np.linspace(0,t0[-1],1000)
    p_new = np.interp(t_new, t0, p)
    hl = t_new[p_new<(p0/2)][0]
    halflife_sim.append(hl)

ax.set_xlabel("time (days)")
ax.set_ylabel("P%",rotation=0,va="center",ha="right")
ax.set_xlim([0, 23])
ax.set_xticks([0,10,20])
ax.set_yscale("log")
ax.set_ylim([1,120])
fig.subplots_adjust(left=0.25,right=0.95,bottom=0.23,top=0.95)
fig.savefig("./figures/plasmid_dynamics_continuous_simulation.png",dpi=300)
fig.savefig("./figures/plasmid_dynamics_continuous_simulation.svg")

popt, pcov = curve_fit(half_life, p0_list, halflife_sim, p0=(0.7, 0.69))
print(popt)
pp = np.linspace(1,100,300)
halflife_theo=half_life(pp,*popt)
ax2.plot(pp,halflife_theo,c='k',zorder=-10,lw=1)
for i,p0 in enumerate(p0_list):
    ax2.scatter(p0,halflife_sim[i],color=colors[i],s=80, linewidth=1, edgecolors='black')

ax2.set_xlim([10, 105])
ax2.set_xticks([50, 100])
ax2.set_ylim([0,14])
ax2.set_yticks([0,6,12])

ax2.set_xlabel("P$_0$%")
ax2.set_ylabel(r"$\tau_{1/2}$"+" (days)")
fig2.subplots_adjust(left=0.25,right=0.95,bottom=0.23,top=0.95)
fig2.savefig("./figures/half_life_dependence_continuous.png",dpi=300)
fig2.savefig("./figures/half_life_dependence_continuous.svg")
plt.show()