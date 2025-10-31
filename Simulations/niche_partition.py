import numpy as np
from equations import func
from scipy.integrate import odeint
import matplotlib.pyplot as plt

colors = ["#8DA0CB","#FC8D62"]

Ns=2
mu=np.array([0.6,0.6])
eta=np.zeros((Ns,Ns))
kappa=5e-5
D=0.2
alpha=np.array([1.1,1.1])
gij=0
gamma=np.array([[1,gij],[gij,1]])
c=1

dt=0.1

init_p=0.3
init = np.array([0.5-init_p, 0.5,init_p, 0])

t1=3*24
t2=4*24
t3=31*24

fig1,ax1=plt.subplots(1,1,figsize=(2.05, 1.9))

args1 = [Ns,mu,eta,kappa,D,alpha,gamma,c,0]
time0 = np.arange(0,t3,dt)
results0 = odeint(func, init, time0, args=(args1,))
p_tot0 = np.sum(results0[:, Ns:],axis=1) / np.sum(results0,axis=1)*100
time0 = time0/24

time1 = np.arange(0,t1,dt)
results1 = odeint(func, init, time1, args=(args1,))
A = np.array([1,0.3])
args2 = [Ns,mu,eta,kappa,D,alpha,gamma,c,A]
time2 = np.arange(time1[-1],t2,dt)
results2 = odeint(func, results1[-1,:], time2, args=(args2,))
time3 = np.arange(time2[-1],t3,dt)
results3 = odeint(func, results2[-1,:], time3, args=(args1,))

results = np.vstack((results1, results2[1:,:], results3[1:,:]))
time = np.hstack((time1, time2[1:], time3[1:]))/24

p_tot = np.sum(results[:, Ns:],axis=1) / np.sum(results,axis=1)*100
p = results[:,Ns:]
s = results[:, 0:Ns]+results[:, Ns:]
tot = np.sum(s, axis=1)

ax1.plot(time0,p_tot0,c=colors[0],linewidth=2,zorder=-10)
ax1.plot(time,p_tot,c=colors[1],linewidth=2,zorder=-10)
ax1.fill_between(x=[t1/24, t2/24], y1=np.ones(2), y2=np.ones(2) * 120, facecolor="#808080", alpha=0.3, zorder=-10)
ax1.set_xticks([0, 10, 20])
ax1.set_xlim([0, 20])
ax1.set_yscale("log")
ax1.set_ylim([1,120])
ax1.set_xlabel("time (days)")
ax1.set_ylabel("P%",rotation=0,va="center",ha="right")
fig1.subplots_adjust(left=0.25,right=0.95,bottom=0.23,top=0.95)
fig1.savefig("./figures/niche_partition_non_mobilizable.png",dpi=300)
fig1.savefig("./figures/niche_partition_non_mobilizable.svg")

fig2,ax2=plt.subplots(1,1,figsize=(2.05, 1.9))
s1 = s[:,0]/tot
p1 = p[:,0]/tot
ax2.plot(time,s1,c="k",linewidth=1.5)
ax2.fill_between(x=time, y1=np.zeros(len(time)), y2=p1, color="#FF7F0E", alpha=1)
ax2.fill_between(x=[t1/24, t2/24], y1=np.zeros(2), y2=np.ones(2), facecolor="#808080", alpha=0.3, zorder=-10)
ax2.set_xlim([0,20])
ax2.set_ylim([0,1])
ax2.set_xticks([0,10,20])
ax2.set_yticks([])
ax2.set_xlabel("time (days)")
ax2.set_ylabel("composition")
fig2.subplots_adjust(left=0.25,right=0.95,bottom=0.23,top=0.95)
fig2.savefig("./figures/niche_partition_composition.png",dpi=300)
fig2.savefig("./figures/niche_partition_composition.svg")
plt.show()
