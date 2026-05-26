import numpy as np
from equations import func
from scipy.integrate import odeint
import matplotlib.pyplot as plt

np.random.seed(42)

cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]
# potential condition 1: stable vs unstable hosts
# Ns=10
# mu=0.7#np.random.uniform(0.5,0.7,Ns)
# stability=np.random.uniform(0,1,Ns)
# stability[0] = 0.1
# eta=10**(stability-4)#np.random.uniform(-4,-3,(Ns,Ns))
# kappa=10**(stability-5)#np.random.uniform(1e-5,1e-4,Ns)#5e-5
# D=0.2
# alpha=1.2 - stability*0.2#np.random.uniform(1.1,1.2,Ns)

# potential condition 2: donor strain has a large fitness disadvantage
Ns=10
mu=np.random.uniform(0.7,1,Ns)
mu[0] = 0.5
eta=10**np.random.uniform(-4,-3,(Ns,Ns))
kappa=5e-5
D=0.2
# alpha=np.random.uniform(1,1.1,Ns)
# alpha[0] = 1.1
alpha = 1.05
gamma = np.random.uniform(0.9,1,(Ns,Ns))
np.fill_diagonal(gamma, 1)
c=1

dt=0.1
time = np.arange(0,30*24,dt)
t0 = time/24
p0_list = np.linspace(0.01,1,20)
args1 = [Ns,mu,eta,kappa,D,alpha,gamma,c,0]

fig, axes = plt.subplots(1, 3, figsize=(6, 1.9))
ax1,ax2,ax3 = axes
fig4, ax4 = plt.subplots(1,1,figsize=(6,4))#figsize=(2.05, 1.9))
# simulation 1: strain 1's clonal dynamics
halflife_clonal = []
for i, p0 in enumerate(p0_list):
    init = np.zeros(2*Ns) # only strain 1 exists in the community
    init[0] = 1-p0
    init[Ns] = p0 # strain 1 as the only donor
    results = odeint(func, init, time, args=(args1,))
    p=np.sum(results[:,Ns:],axis=1)/np.sum(results,axis=1)*100
    ax1.plot(t0,p,c="k",lw=2)

    t_new = np.linspace(0,t0[-1],1000)
    p_new = np.interp(t_new, t0, p)
    hl = t_new[p_new<(p0*100/2)][0]
    halflife_clonal.append(hl)
ax4.plot(p0_list, halflife_clonal, c=colors[0], marker='o', label='clonal dynamics')

# simulation 2: community before antibiotic pulse;
# take strain 1 as the donor
halflife_before_pulse = []
for i, p0 in enumerate(p0_list):
    init = np.zeros(2*Ns)
    init[0] = 0
    init[Ns] = p0 # strain 1 as the only donor
    init[1:Ns] = (1-p0)/(Ns-1) # other strains as recipients
    results = odeint(func, init, time, args=(args1,))
    p=np.sum(results[:,Ns:],axis=1)/np.sum(results,axis=1)*100
    ax2.plot(t0,p,c="k",lw=2)
    t_new = np.linspace(0,t0[-1],1000)
    p_new = np.interp(t_new, t0, p)
    
    hl = t_new[p_new<(p0*100/2)][0]
    halflife_before_pulse.append(hl)
ax4.plot(p0_list, halflife_before_pulse, c=colors[1], marker='o', label='before pulse')

# simulation 3: after antibiotic pulse; 
# assume a uiform distribution of plasmid across different strains
# first simulate the community before pulse (strian 1 as the donor, with p0=0.5)
init = np.zeros(2*Ns)
init[0] = 0
init[Ns] = 0.1 # strain 1 as the only donor
init[1:Ns] = (1-0.1)/(Ns-1) # other strains as recipients
time1 = np.arange(0,2,dt)
results1 = odeint(func, init, time1, args=(args1,))
A = np.ones(Ns)
args2 = [Ns,mu,eta,kappa,D,alpha,gamma,c,A]
time2 = np.arange(time1[-1],4,dt)
results2 = odeint(func, results1[-1,:], time2, args=(args2,))
community_composition = results2[-1,0:Ns]+results2[-1,Ns:2*Ns]
community_composition = community_composition/np.sum(community_composition)
print(alpha)
print(community_composition)
print(np.sum(community_composition))
# for later simulations, fix the comunity composition
# and vary the plasmid abundance
halflife_after_pulse = []
for i, p0 in enumerate(p0_list):
    init = np.zeros(2*Ns)
    init[0:Ns] = (1-p0)*community_composition # community composition fixed, but total abundance varies
    init[Ns:2*Ns] = p0*community_composition 
    results = odeint(func, init, time, args=(args1,))
    p=np.sum(results[:,Ns:],axis=1)/np.sum(results,axis=1)*100
    ax3.plot(t0,p,c="k",lw=2)

    t_new = np.linspace(0,t0[-1],1000)
    p_new = np.interp(t_new, t0, p)

    hl = t_new[p_new<(p0*100/2)][0]
    halflife_after_pulse.append(hl)
ax4.plot(p0_list, halflife_after_pulse, c=colors[2], marker='o', label='after pulse')
ax4.legend()
plt.show()