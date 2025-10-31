import numpy as np
import matplotlib.pyplot as plt

cmap = plt.get_cmap('Set2')
colors = [cmap(i) for i in np.linspace(0, 1, 8)]
def half_life(P0, a, b):
    p0 = P0/100
    tau = np.log((2*a-b*p0)/(a-b*p0))/a
    return tau

time = np.arange(0,13,1)

host = "MG1655"
plasmid="pSC101"
inits = np.array([100, 99, 90, 50, 20])
fig,ax1=plt.subplots(1,1,figsize=(2.05, 1.9))

abundance = np.load(f"./LT_data_py/{plasmid}_mean.npy")

# define a hard floor for plasmid abundance (0.1%) for log scale plot
zero_idx = (abundance<0.1)
abundance[zero_idx]=0.1

n_ic = abundance.shape[0] # number of initial conditions
bio_rep = abundance.shape[1] # number of biological replicates

mean_dyn = np.mean(abundance,axis=1)
std_dyn = np.std(abundance,axis=1,ddof=1)
for j in range(n_ic):
    yk=mean_dyn[j,:]
    ykerr=std_dyn[j,:]
    ax1.plot(time,yk, color=colors[j], linewidth=1.5, zorder=-10)
    ax1.scatter(time,yk, s=80, color=colors[j], linewidth=1, edgecolors='black')
    ax1.errorbar(time,yk,ykerr,marker='None',linewidth=0,elinewidth=1,capsize=2,color="k")

ax1.set_xlim([-0.3, 10.3])
ax1.set_xticks([0, 5, 10])
ax1.set_ylim([1, 120])
ax1.set_yticks([1, 10, 100])
ax1.set_yscale("log")
ax1.set_xlabel("time (days)")
ax1.set_ylabel("P%",rotation=0,va="center",ha="right")
#ax1.text(x=0.1,y=0.1,s=f"MG1655\n{plasmid}",transform=ax1.transAxes)

fig.subplots_adjust(left=0.25,right=0.95,bottom=0.23,top=0.95)
fig.savefig("./figures/pSC101_dynamics.png",dpi=300)
fig.savefig("./figures/pSC101_dynamics.svg")
plt.show()