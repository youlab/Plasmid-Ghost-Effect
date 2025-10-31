import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plasmid = "R388"
Ab = "Trim"
dilution = 100

time=np.array([0,2,3,7,10,15,20,25,30,35])
Conditions=["NoAb","Ab",]
Labels=["no pulse","+Trim",]
Reps=3
colors = ['#8da0cb', '#fc8d62']

fig, ax = plt.subplots(1, 1, figsize=(2.05, 1.9))
simpson = np.zeros((2,Reps, len(time)))
for i in range(len(Conditions)):
    for j in range(Reps):
        df = pd.read_excel("./processed_data/R388_100_dilution_exclude_zeros.xlsx",
                           sheet_name=f"{Conditions[i]}_Rep{j+1}",index_col=0)
        df=np.array(df)
        comp = df/np.sum(df,axis=0)
        div = 1/np.sum(comp**2,axis=0)
        simpson[i,j,:] = div
        ax.plot(time, div, color=colors[i], linewidth=1.5,zorder=0)
        ax.scatter(time, div, s=30, color=colors[i], linewidth=0.5, edgecolors='black',zorder=10)
ax.set_yscale("log")
plt.show()
