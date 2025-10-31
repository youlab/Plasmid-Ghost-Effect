import numpy as np
import pandas as pd
import natsort as ns # install by using 'pip install natsort' #ref: https://pypi.org/project/natsort/
import matplotlib.pyplot as plt
import seaborn as sns

I2=np.arange(1,25,1)
I1=np.arange(1,17,1)
f = pd.read_csv("./raw_data/1/1.tabular", sep="\t")
sorted_id=ns.index_natsorted(f["# Barcode"])
column0=np.array(f["# Barcode"])[sorted_id]
df={"Barcodes":column0}
#align the files by column
Reads=np.zeros((len(I2),len(I1)))
for i, i2 in enumerate(I2):
    for j, i1 in enumerate(I1):
        try:
            f = pd.read_csv("%i/%i.tabular"%(i1,i2), sep="\t")
            counts = np.array(f["Count"], dtype=int)
            sorted_id = ns.index_natsorted(f["# Barcode"])
            df["%i-%i" % (i1, i2)] = counts[sorted_id]  # using np.array to enforce the counts to be numbers, instead of strings
            total = counts[-1]
        except:
            df["%i-%i" % (i1, i2)] = np.zeros(len(column0))
            total = 1
        Reads[i,j]=total
df=pd.DataFrame(df)
df.to_csv("./processed_data/merged_counts.csv")

fig,ax=plt.subplots(1,1,figsize=(6,2.5))
Reads=np.log10(Reads)
#print(Reads[12,10])
sns.heatmap(Reads.T,ax=ax,xticklabels=I2,yticklabels=I1,square=True,cbar_kws={"label":"reads"},vmin=2)

ax.set_xlabel("P7")
ax.set_ylabel("P5")
fig.savefig("./processed_data/merged_counts.png",dpi=300)
plt.show()
