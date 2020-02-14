import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
import seaborn as sns; sns.set()

#curves = ["Aplus.csv","AdLIGODesign.csv","CE1.csv", "AdV_refsens_090427.csv", "AdV_baseline_sensitivity_12May09.csv", "CE2narrow.csv","../Hanford_O2_Strain.csv"]
curves = ["Aplus.csv","AdLIGODesign.csv","../Hanford_O2_Strain.csv"]

fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(10,5))
for curve in curves:
    freqs, noise = np.loadtxt(curve,delimiter=',',unpack=True)
    ax.loglog(freqs,noise,label=curve)

data = np.loadtxt("curves.txt",unpack=True)

ax.loglog(data[0],data[1], label="aLIGO",linestyle="--")
ax.loglog(data[0],data[2], label="A+",linestyle="--")
ax.set_ylim((10**(-24),10**(-20)))

plt.legend()
plt.savefig("noise_curve.pdf")
plt.close()
