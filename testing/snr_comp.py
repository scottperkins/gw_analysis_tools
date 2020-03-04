import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/snr_comp.csv",delimiter=',',unpack=True)
labels = ["SIMPS","GL_NONLOG","GL_LOG"]

temp= []
temp.append(np.abs((data[0]-data[1])/data[0]))
temp.append(np.abs((data[0]-data[2])/data[0]))
temp.append(np.abs((data[0]-data[3])/data[0]))
resids = []
for x in temp:
    resids.append(x[~np.isnan(x)])
bins = np.logspace(np.log10(np.amin(resids)), np.log10(np.amax(resids)),100)

for x in np.arange(len(resids)):
    plt.hist(resids[x],bins=bins,label=labels[x])
plt.legend()
plt.xscale('log')
plt.show()
plt.close()

times = [data[4],data[5],data[6],data[7]];
bins = np.logspace(np.log10(np.amin(times)), np.log10(np.amax(times)),100)
labels = ["GSL"] + labels
for x in np.arange(len(times)):
    plt.hist(times[x],bins=bins,label=labels[x])
plt.legend()
plt.xscale('log')
plt.show()
plt.close()
