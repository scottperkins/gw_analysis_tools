import numpy as np
import matplotlib.pyplot as plt
import gwatpy.gwatpy_plot as gp; gp.set()
import gwatpy.waveform_generator_ext as wf

data = np.loadtxt("data/snr_comp_stellar.csv",delimiter=',',unpack=True)
labels = ["SIMPS","GL_NONLOG","GL_LOG"]
#labels = ["GL_NONLOG","GL_LOG"]

temp= []
#temp.append(np.abs((data[0]-data[1])/data[0]))
temp.append(np.abs((data[0]-data[2])/data[0]))
temp.append(np.abs((data[0]-data[3])/data[0]))
resids = []
fig,ax = plt.subplots(nrows=2,ncols=1)
for x in temp:
    resids.append(x[~np.isnan(x)])
bins = np.logspace(np.log10(np.amin(resids)), np.log10(np.amax(resids)),100)

for x in np.arange(len(resids)):
    ax[0].hist(resids[x],bins=bins,label=labels[x])
ax[0].legend()
ax[0].set_xscale('log')
#plt.show()

times = [data[4],data[5],data[6],data[7]];
#times = [data[4],data[6],data[7]];
bins = np.logspace(np.log10(np.amin(times)), np.log10(np.amax(times)),100)
labels = ["GSL"] + labels
for x in np.arange(len(times)):
    ax[1].hist(times[x],bins=bins,label=labels[x])
ax[1].legend()
ax[1].set_xscale('log')
plt.savefig("plots/snr_comp_stellar.pdf")
plt.close()
