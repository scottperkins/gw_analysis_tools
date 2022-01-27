import numpy as np
import h5py 
from corner import corner 
import matplotlib.pyplot as plt

injections = np.loadtxt("data/injections_PTRJMCMC.csv",delimiter=',')
injections[6] = np.exp(injections[6])
injections[7] = np.exp(injections[7])
print(injections)
dataFile = h5py.File("data/PTRJMCMC_GW_injection_output.hdf5",'r')
first=True
data = None
print(dataFile.keys())
ACs = np.amax(dataFile["MCMC_METADATA"]["AC VALUES"],axis=1)
print(ACs)
ct = 0
for chain in dataFile["MCMC_OUTPUT"].keys():
    if "CHAIN" in chain:
        if first:
            data = np.array(dataFile["MCMC_OUTPUT"][chain])[::ACs[ct]]
            first = False
        else:
            data = np.insert(data,len(data), np.array(dataFile["MCMC_OUTPUT"][chain])[::ACs[ct]],axis=0)
    ct+=1
data = np.array(data)
data[:,7] = np.exp(data[:,7])
data[:,6] = np.exp(data[:,6])
#plt.hist(data,bins=50,histtype='stepfilled')
fig = corner(data)
dim = len(data[0])
axes = np.array(fig.axes).reshape(dim,dim)
for i in np.arange(dim):
    ax = axes[i,i]
    ax.axvline(injections[i])

for yi in np.arange(dim):
    for xi in np.arange(yi):
        ax = axes[yi,xi]
        ax.axvline(injections[xi])
        ax.axhline(injections[yi])
        ax.plot(injections[xi],injections[yi])


plt.savefig("plots/PTRJMCMC_sampler.pdf")
plt.close()


dataFile = h5py.File("data/PTRJMCMC_GW_injectionPrior_output.hdf5",'r')
first=True
data = None
ct = 0
for chain in dataFile["MCMC_OUTPUT"].keys():
    if "CHAIN" in chain:
        if first:
            data = np.array(dataFile["MCMC_OUTPUT"][chain])
            first = False
        else:
            data = np.insert(data,len(data), np.array(dataFile["MCMC_OUTPUT"][chain]),axis=0)
    ct+=1
data = np.array(data)
data[:,7] = np.exp(data[:,7])
data[:,6] = np.exp(data[:,6])
#plt.hist(data,bins=50,histtype='stepfilled')
fig = corner(data)

plt.savefig("plots/PTRJMCMC_sampler_prior.pdf")
plt.close()


