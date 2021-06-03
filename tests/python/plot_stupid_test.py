import numpy as np
import matplotlib.pyplot as plt
import gwatpy.mcmc_routines as gmcmc
from corner import corner
def AC(data):
    h = 0
    N = len(data)
    acs = np.ones(int(N))
    mean = np.mean(data)
    for l in np.arange(int(N)):
        acs[l] = 1./(N-h) * np.sum( (data[h:]-mean)*(data[:N-h]-mean) )
        h +=1
    acs/=acs[0]
    return acs

data = gmcmc.trim_thin_file("data/output_stupid_test.hdf5")
print("Samples: ",len(data))
dim = len(data[0])

fig = corner(data,weights=np.ones(len(data))/len(data))
ax = np.array(fig.axes)
ax = ax.reshape((dim,dim))
for d in np.arange(dim):
    ax[d][d].axhline(1/20)
plt.savefig("plots/stupid_test.pdf")
plt.close()

samples  = 1e4

data = data[::int(len(data)/samples)]

acs  = None
fig, ax = plt.subplots(nrows=dim,ncols=1,figsize=[8,2*dim])
for d in np.arange(dim):
    acs = AC(data[:,d])
    acacs = AC(acs)
    ax[d].plot(acs[1:int(.9*len(acs))],label=str(d)+" AC",alpha=.7)
    ax[d].plot(acacs[1:int(.9*len(acs))],label=str(d)+" ACAC",alpha=.7)
    ax[d].legend()
plt.savefig("plots/AC_stupid_test.pdf")
plt.close()

