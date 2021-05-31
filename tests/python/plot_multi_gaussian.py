import numpy as np
import gwatpy.mcmc_routines as gmcmc
import corner
import matplotlib.pyplot as plt
import emcee

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

lab=np.arange(10)


data = gmcmc.trim_thin_file("data/gaussian_output_0.hdf5",ac=None,trim=None,recalc_ac=True)
dim = len(data[0])
print("Samples: ",np.std(data[:,0]))
L1 = len(data)
print(len(data))
fig = corner.corner(data,show_titles=True, labels=lab)
plt.savefig("plots/gaussian_mcmc.pdf")
plt.close()

for x in np.arange(len(data[0])):
    print(emcee.autocorr.integrated_time(data[:,x],tol=0))
    plt.plot(data[:,x])
    plt.savefig("plots/gwat_trace_{}.pdf".format(x))
    #plt.show()
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
plt.savefig("plots/AC_gaussian.pdf")
plt.close()

#data = gmcmc.trim_thin_file("data/gaussian_output_1_small_E.hdf5",ac=None,trim=None,recalc_ac=False)
#print(np.std(data[:,0]))
#L1 = len(data)
#print(len(data))
#fig = corner.corner(data,show_titles=True, labels=lab)
##plt.savefig("plots/gaussian_mcmc_large_E.pdf")
##plt.close()
#
#data = gmcmc.trim_thin_file("data/gaussian_output_1.hdf5",ac=None,trim=None,recalc_ac=False)
#print(np.std(data[:,0]))
#L2 = len(data)
#print(len(data))
#corner.corner(data,show_titles=True, fig=fig,labels=lab,color='blue',weights=np.ones(len(data))*L1/L2)
#plt.savefig("plots/gaussian_mcmc_ensemble_comparison.pdf")
#plt.close()

#data = gmcmc.trim_thin_file("data/gaussian_output_19.hdf5",ac=None,trim=None,recalc_ac=False)
#print(np.std(data[:,0]))
#corner.corner(data,fig=fig,show_titles=True, labels=lab)
#plt.savefig("plots/gaussian_mcmc_comb_small_E.pdf")
#plt.close()

#lab=["w","x","y","z"]
#fig = corner.corner(data,show_titles=True, labels=lab)
#plt.savefig("plots/gaussian_mcmc_19.pdf")
#plt.close()

#data = gmcmc.trim_thin_file("data/gaussian_output_0.hdf5",ac=None,trim=None,recalc_ac=False)
#for x in np.arange(19):
#    t = gmcmc.trim_thin_file("data/gaussian_output_{}.hdf5".format(x+1),ac=None,trim=None,recalc_ac=False)
#    data = np.insert(data, [-1],t, axis=0)
#fig = corner.corner(data,show_titles=True, labels=lab)
#plt.savefig("plots/gaussian_mcmc_full_small_E.pdf")
#plt.close()
