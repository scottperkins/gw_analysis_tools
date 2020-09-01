
import numpy as np
import gwatpy.mcmc_routines as gmcmc
import gwatpy.plot as gpp; gpp.set()
import corner
import matplotlib.pyplot as plt
import h5py
import scipy


def fn(x, xpts):
    return x[0]*np.sin(x[2]*xpts + x[1])\
    *1./np.sqrt(2.*x[3]*x[3]) \
    * np.exp(-.5 * (x[4]-xpts)*(x[4]-xpts)/(x[3]*x[3])) * (1+ scipy.special.erf( (x[5]*(xpts-x[4])/np.sqrt(2))))
def fn_sym(x, xpts):
    return x[0]*np.sin(x[2]*xpts + x[1])\
    *1./np.sqrt(2.*x[3]*x[3]) \
    * np.exp(-.5 * (x[4]-xpts)*(x[4]-xpts)/(x[3]*x[3]))

injections = np.loadtxt("data/RJ_injections_sin.csv",delimiter=',')
max_dim = int(len(injections))
injected_data = np.loadtxt("data/RJ_sin_data_injection.csv",delimiter=',')
injected_data_clean = np.loadtxt("data/RJ_sin_data_pure_injection.csv",delimiter=',')


plt.plot(injected_data)
#plt.show()
plt.savefig("plots/RJ_sin_injected_data_raw.pdf")
plt.close()

plt.plot(injected_data)
plt.plot(injected_data_clean)
#plt.show()
plt.savefig("plots/RJ_sin_injected_data.pdf")
plt.close()

#data = gmcmc.trim_thin_file("data/chains_RJ_sin.hdf5",trim=None,ac=None)
data = np.loadtxt("data/chains_RJ_sin.hdf5",delimiter=',')
status =data[1000::100,max_dim:]
data =data[1000::100,:max_dim]
symm_mask = (status[:,-1] ==0)
print("Sym: ",np.sum(symm_mask))
nonsymm_mask = (status[:,-1] ==1)
print("NonSym: ",np.sum(nonsymm_mask))
data_sym = data[symm_mask,:-1]
data_nonsym = data[nonsymm_mask]
print(np.shape(data))
print(np.shape(data_sym))
print(np.shape(data_nonsym))
print(np.shape(status))
print(np.shape(status))
print("Length: ",len(data))
#for x in range(len(data[0])):
#    plt.plot(data[:,x])
#    plt.show()
#    plt.close()
xpts = np.arange(0,len(injected_data))


#labels = ["A",r"$\phi_0$","f"]
#labels = ["A",r"$\phi_0$","f",r"$\sigma$",r"$t_0$"]
labels = ["A",r"$\phi_0$","f",r"$\sigma$",r"$t_0$",r"$\alpha$"]
#fig,ax = plt.subplots(nrows=5,ncols=5,figsize=[8,8])
fig,ax = plt.subplots(nrows=6,ncols=6,figsize=[8,8])
fig = corner.corner(data_nonsym,labels=labels,quantiles=[.16,.84],bins=30,show_titles=True,fig=fig)
dim = len(injections)
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
plt.savefig("plots/mcmc_RJ_sin_nonsym.pdf")
plt.close()
#exit()

#labels = ["A",r"$\phi_0$","f"]
#labels = ["A",r"$\phi_0$","f",r"$\sigma$",r"$t_0$"]
labels = ["A",r"$\phi_0$","f",r"$\sigma$",r"$t_0$"]
#fig,ax = plt.subplots(nrows=5,ncols=5,figsize=[8,8])
fig,ax = plt.subplots(nrows=5,ncols=5,figsize=[8,8])
fig = corner.corner(data_sym,labels=labels,quantiles=[.16,.84],bins=30,show_titles=True,fig=fig)
dim = len(injections[:-1])
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
plt.savefig("plots/mcmc_RJ_sin_sym.pdf")
plt.close()


ct = 0
plt.plot(injected_data)
for x in data_nonsym[::100]:
    #plt.plot(x[0]*np.sin(x[2]*xpts + x[1]),alpha=.1,color='red')
    plt.plot(fn(x,xpts),alpha=.01,color='red')
    ct+=1
plt.savefig("plots/RJ_sin_reconstruction_full_data_nonsym.pdf")
plt.close()

ct = 0
plt.plot(injected_data)
for x in data_sym[::100]:
    #plt.plot(x[0]*np.sin(x[2]*xpts + x[1]),alpha=.1,color='red')
    plt.plot(fn_sym(x,xpts),alpha=.01,color='red')
    ct+=1
plt.savefig("plots/RJ_sin_reconstruction_full_data_sym.pdf")
plt.close()

ct = 0
for x in data_nonsym[::10]:
    #plt.plot(x[0]*np.sin(x[2]*xpts + x[1]),alpha=.1,color='red')
    plt.plot(fn(x,xpts),alpha=.01,color='red')
    #plt.plot(\
    #    x[0]*np.sin(x[2]*xpts + x[1])\
    #    *1./np.sqrt(2.*x[3]*x[3]) \
    #    * np.exp(-.5 * (x[4]-xpts)*(x[4]-xpts)/(x[3]*x[3]))\
    #    ,alpha=.01,color='red')
    ct+=1
plt.plot(injected_data_clean)
plt.savefig("plots/RJ_sin_reconstruction_clean_data_nonsym.pdf")
plt.close()

ct = 0
for x in data_sym[::10]:
    #plt.plot(x[0]*np.sin(x[2]*xpts + x[1]),alpha=.1,color='red')
    plt.plot(fn_sym(x,xpts),alpha=.01,color='red')
    #plt.plot(\
    #    x[0]*np.sin(x[2]*xpts + x[1])\
    #    *1./np.sqrt(2.*x[3]*x[3]) \
    #    * np.exp(-.5 * (x[4]-xpts)*(x[4]-xpts)/(x[3]*x[3]))\
    #    ,alpha=.01,color='red')
    ct+=1
plt.plot(injected_data_clean)
plt.savefig("plots/RJ_sin_reconstruction_clean_data_sym.pdf")
plt.close()

#xmedian = np.asarray([ np.percentile(y,50) for y in data.T])
#print(xmedian)
#fig,ax = plt.subplots(nrows=2,ncols=1)
#ax[0].plot(injected_data - fn(xmedian, xpts))
#ax[1].hist(injected_data - fn(xmedian, xpts), bins=np.linspace(-5,5,100))
#plt.savefig("plots/RJ_sin_residual.pdf")
#plt.close()



