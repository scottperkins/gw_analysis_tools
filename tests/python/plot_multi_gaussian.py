import numpy as np
import gwatpy.mcmc_routines as gmcmc
import corner
import matplotlib.pyplot as plt
import emcee

colors = ['b', 'g', 'r', 'c', 'm', 'y', 'm', 'y', 'y', 'm', 'y']

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


#############################################
#true_cov = np.loadtxt("data/multi_gaussian_like_cov.csv",delimiter=',')
#true_mean = np.loadtxt("data/multi_gaussian_like_mean.csv",delimiter=',')
#true_cov_prior = np.loadtxt("data/multi_gaussian_prior_cov.csv",delimiter=',')
#true_mean_prior = np.loadtxt("data/multi_gaussian_prior_mean.csv",delimiter=',')
#print(true_mean)
#print(true_cov)
#
#cov1 = np.linalg.inv((np.linalg.inv(true_cov_prior)+np.linalg.inv(true_cov)))
##mean1 = np.dot(true_mean,np.matmul(cov1,np.linalg.inv(true_cov))) +np.dot(true_mean_prior,np.matmul(cov1,np.linalg.inv(true_cov_prior))) 
#mean1 = np.dot(np.matmul(cov1,np.linalg.inv(true_cov)),true_mean) +np.dot(np.matmul(cov1,np.linalg.inv(true_cov_prior)),true_mean_prior) 
#print(mean1)
#
##data_true1 = np.random.multivariate_normal(true_mean,true_cov,size=int(1e4))
#data_true1 = np.random.multivariate_normal(mean1,cov1,size=int(1e4))
#
#mean2 = np.dot(np.matmul(cov1,np.linalg.inv(true_cov)),-1*true_mean) +np.dot(np.matmul(cov1,np.linalg.inv(true_cov_prior)),true_mean_prior) 
#print(mean2)
#
##data_true2 = np.random.multivariate_normal(-1*true_mean,true_cov,size=int(1e4))
#data_true2 = np.random.multivariate_normal(mean2,cov1,size=int(1e4))
#
##data_true = np.insert(data_true1,0,data_true2,axis=0)
#data_true = np.random.multivariate_normal(mean1+mean2,cov1+cov1,size=int(1e4))
#############################################

data_true = np.loadtxt("data/bilby/multivariate_gaussian_prior_samples.txt",skiprows=1)
#data_true = np.loadtxt("data/bilby/emcee_multivariate_gaussian_prior/chain.dat",skiprows=1)
#data_true = data_true[data_true[:,0] == 0,1:-2]
#print("True samples: ",len(data_true))
#ct = 0
#for x in data_true.T:
#    plt.plot(AC(x))
#    plt.savefig("plots/t{}.pdf".format(ct))
#    plt.close()
#    ct+=1

lab=np.arange(10)
data = gmcmc.trim_thin_file("data/gaussian_output_0.hdf5",ac=None,trim=None,recalc_ac=False)
#plt.hist(abs(data[:,0]),bins=50)
#plt.savefig("plots/temp.pdf")
#plt.close()
dim = len(data[0])
print("Samples: ",np.std(data[:,0]))
L1 = len(data)
print(len(data))

ranges = [(-10,10) for x in np.arange(dim)]
#ranges = [(0,10),(-4,4),(-1,8),(-4,4),(-4,8),(0,8)]

fig = corner.corner(data,show_titles=True, labels=lab,bins=50, weights=np.ones(len(data))/len(data),color= colors[0],range=ranges)


fig = corner.corner(data_true,fig=fig,show_titles=True, labels=lab,bins=50, weights=np.ones(len(data_true))/len(data_true),color= colors[1],range=ranges)

plt.savefig("plots/true_gaussian_mcmc.pdf")
plt.close()
exit()

data_full = np.copy(data)

for l in np.arange(9):
#for l in [2]:
    data = gmcmc.trim_thin_file("data/gaussian_output_{}.hdf5".format(l+1),ac=None,trim=None,recalc_ac=False)
    data_full = np.insert(data_full, -1, np.copy(data),axis=0)
    dim = len(data[0])
    L1 = len(data)
    print("Samples: ",L1)
    fig = corner.corner(data,show_titles=True, fig = fig,labels=lab, bins=50,weights=np.ones(len(data))/len(data),color= colors[l+1],range=ranges)
plt.savefig("plots/gaussian_mcmc.pdf")
plt.close()

for x in np.arange(len(data_full[0])):
    print(emcee.autocorr.integrated_time(data_full[:,x],tol=0))
    plt.plot(data_full[:,x])
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
