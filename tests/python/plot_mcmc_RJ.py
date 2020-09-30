import matplotlib as mpl
import corner
import matplotlib.pyplot as plt
import numpy as np
import gwatpy.mcmc_routines as gmcmc
import h5py
import scipy
import emcee
#import gwatpy.util as gpu
#import gwatpy.gwatpy_plot as gp; gp.set()
#from phenompy.utilities import calculate_mass1, calculate_mass2

#data = np.loadtxt("data/injection_output.csv",delimiter=',')
#data = np.loadtxt("data/experiment_output.csv",delimiter=',')
#data = np.loadtxt("data/test_output.csv",delimiter=',')
injections = np.loadtxt("data/injections.csv",delimiter=',',unpack=True)
f = h5py.File("data/RJ_injection_output.csv",'r')


data,status,model = gmcmc.RJPTMCMC_unpack_file("data/RJ_injection_output.csv")
for x in data:
    x[6] = np.exp(x[6])
    x[7] = np.exp(x[7])
GRdim = 11
mask_GR = np.sum(status[:,GRdim:],axis=1) == 0
mask_nonGR = np.logical_not( mask_GR)
data_GR = data[mask_GR]
data_GR = data_GR[:,:GRdim]
data_nonGR = data[mask_nonGR]
data_nonGR = data_nonGR[:]
print(np.shape(data_GR))
print(np.shape(data_nonGR))
print("non Isolated")
print(np.sum((status[:,-1] == 1) ) )
print(np.sum((status[:,-2] == 1) ) )
print(np.sum((status[:,-3] == 1) ) )
print("Isolated")
print(np.sum((status[:,-1] == 1) & (status[:,-2] == 0)& (status[:,-3] == 0) ))
print(np.sum((status[:,-1] == 0) & (status[:,-2] == 1)& (status[:,-3] == 0) ))
print(np.sum((status[:,-1] == 0) & (status[:,-2] == 0)& (status[:,-3] == 1) ))

#for x in np.arange(len(data_GR[0])):
#    plt.plot(data_GR[:,x])
#    plt.show()
#    plt.close()

numMods =  np.sum(status[:,GRdim:],axis=1)
print(numMods)
plt.hist(numMods)
plt.savefig('plots/mod1.pdf')
#plt.show()
plt.close()

plt.hist(data[:,-1],bins=50,density=True)
plt.savefig('plots/mod2.pdf')
#plt.show()
plt.close()

plt.hist(data[:,-2],bins=50,density=True)
plt.savefig('plots/mod3.pdf')
#plt.show()
plt.close()

plt.hist(data[:,-3],bins=50,density=True)
plt.savefig('plots/mod4.pdf')
#plt.show()
plt.close()

#for x in range(len(data_GR[0])):
#    print(emcee.autocorr.integrated_time(data_GR[:,x]))

#chains = list(f["MCMC_OUTPUT/LOGL_LOGP"])
#for y in range(len(list(f["MCMC_OUTPUT/LOGL_LOGP"]))):
#    print(chains[y])
#    lllpdata = f["MCMC_OUTPUT/LOGL_LOGP"][chains[y]]
#    def _line(x,m,b):
#        return x*m +b
#    x = np.arange(len(lllpdata[:,0]))
#    popt,pcov = scipy.optimize.curve_fit(_line,x,lllpdata[:,0])
#    print(popt)
#    plt.plot(lllpdata[:,0],label="{} \ {}".format(y,f["MCMC_METADATA"]["CHAIN TEMPERATURES"][y]))
#    plt.plot(x,_line(x,*popt))
#    plt.plot(lllpdata[:,1],label='P')
#    plt.plot(lllpdata[:,1]+lllpdata[:,0],label='P+L')
#    plt.legend()
#    plt.show()
#    plt.close()
#    for x in np.arange(len(f["MCMC_OUTPUT/"+chains[y]][0])):
#        plt.plot(f["MCMC_OUTPUT/"+chains[y]][:,x],label=chains[y]+"-"+str(x))
#        #plt.plot(testdata2[:,x])
#        plt.legend()
#        plt.show()
#        plt.close()
#testdata=f["MCMC_OUTPUT"]["CHAIN 0"]
#for x in np.arange(len(testdata[0])):
#    plt.plot(testdata[:,x])
#    plt.show()
#    plt.close()
#exit()
injections[6] = np.exp(injections[6])
injections[7] = np.exp(injections[7])
ndim, nsamples = 15, len(data) 
labels = [r"$\alpha$",r"$\sin(\delta)$",r"$\psi$",r"$\cos(\iota)$","$\phi_{ref}$","$t_c$",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$a_{1}$",r"$a_2$",r"$\cos \theta_1$",r"$\cos \theta_2$",r"$\phi_1$",r"$\phi_2$"]
figure = corner.corner(data_GR, labels=labels,quantiles=[.1,.5,.9], show_titles=True)
dimGR = len(data_GR[0])
axes = np.array(figure.axes).reshape(dimGR,dimGR)
for i in np.arange(dimGR):
    ax = axes[i,i]
    ax.axvline(injections[i])

for yi in np.arange(dimGR):
    for xi in np.arange(yi):
        ax = axes[yi,xi]
        ax.axvline(injections[xi])
        ax.axhline(injections[yi])
        ax.plot(injections[xi],injections[yi])

plt.savefig("plots/mcmc_injection_RJ.pdf")
#plt.savefig("plots/mcmc_experiment.pdf")
plt.close()

#labels = [r"$\alpha$",r"$\sin(\delta)$",r"$\psi$",r"$\cos(\iota)$","$\phi_{ref}$","$t_c$",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$a_{1}$",r"$a_2$",r"$\cos \theta_1$",r"$\cos \theta_2$",r"$\phi_1$",r"$\phi_2$",r'$\delta \phi$']
#figure = corner.corner(data_nonGR, labels=labels,quantiles=[.1,.5,.9], show_titles=True)
#dimnonGR = len(data_nonGR[0])
#axes = np.array(figure.axes).reshape(dimnonGR,dimnonGR)
#for i in np.arange(dimnonGR):
#    ax = axes[i,i]
#    ax.axvline(injections[i])
#
#for yi in np.arange(dimnonGR):
#    for xi in np.arange(yi):
#        ax = axes[yi,xi]
#        ax.axvline(injections[xi])
#        ax.axhline(injections[yi])
#        ax.plot(injections[xi],injections[yi])
#
#plt.savefig("plots/mcmc_injection_RJ_nonGR.pdf")
##plt.savefig("plots/mcmc_experiment.pdf")
#plt.close()
