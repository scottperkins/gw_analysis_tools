import corner
import matplotlib.pyplot as plt
import numpy as np
from phenompy.utilities import calculate_mass1, calculate_mass2

data = np.loadtxt("data/mcmc_output.csv",delimiter=',')
#data = data[:-100]
#chirpmasses = [x[1] for x in data]
#plt.plot(chirpmasses)
#plt.show()
#plt.close()
ndim, nsamples = 4, len(data) 
#labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
labels = [r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]

figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing.pdf")
plt.close()
##############################################################
data = np.loadtxt("data/mcmc_output_hot.csv",delimiter=',')
#data = data[:-100]
#chirpmasses = [x[1] for x in data]
#plt.plot(chirpmasses)
#plt.show()
#plt.close()
ndim, nsamples = 4, len(data) 
#labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
labels = [r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]

figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_hot.pdf")
plt.close()

##############################################################
autocorr = np.loadtxt("data/auto_corr_mcmc.csv",delimiter=',')
lengths = autocorr[0]
autocorr = autocorr[1:]
names = ["CHIRP","ETA","CHI1","CHI2"]
for i in np.arange(len(autocorr)):
    plt.plot(lengths,autocorr[i], label=names[i])
plt.legend()
plt.savefig("autocorr_testing.pdf")
plt.close()

##############################################################
data = np.loadtxt("data/mcmc_output.csv",delimiter=',')
def chieff(m1, m2, spin1, spin2,):
    return (m1*spin1 + m2*spin2)/(m1+m2)
datatransform= []
for x in data:
    datatransform.append([calculate_mass1(x[0],x[1]),calculate_mass2(x[0],x[1] ), chieff(calculate_mass1(x[0],x[1]),calculate_mass1(x[0],x[1]),x[2],x[3])])
datatransform = np.asarray(datatransform);

ndim, nsamples = 3, len(datatransform) 
#labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
labels = [r"$M_1$",r"$M_2$",r"$\chi_{eff}$"]

figure = corner.corner(datatransform, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_transform.pdf")
plt.close()
