import corner
import matplotlib.pyplot as plt
import numpy as np
from phenompy.utilities import calculate_mass1, calculate_mass2

data = np.loadtxt("data/mcmc_output_Pv2.csv",delimiter=',')
#data = data[:-100]
#chirpmasses = [x[1] for x in data]
#plt.plot(chirpmasses)
#plt.show()
#plt.close()
ndim, nsamples = 13, len(data) 
#labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
labels = [r"$cos \iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$",r"$\theta_1$", r"$\theta_2$",r"$\phi_1$", r"$\phi_2$",r"$\phi_{ref}$"]

dataplot=[]
for x in data:
    dataplot.append(x)
    dataplot[-1][3] = np.exp(dataplot[-1][3])
    dataplot[-1][4] = np.exp(dataplot[-1][4])
figure = corner.corner(dataplot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_Pv2.pdf")
plt.close()

##############################################################
#labels = [r"$cos \iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$",r"$\theta_1$", r"$\theta_2$",r"$\phi_1$", r"$\phi_2$",r"$\phi_{ref}$"]
#
#dataplot=[]
#for x in data:
#    dataplot.append(x)
#    dataplot[-1][3] = np.exp(dataplot[-1][3])
#    dataplot[-1][4] = np.exp(dataplot[-1][4])
#figure = corner.corner(dataplot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
#plt.savefig("mcmc_testing_Pv2.pdf")
#plt.close()
##############################################################
#data = np.loadtxt("data/mcmc_output_Pv2_hot.csv",delimiter=',')
#data = data[:-100]
#chirpmasses = [x[1] for x in data]
#plt.plot(chirpmasses)
#plt.show()
#plt.close()
#ndim, nsamples = 7, len(data) 
#labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
#labels = [r"$cos J_N$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$",r"$cos \theta_1$", r"$cos \theta_2$"]

#figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
#plt.savefig("mcmc_testing_Pv2_hot.pdf")
#plt.close()

##############################################################
autocorr = np.loadtxt("data/auto_corr_mcmc_Pv2.csv",delimiter=',')
lengths = autocorr[0]
autocorr = autocorr[1:]
for i in np.arange(len(autocorr)):
    plt.plot(lengths,autocorr[i], label=labels[i])
plt.legend()
plt.savefig("autocorr_testing_Pv2.pdf")
plt.close()

##############################################################
#data = np.loadtxt("data/mcmc_output_Pv2.csv",delimiter=',')
#def chieff(m1, m2, spin1, spin2,):
#    return (m1*spin1 + m2*spin2)/(m1+m2)
#datatransform= []
#for x in data:
#    chirpm = x[1]
#    symmratio = x[2]
#    m1 = calculate_mass1(chirpm, symmratio)
#    m2 = calculate_mass2(chirpm, symmratio)
#    datatransform.append([m1,m2, chieff(m1,m2,x[3]*x[5],x[4]*x[6])])
#datatransform = np.asarray(datatransform)
#
#ndim, nsamples = 3, len(datatransform) 
##labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
#labels = [r"$M_1$",r"$M_2$",r"$\chi_{eff}$"]
#
#figure = corner.corner(datatransform, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
#plt.savefig("mcmc_testing_transform_Pv2.pdf")
#plt.close()
