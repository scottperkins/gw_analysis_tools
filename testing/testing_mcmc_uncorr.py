import matplotlib as mpl
#mpl.use("pdf")
import corner
import matplotlib.pyplot as plt
import numpy as np
from phenompy.utilities import calculate_mass1, calculate_mass2

data = np.loadtxt("data/mcmc_output_uncorr_P.csv",delimiter=',')
#data = data[2000:]
data_thinned = []
for x in np.arange(len(data)):
    if x%1 ==0:
        data_thinned.append(data[x]) 
#for i in np.arange(8):
#    parameter = [x[i] for x in data]
#    plt.plot(parameter)
#    plt.show()
#    plt.close()
ndim, nsamples = 11, len(data) 
labels = [r"$\alpha$",r"$\sin(\delta)$",r"$\psi$",r"$\cos(\iota)$","$\phi_{ref}$","$t_c$",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$",r"$\chi_p$",r"$\phi_p$"]
#labels = [r"$\alpha$",r"$\sin(\delta)$",r"$\psi$",r"$\cos(\iota)$","$\phi_{ref}$","$t_c$",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
for x in data:
    x[6] = np.exp(x[6])
    x[7] = np.exp(x[7])

figure = corner.corner(data_thinned, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
#plt.savefig("mcmc_testing_D_uncorr.pdf")
plt.savefig("mcmc_testing_P_uncorr.pdf")
plt.close()
##############################################################
#data = np.loadtxt("data/mcmc_output_DFull_hot.csv",delimiter=',')
##data = data[:-100]
##chirpmasses = [x[1] for x in data]
##plt.plot(chirpmasses)
##plt.show()
##plt.close()
#ndim, nsamples = 8, len(data) 
##labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
#labels = [r"$sos\iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
#
#figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
#plt.savefig("mcmc_testing_DFull_hot.pdf")
#plt.close()
#
###############################################################
#autocorr = np.loadtxt("data/auto_corr_mcmc_DFull.csv",delimiter=',')
#lengths = autocorr[0]
#autocorr = autocorr[1:]
#
#labels = [r"$cos\iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
#for i in np.arange(len(autocorr)):
#    plt.plot(lengths,autocorr[i], label=labels[i])
#plt.legend()
#plt.savefig("autocorr_testing_DFull.pdf")
#plt.close()
#
###############################################################
#data = np.loadtxt("data/mcmc_output_DFull.csv",delimiter=',')
#def chieff(m1, m2, spin1, spin2,):
#    return (m1*spin1 + m2*spin2)/(m1+m2)
#datatransform= []
#for x in data:
#    chirpm = x[4]
#    symmratio = x[5]
#    m1 = calculate_mass1(chirpm, symmratio)
#    m2 = calculate_mass2(chirpm, symmratio)
#    datatransform.append([m1,m2, chieff(m1,m2,x[6],x[7])])
#datatransform = np.asarray(datatransform)
#
#ndim, nsamples = 3, len(datatransform) 
##labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
#labels = [r"$M_1$",r"$M_2$",r"$\chi_{eff}$"]
#
#figure = corner.corner(datatransform, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
#plt.savefig("mcmc_testing_transform_DFull.pdf")
#plt.close()