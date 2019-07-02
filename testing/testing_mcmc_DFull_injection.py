import matplotlib as mpl
#mpl.use("pdf")
import corner
import matplotlib.pyplot as plt
import numpy as np
from phenompy.utilities import calculate_mass1, calculate_mass2
burnin = 0
datasets = []
numsets = 5
for i in np.arange(numsets):
    datasets.append( np.loadtxt("data/mcmc_output_injection{}.csv".format(i+1),delimiter=','))
data = []
for x in datasets:
    for y in x:
        data.append(y)
data = data[burnin:]
dataplot = []
for x in data:
    dataplot.append(x)
    dataplot[-1][-1] = (dataplot[-1][-1])**(1./4.)*(3e5) 
#data = []
#for i in np.arange(len(dataload)):
#    if (i%10==0):
#        data.append(dataload[i]) 
#data = data[1000:]
for i in np.arange(9):
    parameter = [x[i] for x in data]
    plt.plot(parameter)
    plt.show()
    plt.close()
ndim, nsamples = 9, len(data) 
#labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
#labels = [r"$cos\iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
labels = [r"$cos\iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$", r'\sqrt{\alpha}']

figure = corner.corner(dataplot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_injection.pdf")
plt.close()
##############################################################
#dataload = np.loadtxt("data/mcmc_output_injection_hot.csv",delimiter=',')
#data = dataload[burnin:]
#dataplot = []
#for x in data:
#    dataplot.append(x)
#    dataplot[-1][-1] = (dataplot[-1][-1])**(1./4.)*(3e5) 
##data = data[:-100]
##chirpmasses = [x[1] for x in data]
##plt.plot(chirpmasses)
##plt.show()
##plt.close()
#ndim, nsamples = 8, len(data) 
#
#figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
#plt.savefig("mcmc_testing_injection_hot.pdf")
#plt.close()
#
################################################################
for i in np.arange(numsets):
    autocorr = np.loadtxt("data/auto_corr_mcmc_injection{}.csv".format(i+1),delimiter=',')
    lengths = autocorr[0]
    autocorr = autocorr[1:]
    
    labels = [r"$cos\iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$", r'\sqrt{\alpha}']
    for i in np.arange(len(autocorr)):
        plt.plot(lengths,autocorr[i], label=labels[i])
    plt.legend()
    #plt.savefig("autocorr_testing_injection.pdf")
    #plt.show()
    plt.close()
#
###############################################################
#data = np.loadtxt("data/mcmc_output_injection.csv",delimiter=',')
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
#plt.savefig("mcmc_testing_transform_injection.pdf")
#plt.close()
