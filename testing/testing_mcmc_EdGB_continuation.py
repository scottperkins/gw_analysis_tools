import matplotlib as mpl
#mpl.use("pdf")
import corner
import matplotlib.pyplot as plt
import numpy as np
from phenompy.utilities import calculate_mass1, calculate_mass2, mpc
burnin = False
burninlength = 50000
datasets = []
numsets = 4
for i in np.arange(numsets):
    datasets.append( np.loadtxt("data/mcmc_output_EdGB{}.csv".format(i+1),delimiter=','))
data = []
for x in datasets:
    for y in x:
        data.append(y)
if burnin:
    data = data[burninlength:]
for i in np.arange(9):
    parameter = [x[i] for x in data]
    plt.plot(parameter)
    if i == 8:
        plt.yscale("log")
    plt.show()
    plt.close()
ndim, nsamples = 9, len(data) 
#labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
labels = [r"$cos\iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$",r"$\sqrt{\alpha}$"]
dataplot = []
for x in data:
    dataplot.append(x)
    dataplot[-1][-1]=(dataplot[-1][-1])**(1./4.)*(3e5)
dataplot = np.asarray(dataplot)
figure = corner.corner(dataplot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_EdGB.pdf")
plt.close()
alphahist = []
i = 0
for x in dataplot:
    if x[-1]>10:
        i+=1
    alphahist.append(x[-1])
alphahist = np.asarray(alphahist)
print("Minimum root alpha: ",alphahist.min())
print(i)
plt.hist(alphahist,bins=100,density=True)
plt.savefig("alpha_hist_EdGB.pdf")
plt.close()
##############################################################
datasets = []
for i in np.arange(numsets):
    datasets.append( np.loadtxt("data/mcmc_output_EdGB_hot{}.csv".format(i+1),delimiter=','))
data = []
for x in datasets:
    for y in x:
        data.append(y)
if burnin:
    data = data[burninlength:]
#data = data[:30000]
#data = data[:-100]
#chirpmasses = [x[1] for x in data]
#plt.plot(chirpmasses)
#plt.show()
#plt.close()
ndim, nsamples = 9, len(data) 
#labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]

dataplot = []
for x in data:
    dataplot.append(x)
    dataplot[-1][-1]=(dataplot[-1][-1])**(1./4.)*(3e5)
figure = corner.corner(dataplot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_EdGB_hot.pdf")
plt.close()
#
###############################################################
#autocorr = np.loadtxt("data/auto_corr_mcmc_EdGB.csv",delimiter=',')
#lengths = autocorr[0]
#autocorr = autocorr[1:]
#
#for i in np.arange(len(autocorr)):
#    plt.plot(lengths,autocorr[i], label=labels[i])
#plt.legend()
#plt.savefig("autocorr_testing_EdGB.pdf")
#plt.close()
#
###############################################################
#data = np.loadtxt("data/mcmc_output_EdGB.csv",delimiter=',')
#if burnin:
#    data = data[burninlength:]
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
#plt.savefig("mcmc_testing_transform_EdGB.pdf")
#plt.close()
