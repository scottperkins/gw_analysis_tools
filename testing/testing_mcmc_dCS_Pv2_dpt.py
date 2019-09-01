import corner
import matplotlib.pyplot as plt
import numpy as np
from phenompy.utilities import calculate_mass1, calculate_mass2

burn = False
burnlength = 80000
num_files = 5
data =[] 
for i in np.arange(num_files):
    data.append(np.loadtxt("data/mcmc_output_dCS_Pv2_{}.csv".format(i),delimiter=','))
data_unpack = []
for x in data:
    for y in x:
        data_unpack.append(y)
if burn:
    data_unpack = data_unpack[burnlength:]
labels = [r"$cos \iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$",r"$\theta_1$", r"$\theta_2$",r"$\phi_1$", r"$\phi_2$",r"$\phi_{ref}$", r"$\psi$",r"$\sqrt{\alpha}$"]

dataplot=[]
for x in data_unpack:
    dataplot.append(x)
    dataplot[-1][3] = np.exp(dataplot[-1][3])
    dataplot[-1][4] = np.exp(dataplot[-1][4])
figure = corner.corner(dataplot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_dCS_Pv2.pdf")
plt.close()

##############################################################
data =[] 
for i in np.arange(num_files):
    data.append(np.loadtxt("data/mcmc_output_dCS_Pv2_{}.csv".format(i),delimiter=','))
data_unpack = []
for x in data:
    for y in x:
        data_unpack.append(y)
if burn:
    data_unpack = data_unpack[burnlength:]
labelstrans = [r"$cos \iota$",r"RA",r"DEC",r"$D_L$",r"$M_1$",r"$M_2$",r"$\chi_{1}$",r"$\chi_{2}$",r"$cos(\theta_1)$",r"$cos(\theta_2)$",r"$\sqrt{\alpha}$"]
#
dataplot=[]
for x in data_unpack:
    dataplot.append(np.ones(11))
    dataplot[-1][0] = x[0]
    dataplot[-1][1] = x[1]
    dataplot[-1][2] = x[2]
    dataplot[-1][3] = np.exp(x[3])
    dataplot[-1][4] = calculate_mass1(np.exp(x[4]),x[5])
    dataplot[-1][5] = calculate_mass2(np.exp(x[4]),x[5])
    dataplot[-1][6] = x[6]
    dataplot[-1][7] = x[7]
    dataplot[-1][8]= np.cos(x[8])
    dataplot[-1][9]= np.cos(x[9])
    dataplot[-1][10]= x[14]
dataplot = np.asarray(dataplot)
figure = corner.corner(dataplot, labels=labelstrans,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_dCS_Pv2_transformed.pdf")
plt.close()
##############################################################
#data = np.loadtxt("data/mcmc_output_ppE_Pv2_hot.csv",delimiter=',')
#ndim, nsamples = 7, len(data) 
#figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
#plt.savefig("mcmc_testing_ppE_Pv2_hot.pdf")
#plt.close()

##############################################################
#autocorr = np.loadtxt("data/auto_corr_mcmc_dCS_Pv2.csv",delimiter=',')
#lengths = autocorr[0]
#autocorr = autocorr[1:]
#for i in np.arange(len(autocorr)):
#    plt.plot(lengths,autocorr[i], label=labels[i])
#plt.legend()
#plt.savefig("autocorr_testing_dCS_Pv2.pdf")
#plt.close()

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
