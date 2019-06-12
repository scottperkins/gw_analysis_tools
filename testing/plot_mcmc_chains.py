import matplotlib as mpl
mpl.use("pdf")
import corner
import matplotlib.pyplot as plt
import numpy as np

from phenompy.utilities import calculate_mass1, calculate_mass2, mpc
import sys

data_dir = "data/"
plot_dir = "plots/"
data_base_file = "mcmc_output_"
autocorr_base_file = "auto_corr_"
extension = ".csv"
event = sys.argv[1]
theory = sys.argv[2]
if theory =="dCS" or theory =="EdGB":
    labels = [r"$cos\iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$",r"$\sqrt{\alpha}$"]
if theory =="ppE_Inspiral" or theory =="ppE_IMR":
    labels = [r"$cos\iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$",r"$\beta$"]
if theory =="GR":
    labels = [r"$cos\iota$",r"RA",r"DEC",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
burn = False
burn_length =10000
data = np.loadtxt(data_dir+data_base_file+event+"_"+theory+extension,delimiter=',')
if burn:
    data = data[burn_length:]
#for i in np.arange(9):
#    parameter = [x[i] for x in data]
#    plt.plot(parameter)
#    if i == 8:
#        plt.yscale("log")
#    plt.show()
#    plt.close()
#ndim, nsamples = 9, len(data) 

dataplot = []
for x in data:
    dataplot.append(x)
    if theory == "dCS" or theory=="EdGB":
        dataplot[-1][-1]=(dataplot[-1][-1])**(1./4.)*(3e5)
dataplot = np.asarray(dataplot)
figure = corner.corner(dataplot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig(plot_dir+"corner_"+event+"_"+theory+".pdf")
plt.close()
if theory =="dCS" or theory=="EdGB":
    alphahist = []
    for x in dataplot:
        alphahist.append(x[-1])
    alphahist = np.asarray(alphahist)
    plt.hist(alphahist,bins=100,density=True)
    plt.savefig(plot_dir+"alpha_hist_"+event+"_"+theory+".pdf")
    plt.close()
##############################################################
data = np.loadtxt(data_dir+data_base_file+event+"_"+theory+"_hot"+extension,delimiter=',')
if burn:
    data = data[burn_length:]
#ndim, nsamples = 9, len(data) 
dataplot = []
for x in data:
    dataplot.append(x)
    if theory == "dCS" or theory=="EdGB":
        dataplot[-1][-1]=(dataplot[-1][-1])**(1./4.)*(3e5)
figure = corner.corner(dataplot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig(plot_dir+"corner_"+event+"_"+theory+"_hot.pdf")
plt.close()

##############################################################
autocorr = np.loadtxt(data_dir+autocorr_base_file+event+"_"+theory+extension,delimiter=',')
lengths = autocorr[0]
autocorr = autocorr[1:]

for i in np.arange(len(autocorr)):
    plt.plot(lengths,autocorr[i], label=labels[i])
plt.legend()
plt.savefig(plot_dir+"autocorr_"+event+theory+".pdf")
plt.close()

##############################################################
data = np.loadtxt(data_dir+data_base_file+event+"_"+theory+extension,delimiter=',')
if burn:
    data = data[burn_length:]
def chieff(m1, m2, spin1, spin2,):
    return (m1*spin1 + m2*spin2)/(m1+m2)
datatransform= []
for x in data:
    chirpm = x[4]
    symmratio = x[5]
    m1 = calculate_mass1(chirpm, symmratio)
    m2 = calculate_mass2(chirpm, symmratio)
    datatransform.append([m1,m2, chieff(m1,m2,x[6],x[7])])
datatransform = np.asarray(datatransform)

#ndim, nsamples = 3, len(datatransform) 
#labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
labels = [r"$M_1$",r"$M_2$",r"$\chi_{eff}$"]

figure = corner.corner(datatransform, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig(plot_dir+"corner_transformed_"+event+"_"+theory+".pdf")
plt.close()
