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
##############################################################
autocorr = np.loadtxt(data_dir+autocorr_base_file+event+"_"+theory+extension,delimiter=',')
lengths = autocorr[0]
autocorr = autocorr[1:]

for i in np.arange(len(autocorr)):
    plt.plot(lengths,autocorr[i], label=labels[i])
plt.legend()
plt.savefig(plot_dir+"autocorr_"+event+"_"+theory+".pdf")
plt.close()

##############################################################
