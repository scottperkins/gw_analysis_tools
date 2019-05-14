import corner
import matplotlib.pyplot as plt
import numpy as np

#strain = np.loadtxt("data/gw150914_data.csv",delimiter=',',unpack=True)
#psd = np.loadtxt("data/gw150914_psd.csv",delimiter=',',unpack=True)
#freq = np.loadtxt("data/gw150914_freq.csv",delimiter=',',unpack=True)
#
#plt.loglog(freq,psd)
#plt.loglog(freq,strain[0]*strain[0]+strain[1]*strain[1])
##plt.show()
##plt.close()
#
#strain = np.loadtxt("data/gw_150914_data.csv",delimiter=',',unpack=True)
#psd = np.loadtxt("data/gw_150914_psd.csv",delimiter=',',unpack=True)
#freq = np.loadtxt("data/gw_150914_freq.csv",delimiter=',',unpack=True)
#
##plt.loglog(freq,psd)
##plt.loglog(freq,strain[0]*strain[0]+strain[1]*strain[1])
##plt.show()
##plt.close()
#
data = np.loadtxt("data/mcmc_output_7dim.csv",delimiter=',')
##data = data[:-100]
#dls = [x[0] for x in data]
#plt.plot(dls)
##plt.show()
#plt.close()
#
#tc = [x[1] for x in data]
#plt.plot(tc)
##plt.show()
#plt.close()
#
#phic = [x[2] for x in data]
#plt.plot(phic)
##plt.show()
#plt.close()
#
#chirpmasses = [x[3] for x in data]
#plt.plot(chirpmasses)
##plt.show()
#plt.close()
#
#etas = [x[4] for x in data]
#plt.plot(etas)
##plt.show()
#plt.close()
#
#chi1 = [x[5] for x in data]
#plt.plot(chi1)
##plt.show()
#plt.close()
#
#chi2 = [x[6] for x in data]
#plt.plot(chi2)
##plt.show()
#plt.close()
# #params = [dls,tc, phic,chirpmasses,etas,chi1,chi2]
#for x in np.arange(len(params)):
#    params[x] = np.asarray(params[x])/np.mean(params[x])
#    plt.plot(params[x],label=x)
#plt.legend()
#plt.show()
#plt.close()

ndim, nsamples = 7, len(data) 
labels = [r"$D_{L}$",r"$t_c$",r"$\phi_c$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]

figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_7dim.pdf")
plt.close()

autocorr = np.loadtxt("data/auto_corr_mcmc_7dim.csv",delimiter=',')
lengths = autocorr[0]
autocorr = autocorr[1:]
names = ["DL","tc","phic","CHIRP","ETA","CHI1","CHI2"]
for i in np.arange(len(autocorr)):
    plt.plot(lengths,autocorr[i], label=names[i])
plt.legend()
plt.savefig("autocorr_testing_7dim.pdf")
plt.close()
