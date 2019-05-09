import corner
import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("data/mcmc_output.csv",delimiter=',')
ndim, nsamples = 5, len(data) 
labels = [r"$D_{L}$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]

figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing.pdf")
plt.close()

autocorr = np.loadtxt("data/auto_corr_mcmc.csv",delimiter=',')
lengths = autocorr[0]
autocorr = autocorr[1:]
for i in autocorr:
    plt.plot(lengths,i)
plt.savefig("autocorr_testing.pdf")
plt.close()
