import numpy as np
import matplotlib.pyplot as plt
import corner 

data = np.loadtxt("data/mcmc_output.csv",delimiter=',',unpack=True)
burnin =1000
#data = data[burnin:]
plt.plot(data[0][burnin:])
plt.plot(data[1][burnin:])
plt.show()
plt.close()
#plt.hist(data,bins=100,density=True)
#x = np.linspace(-3,3)
data = np.loadtxt("data/mcmc_output.csv",delimiter=',')
ndim, nsamples = 2, len(data) 
labels = [r"x",r"y"]

figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("neil_mcmc_testing.pdf")
plt.close()


autocorr = np.loadtxt("data/auto_corr_mcmc.csv",delimiter=',')
lengths = autocorr[0]
autocorr = autocorr[1:]
for i in autocorr:
    plt.plot(lengths,i)
plt.savefig("neil_autocorr_testing.pdf")
plt.close()
