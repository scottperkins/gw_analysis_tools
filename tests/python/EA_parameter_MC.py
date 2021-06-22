import numpy as np
import matplotlib.pyplot as plt

iterations = 13
n_bins = 50
for i in np.arange(iterations):
    data = np.loadtxt("data/EA_parameter_MC.csv".format(i),delimiter=',',unpack=True)
    #print(data[i])

#print(data[0])
plt.hist(data[0], n_bins)
plt.title(r'$c_a$')
plt.savefig("plots/hist_ca")
plt.close()

plt.hist(data[1], n_bins)
plt.title(r'$c_\theta$')
plt.savefig("plots/hist_ctheta")
plt.close()

plt.hist(data[2], n_bins)
plt.title(r'$c_\omega$')
plt.savefig("plots/hist_cw")
plt.close()

plt.hist(data[3], n_bins)
plt.title(r'$c_\sigma$')
plt.savefig("plots/hist_csigma")
plt.close()

bins = np.linspace(.5,1.5,n_bins)
#print(bins)
plt.hist(data[4], bins=bins)
plt.title(r'$c_T$')
plt.savefig("plots/hist_cT")
plt.close()

plt.hist(data[5], n_bins)
plt.title(r'$c_S$')
plt.savefig("plots/hist_cS")
plt.close()

plt.hist(data[6], n_bins)
plt.title(r'$c_V$')
plt.savefig("plots/hist_cV")
plt.close()

plt.hist(data[7], n_bins)
plt.title(r'$\alpha_{ppE,0}^{(2,T)}$')
plt.savefig("plots/hist_alphaT")
plt.close()

plt.hist(data[8], n_bins)
plt.title(r'$g_{b1}$')
plt.savefig("plots/hist_gb1")
plt.close()

plt.hist(data[9], n_bins)
plt.title(r'$a_{bL}$')
plt.savefig("plots/hist_abL")
plt.close()

plt.hist(data[10], n_bins)
plt.title(r'$g_{X1}$')
plt.savefig("plots/hist_gX1")
plt.close()
