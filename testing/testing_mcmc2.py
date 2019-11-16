import numpy as np
import corner
import matplotlib.pyplot as plt

data = np.genfromtxt("data/mcmc_pv2_chain.csv",delimiter=",")
print(np.shape(data))
labels = [r"RA",r"DEC","psi","phiref","iota",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$","chip","phip"]
figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.show()
plt.close()
