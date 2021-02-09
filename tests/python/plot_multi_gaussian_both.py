import numpy as np 
import gwatpy.mcmc_routines as gmcmc
import matplotlib.pyplot as plt
from corner import corner

gdat = gmcmc.trim_thin_file("data/gaussian_output_0_.hdf5",trim=None,ac=None)
edat = np.loadtxt("data/emcee_samples_multi_gaussian.csv",delimiter=',')
labels = np.arange(len(edat[0]))
fig = corner(edat, labels=labels,show_titles=True)
fig = corner(gdat,fig=fig, weights=np.ones(len(gdat))*len(edat)/len(gdat),color='blue')
plt.savefig("plots/gaussian_combined.pdf")
plt.close()
