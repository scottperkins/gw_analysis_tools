import matplotlib as mpl
#mpl.use("pdf")
import corner
import matplotlib.pyplot as plt
import numpy as np
import h5py
import gwatpy.mcmc_routines as gmcmc
labels = ["X"]
data = gmcmc.trim_thin_file("data/mcmc_output_1d.hdf5",trim=None,ac=None)
data_plot = np.asarray(data)
figure = corner.corner(data_plot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)

plt.savefig("plots/mcmc_1d.pdf")
plt.close()
