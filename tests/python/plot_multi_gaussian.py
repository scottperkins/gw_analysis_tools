import numpy as np
import gwatpy.mcmc_routines as gmcmc
import corner
import matplotlib.pyplot as plt

data = gmcmc.trim_thin_file("data/gaussian_output_0.hdf5",ac=None,trim=None,recalc_ac=False)
print(np.std(data[:,0]))
lab=["w","x","y","z"]
fig = corner.corner(data,show_titles=True, labels=lab)
plt.savefig("plots/gaussian_mcmc.pdf")
plt.close()
