import numpy as np
import gwatpy.mcmc_routines as gmcmc
import corner
import matplotlib.pyplot as plt

lab=["w","x","y","z"]

data = gmcmc.trim_thin_file("data/gaussian_output_0.hdf5",ac=None,trim=None,recalc_ac=False)
print(np.std(data[:,0]))
fig = corner.corner(data,show_titles=True, labels=lab)
#plt.savefig("plots/gaussian_mcmc_0.pdf")
#plt.close()

data = gmcmc.trim_thin_file("data/gaussian_output_19.hdf5",ac=None,trim=None,recalc_ac=False)
print(np.std(data[:,0]))
corner.corner(data,fig=fig,show_titles=True, labels=lab)
plt.savefig("plots/gaussian_mcmc_comb_small_E.pdf")
plt.close()

#lab=["w","x","y","z"]
#fig = corner.corner(data,show_titles=True, labels=lab)
#plt.savefig("plots/gaussian_mcmc_19.pdf")
#plt.close()

data = gmcmc.trim_thin_file("data/gaussian_output_0.hdf5",ac=None,trim=None,recalc_ac=False)
for x in np.arange(19):
    t = gmcmc.trim_thin_file("data/gaussian_output_{}.hdf5".format(x+1),ac=None,trim=None,recalc_ac=False)
    data = np.insert(data, [-1],t, axis=0)
fig = corner.corner(data,show_titles=True, labels=lab)
plt.savefig("plots/gaussian_mcmc_full_small_E.pdf")
plt.close()
