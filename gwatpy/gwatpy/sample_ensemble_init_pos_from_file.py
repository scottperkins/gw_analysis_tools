import numpy as np
import matplotlib.pyplot as plt
import gwatpy.mcmc_routines as gmcmc
import scipy.stats as sstats
from corner import corner

dim = 4
chains= 40
ensemble_chains= 10


old_data = gmcmc.trim_thin_file("data/test2.hdf5",ac=None,trim=None)

kernel = sstats.gaussian_kde(old_data.T)
new_data = kernel.resample(chains)
new_data = new_data.T
for x in new_data:
    if x[1] > .25:
        x[1] = .24
    if x[2] > 0.95 or x[2] < -0.95:
        x[2] = 0
    if x[3] > 0.95 or x[3] < -0.95:
        x[3] = 0
np.savetxt("data/sample_ensemble_init_pos.csv",new_data,delimiter=',',fmt='%.2f')

fig = corner(new_data)
plt.savefig("plots/ensemble_initial_pos.pdf")
plt.close()
