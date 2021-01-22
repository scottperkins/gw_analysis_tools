import gwatpy.mcmc_routines as gmcmc
import numpy as np
import matplotlib.pyplot as plt

#fig = gmcmc.plot_convergence("data/continue_output_injection.hdf5",ac=None,trim=None)
fig = gmcmc.plot_convergence("data/test.hdf5",ac=None,trim=None)
ct = 0
for f in fig:
    plt.figure(f.number)
    plt.savefig("plots/convergence_testing_{}.pdf".format(ct))
    ct+=1
plt.close()
