import gwatpy.mcmc_routines as gmcmc
import matplotlib.pyplot as plt
import numpy as np

filename = "data/light.hdf5"
#filename = "data/heavy.hdf5"
psd_file = "/home/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt"
#filename = "data/test.hdf5"
#psd_file = "/home/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt"

generation_method_base = "PNSeries_ppE_IMRPhenomD_Inspiral"
#generation_method_base = "ppE_IMRPhenomD_Inspiral"
#generation_method_base = "IMRPhenomD"
mod_kwargs = {"ppE_Nmod":6,"bppe":[-7,-5,-4,-3,-2,-1]}
#mod_kwargs = {"ppE_Nmod":6,"bppe":[0,0,0,0,0,0]}
#mod_kwargs = {}
threads = 10
inject = np.loadtxt("data/injections.csv",dtype=np.float64)
inject_stat = np.ones(len(inject))
#inject = None
#inject_stat = None
xlim = [2.8,3.1]
#xlim = None

fig = gmcmc.plot_bayesogram(filename, psd_file, "Hanford", generation_method_base, psd_column=0,threads=threads, xlim = xlim,injection=inject, injection_status=inject_stat,gmst=2.45682,mod_struct_kwargs=mod_kwargs)
plt.savefig("plots/bayesogram.pdf")
plt.close()

fig = gmcmc.plot_injection(inject, inject_stat, psd_file, "Hanford", generation_method_base, psd_column=0,threads=threads, xlim = xlim,gmst=2.45682,mod_struct_kwargs=mod_kwargs)
plt.savefig("plots/injection.pdf")
plt.close()
