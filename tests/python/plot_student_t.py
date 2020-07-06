import matplotlib as mpl
#mpl.use("pdf")
import corner
import matplotlib.pyplot as plt
import numpy as np
import h5py
#import gwatpy.util as gpu
#import gwatpy.gwatpy_plot as gp; gp.set()
#from phenompy.utilities import calculate_mass1, calculate_mass2
import gwatpy.mcmc_routines as gmcmc

#data = np.loadtxt("data/mcmc_output.csv",delimiter=',')
#f = h5py.File("test_flat_standard.hdf5","r")
#data = f["THINNED_MCMC_OUTPUT"]["THINNED FLATTENED CHAINS"]
data = gmcmc.trim_thin_file("data/mcmc_output.csv",ac=17*5)

labels = ["X","Y"]
data_plot=[]
for x in np.arange(len(data)):
    if x% 1 == 0:
        data_plot.append(data[x])
data_plot = np.asarray(data_plot)
figure = corner.corner(data_plot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)

plt.savefig("plots/mcmc_student_t.pdf")
plt.close()
