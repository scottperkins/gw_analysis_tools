import matplotlib as mpl
#mpl.use("pdf")
import corner
import matplotlib.pyplot as plt
import numpy as np
import h5py
import gwatpy.mcmc_routines as gmcmc
labels = ["X","Y"]
data = gmcmc.trim_thin_file("data/mcmc_output.hdf5",trim=None,ac=None)
ct=0
for x in data.T:
    plt.plot(x)
    plt.savefig("plots/trace_rosen_{}.pdf".format(ct))
    plt.close()
    ct+=1
data_plot = np.asarray(data)
figure = corner.corner(data_plot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)

plt.savefig("plots/mcmc_student_t.pdf")
plt.close()
