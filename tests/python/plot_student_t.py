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
f = h5py.File("data/mcmc_output.csv","r")
#temps = f["MCMC_METADATA"]["CHAIN TEMPS"]
temps =f["MCMC_METADATA"]["CHAIN TEMPERATURES"][:] 
chains =[]
chains_cold =[]
LL =[]
for x in np.arange(len(f["MCMC_OUTPUT"])-1):
    chains.append(f["MCMC_OUTPUT"]["CHAIN "+str(int(x))][:])
    LL.append(f["MCMC_OUTPUT"]["LOGL_LOGP"]["CHAIN "+str(int(x))][:])
    if(temps[x] ==1):
        chains_cold.append(f["MCMC_OUTPUT"]["CHAIN "+str(int(x))][:])
chains=np.asarray(chains)
LL=np.asarray(LL)
print(np.shape(LL))
#for x in np.arange(len(chains)):
#    if(temps[x] == 1e14):
#        plt.plot(chains[x][:,0])
#        plt.show()
#        plt.close()
#exit()
#for x in np.arange(len(LL)):
#    if(temps[x] == 1e14):
#        plt.plot(LL[x])
#        plt.show()
#exit()

flat = chains_cold[0]
for x in chains_cold[1:]:
    flat = np.concatenate((flat, x))
labels = ["X","Y"]
figure = corner.corner(flat[::20], labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("plots/mcmc_student_t.pdf")
plt.close()

exit()









data = gmcmc.trim_thin_file("data/mcmc_output.csv",trim=None,ac=18)
print(np.shape(data))
labels = ["X","Y"]
data_plot=[]
for x in np.arange(len(data)):
    if x% 1 == 0:
        data_plot.append(data[x])
data_plot = np.asarray(data_plot)
figure = corner.corner(data_plot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)

plt.savefig("plots/mcmc_student_t.pdf")
plt.close()
