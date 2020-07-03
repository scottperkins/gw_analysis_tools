import numpy as np
import h5py
import matplotlib.pyplot as plt
import corner
import gwatpy.mcmc_routines as gmcmc

#data = np.loadtxt("data/mcmc_output_RB.csv",delimiter=',')
#f = h5py.File("data/mcmc_output_RB.hdf5",'r')
#ff = h5py.File("data/mcmc_output_RB_trimmed.hdf5",'r')
#print(f.keys())
#print(f["MCMC_METADATA"]["CHAIN TEMPERATURES"][:])
#print(f["MCMC_METADATA"]["SUGGESTED TRIM LENGTHS"][:])
#print(f["MCMC_METADATA"]["AC VALUES"][:])
#ff = h5py.File("test.hdf5",'r')
#datasuper = ff["MCMC_OUTPUT"]["CHAIN 0"]
#data_old = np.loadtxt("data/mcmc_output_RB.csv",delimiter=',',unpack=True)
#print(len(ff["MCMC_OUTPUT"]["CHAIN 0"]))
#exit()
#data = f["THINNED_MCMC_OUTPUT"]["THINNED FLATTENED CHAINS"]
#data = f["MCMC_OUTPUT"]["CHAIN 0"]

#f = h5py.File("data/mcmc_output_RB.hdf5",'r')
#print(f["MCMC_OUTPUT"].keys())
#dataLL = f["MCMC_OUTPUT"]["LOGL_LOGP"]
#for x in dataLL:
#    print(x)
#    dat = f["MCMC_OUTPUT"]["LOGL_LOGP"][x]
#    print(np.shape(dat))
#    plt.plot(dat[:,0],label="L")
#    plt.plot(dat[:,1],label="P")
#    plt.legend()
#    plt.show()
#    plt.close()

data = gmcmc.trim_thin_file("data/mcmc_output_RB.hdf5",trim=None,ac=None)
print(np.shape(data))
#exit()
#dataTrim = ff["MCMC_OUTPUT"]["CHAIN 0"]
dataT = np.transpose(data)
#dataTTrim = np.transpose(dataTrim)
for x in range(len(dataT)):
    plt.plot(dataT[x])
    #plt.plot(dataTTrim[x])
    #plt.plot(data_old[x])
    plt.show() 
    plt.close() 
fig = corner.corner(data)
plt.savefig("plots/mcmc_RB.pdf")
plt.close()
#fig = corner.corner(dataTrim)
#plt.savefig("plots/mcmc_RB_trim.pdf")
#plt.close()

