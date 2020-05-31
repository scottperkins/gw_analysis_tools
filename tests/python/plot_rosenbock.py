import numpy as np
import h5py
import matplotlib.pyplot as plt
import corner

#data = np.loadtxt("data/mcmc_output_RB.csv",delimiter=',')
f = h5py.File("test_flat.hdf5",'r')
ff = h5py.File("test.hdf5",'r')
datasuper = ff["MCMC_OUTPUT"]["CHAIN 0"]
data_old = np.loadtxt("data/mcmc_output_RB.csv",delimiter=',',unpack=True)
#print(len(ff["MCMC_OUTPUT"]["CHAIN 0"]))
#exit()
data = f["THINNED_MCMC_OUTPUT"]["THINNED FLATTENED CHAINS"]
dataT = np.transpose(data)
for x in range(len(dataT)):
    plt.plot(dataT[x])
    plt.plot(data_old[x])
    plt.show() 
    plt.close() 
fig = corner.corner(data)
plt.savefig("plots/mcmc_RB.pdf")
plt.close()

