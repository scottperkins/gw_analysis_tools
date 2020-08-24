import numpy as np
import h5py
import matplotlib.pyplot as plt
import corner
import gwatpy.mcmc_routines as gmcmc
import scipy

#data = np.loadtxt("data/mcmc_output_RB.csv",delimiter=',')
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

f = h5py.File("data/mcmc_output_RB.hdf5",'r')


chains = list(f["MCMC_OUTPUT/LOGL_LOGP"])
testdata = f["MCMC_OUTPUT/CHAIN 0"]
for y in range(len(list(f["MCMC_OUTPUT/LOGL_LOGP"]))):
    print(chains[y])
    lllpdata = f["MCMC_OUTPUT/LOGL_LOGP"][chains[y]]
    def _line(x,m,b):
        return x*m +b
    x = np.arange(len(lllpdata[:,0]))
    popt,pcov = scipy.optimize.curve_fit(_line,x,lllpdata[:,0])
    print(popt)
    plt.plot(lllpdata[:,0],label="{} \ {}".format(y,f["MCMC_METADATA"]["CHAIN TEMPERATURES"][y]))
    plt.plot(x,_line(x,*popt))
    plt.plot(lllpdata[:,1],label='P')
    plt.plot(lllpdata[:,1]+lllpdata[:,0],label='P+L')
    plt.legend()
    plt.show()
    #plt.savefig("plots/ll_lp_d.pdf")
    plt.close()
    for x in np.arange(len(f["MCMC_OUTPUT/"+chains[y]][0])):
        plt.plot(f["MCMC_OUTPUT/"+chains[y]][:,x],label=chains[y]+"-"+str(x))
        #plt.plot(testdata2[:,x])
        plt.legend()
        plt.show()
        plt.close()
#for x in np.arange(len(testdata[0])):
#    plt.plot(testdata[:,x])
#    #plt.plot(testdata2[:,x])
#    plt.show()
#    plt.close()
exit()





data = gmcmc.trim_thin_file("data/mcmc_output_RB.hdf5",trim=None,ac=None)
print(np.shape(data))
fig = corner.corner(data)
plt.savefig("plots/mcmc_RB.pdf")
plt.close()
#fig = corner.corner(dataTrim)
#plt.savefig("plots/mcmc_RB_trim.pdf")
#plt.close()

