import ctypes
import gwatpy.config as cf
import h5py
import numpy as np

rlib = ctypes.cdll.LoadLibrary(cf.LIB)

def trim_thin_file(filename,trim=None, ac=None):
    f = h5py.File(filename,'r')
    chains = list(f["MCMC_OUTPUT"].keys())
    chains_N = len(chains)
    trim_local=trim 
    ac_local=ac 
    if trim is None :
        trim_local = f["MCMC_METADATA"]["SUGGESTED TRIM LENGTHS"][0]
    if ac is None:
        ac_local = np.amax(f["MCMC_METADATA"]["AC VALUES"][0][:])
    #print("CHAIN 0 shape: ",np.shape(f["MCMC_OUTPUT"][chains[0]]))
    data = f["MCMC_OUTPUT"][chains[0]][int(trim_local)::int(ac_local),:]
    for x in range(chains_N-1):
        if trim is None :
            trim_local = f["MCMC_METADATA"]["SUGGESTED TRIM LENGTHS"][x+1]
        if ac is None:
            ac_local = np.amax(f["MCMC_METADATA"]["AC VALUES"][x+1][:])
        data = np.insert(data,-1, f["MCMC_OUTPUT"][chains[x+1]][int(trim_local)::int(ac_local),:],axis=0)
    #print("data shape",np.shape(data))
    return data
