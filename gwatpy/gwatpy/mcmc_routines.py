import ctypes
import gwatpy.config as cf
import h5py
import numpy as np
import matplotlib.pyplot as plt
import emcee

rlib = ctypes.cdll.LoadLibrary(cf.LIB)

def plot_convergence(filename,trim=None,ac=None):
    local_trim = trim
    if trim is None:
        local_trim  = 0
    local_ac = ac
    if ac is None:
        local_ac  = 1
    f = h5py.File(filename,'r')
    chains = list(f["MCMC_OUTPUT"].keys())
    chains_N = len(chains)
    data = f["MCMC_OUTPUT"][chains[-1]][local_trim::local_ac]

    fig = plt.figure()
    ax = fig.add_subplot()
    alpha = .6;
    step = int(len(data)/50)
    for x in np.arange(len(data[0])):
    #for x in np.arange(1):
        means = []
        variances = []
        pts = []
        meanT = np.mean(data[:,x])
        varT = np.var(data[:,x])
        for y in range(50):
            print(x,y)
            mean = abs(np.mean(data[y*step: (y+1)*step,x])/meanT)
            means.append(mean)
            var = np.var(data[y*step: (y+1)*step,x])/varT
            variances.append(var)
            #frac_diff = abs((mean-data[:,x])/mean)
            pt = (y*step+(y+1)*step)/2
            pts.append(pt)
            ax.scatter(pt,mean,alpha=alpha,color="black")
            ax.scatter(pt,var,alpha=alpha,color="blue")
        #ax.set_yscale('log')
        plt.plot(pts,means,color="black")
        plt.plot(pts,variances,color="blue")
    return fig 

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
    print("trim: ",trim_local)
    print("ac: ",ac_local)
    data = f["MCMC_OUTPUT"][chains[0]][int(trim_local)::int(ac_local),:]
    for x in range(chains_N-1):
        if trim is None :
            trim_local = f["MCMC_METADATA"]["SUGGESTED TRIM LENGTHS"][x+1]
        if ac is None:
            ac_local = np.amax(f["MCMC_METADATA"]["AC VALUES"][x+1][:])
        acs=[]
        for y in range(len(data[0])):
            acs.append(emcee.autocorr.integrated_time(f["MCMC_OUTPUT"][chains[x+1]][int(trim_local)::int(ac_local),y])[0])
        print(chains[x],np.amax(acs))
        data = np.insert(data,-1, f["MCMC_OUTPUT"][chains[x+1]][int(trim_local)::int(ac_local),:],axis=0)
    #print("data shape",np.shape(data))
    return data
