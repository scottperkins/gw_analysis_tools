import ctypes
import gwatpy.config as cf
import h5py
import numpy as np
import matplotlib.pyplot as plt
import emcee
import gwatpy.waveform_generator as gwg
import gwatpy.util as gpu
from scipy.signal.windows import tukey

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

def trim_thin_file(filename,trim=None, ac=None, recalc_ac=False):
    f = h5py.File(filename,'r')
    chains = list(f["MCMC_OUTPUT"].keys())
    chains_N = len(chains)
    trim_local=trim 
    ac_local=ac 
    if trim is None :
        trim_local = f["MCMC_METADATA"]["SUGGESTED TRIM LENGTHS"][0]
    if ac is None:
        aclist = []
        for x in np.arange(len(f["MCMC_METADATA"]["AC VALUES"])):
            aclist.append(np.amax(f["MCMC_METADATA"]["AC VALUES"][x][:]))
        ac_local = np.mean(aclist)
    #print("CHAIN 0 shape: ",np.shape(f["MCMC_OUTPUT"][chains[0]]))
    print("trim: ",trim_local)
    print("ac: ",ac_local)
    data = f["MCMC_OUTPUT"][chains[0]][int(trim_local)::int(ac_local),:]
    for x in range(chains_N-1):
        if( "CHAIN" in chains[x+1]):
            if trim is None :
                trim_local = f["MCMC_METADATA"]["SUGGESTED TRIM LENGTHS"][x+1]
            if ac is None:
                ac_local = np.amax(f["MCMC_METADATA"]["AC VALUES"][x+1][:])
            if recalc_ac:
                acs=[]
                for y in range(len(data[0])):
                    acs.append(emcee.autocorr.integrated_time(f["MCMC_OUTPUT"][chains[x+1]][int(trim_local)::int(ac_local),y],tol=0)[0])
                print(chains[x],np.amax(acs),np.argmax(acs))
            data = np.insert(data,-1, f["MCMC_OUTPUT"][chains[x+1]][int(trim_local)::int(ac_local),:],axis=0)
    #print("data shape",np.shape(data))
    #data = data[::chains_N]
    return data
    
def RJPTMCMC_unpack_file(filename):
    f = h5py.File(filename,'r')
    chains = list(f["MCMC_OUTPUT"].keys())
    chains_N = len(chains)
    data = f["MCMC_OUTPUT"][chains[0]]
    status = f["MCMC_OUTPUT/STATUS"][chains[0]]
    for x in range(chains_N-1):
        if( "CHAIN" in chains[x+1]):
            data = np.insert(data,-1, f["MCMC_OUTPUT"][chains[x+1]],axis=0)
            status = np.insert(status,-1, f["MCMC_OUTPUT/STATUS"][chains[x+1]],axis=0)
    return data, status

#Thanks to Neil for the term 'Bayesogram..'
def plot_bayesogram(filename, freqs,detector, generation_method_base, generation_method_extended=None,min_dim= 0):

    psd_in = np.loadtxt("/home/sperkins/Downloads/LOSC_data/GW170729/GWTC1_GW170729_PSDs.dat.txt",skiprows=1)
    freqs = psd_in[:,0]
    psd = psd_in[:,1]

    fig,ax = plt.subplots(nrows=1,ncols=1)
    if(generation_method_extended is not None):
        data,status = RJPTMCMC_unpack_file(filename)     
        for i in np.arange(len(data)):
            if(np.sum(status[i]) > min_dim):
                gparam = map_method_output(data[i],status[i], generation_method_extended,min_dim )
            else:
                gparam = map_method_output(data[i],status[i], generation_method_base,min_dim )
    else:
        N = 100
        data = trim_thin_file(filename)     
        data_sub_indices = np.random.choice(np.linspace(0,len(data)-1,dtype=np.int),N)
        data_sub = data[data_sub_indices,:] 
        status = np.ones(len(data[0]))
        df = freqs[1]-freqs[0]
        T = 1./(df)
        dt = T / len(freqs)
        times = np.linspace(0,T, len(freqs))
        responses = np.zeros( (N,len(freqs)))
        #parallelize in the future
        for i in np.arange(N):
            gparam = map_method_output(data_sub[i],status, generation_method_base ,len(data_sub[i]))
            response = gwg.response_generator(freqs,detector,generation_method_base, gparam)
            response*= np.exp( 1j*(T-data_sub[i,5] )*2*np.pi*freqs)
            window = tukey(len(response),4*2/T)
            response_t = np.fft.ifft(response*window/psd**.5)*df
            responses[i] = response_t
        for i in np.arange(N):
            ax.plot(times,responses[i],alpha=.1,color='blue' )
            #ax.hexbin(times,responses[i],mincnt=1 )
    return fig

def map_method_output(data, status, method, min_dim):
    kwargs = {};
    if("PhenomPv2" in method):
        m1 = gpu.calculate_mass1(np.exp(data[7]),data[8])
        m2 = gpu.calculate_mass2(np.exp(data[7]),data[8])
        kwargs["mass1"] = m1
        kwargs["mass2"] = m2
        kwargs["RA"] = data[0]
        kwargs["DEC"] = np.arcsin(data[1])
        kwargs["psi"] = data[2]
        kwargs["incl_angle"] = np.arccos(data[3])
        kwargs["phiRef"] = data[4]
        kwargs["tc"] = data[5]
        kwargs["Luminosity_Distance"] = np.exp(data[6])
        a1 = data[9]
        a2 = data[10]
        ctilt1 = data[11]
        stilt1 = np.sqrt(1-ctilt1**2)
        ctilt2 = data[12]
        stilt2 = np.sqrt(1-ctilt2**2)
        phi1 = data[13]
        phi2 = data[14]
        kwargs["spin1"] = [ a1*stilt1*np.cos(phi1), a1*stilt1*np.sin(phi1), a1*ctilt1]
        kwargs["spin2"] = [ a2*stilt2*np.cos(phi2), a2*stilt2*np.sin(phi2), a1*ctilt2]
    elif("PhenomD" in method):
        m1 = gpu.calculate_mass1_py(np.exp(data[7]),data[8])
        m2 = gpu.calculate_mass2_py(np.exp(data[7]),data[8])
        kwargs["mass1"] = m1
        kwargs["mass2"] = m2
        kwargs["RA"] = data[0]
        kwargs["DEC"] = np.arcsin(data[1])
        kwargs["psi"] = data[2]
        kwargs["incl_angle"] = np.arccos(data[3])
        kwargs["phiRef"] = data[4]
        #kwargs["tc"] = data[5]
        kwargs["tc"] = 0
        kwargs["Luminosity_Distance"] = np.exp(data[6])
        kwargs["spin1"] = [ 0,0, data[9]]
        kwargs["spin2"] = [ 0,0, data[10]]
    gparams =gwg.gen_params(**kwargs) 
    return gparams
