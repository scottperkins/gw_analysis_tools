import ctypes
import gwatpy.config as cf
import h5py
import numpy as np
import matplotlib.pyplot as plt
import emcee
import gwatpy.waveform_generator as gwg
import gwatpy.util as gpu
from scipy.signal.windows import tukey
import multiprocessing as mp
from functools import partial

rlib = ctypes.cdll.LoadLibrary(cf.LIB)
##########################################################
rlib.MCMC_modification_struct_py.argtypes = \
    [
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_int,
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_bool,
    ctypes.c_bool\
    ]
rlib.MCMC_modification_struct_py.restype = ctypes.c_void_p

rlib.MCMC_modification_struct_py_destructor.argtypes = [ctypes.c_void_p]
rlib.MCMC_modification_struct_py_destructor.restype = ctypes.c_void_p

rlib.repack_parameters_py.argtypes = [\
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.c_int
    ]
rlib.repack_parameters_py.restype = None

rlib.MCMC_prep_params_py.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_void_p,    
    ctypes.c_int, 
    ctypes.c_char_p,
    ctypes.c_void_p
    ]
rlib.MCMC_prep_params_py.restype = ctypes.c_char_p

##########################################################

def MCMC_prep_params_py(param,  gen_params, dim,generation_method, mod_struct):
    p_type = ctypes.c_double * (len(param)) 
    tp = (ctypes.c_double * (len(param)) )()
     
    return_meth = rlib.MCMC_prep_params_py(
        p_type(* ( np.ascontiguousarray(param,dtype=ctypes.c_double))), 
        tp, 
        gen_params.obj, 
        dim,
        generation_method.encode('utf-8'), 
        mod_struct.obj
        ).decode()
    return return_meth, np.asarray(tp)

def repack_parameters_py(param,  gen_params, generation_method, dim):
    p_type = ctypes.c_double * (len(param)) 
    return rlib.repack_parameters_py(
        p_type(* ( np.ascontiguousarray(param,dtype=ctypes.c_double))), gen_params.obj, 
        generation_method.encode('utf-8'), 
        dim)

class MCMC_modification_struct_py(object):
    ppE_Nmod = 0
    bppe = []
    gIMR_Nmod_phi = 0
    gIMR_phii = []
    gIMR_Nmod_sigma = 0
    gIMR_sigmai = []
    gIMR_Nmod_beta = 0
    gIMR_betai = []
    gIMR_Nmod_alpha = 0
    gIMR_alphai = []
    NSflag1 = False
    NSflag2 = False
    def __init__(self, **kwargs):
        if "ppE_Nmod" in kwargs:
            self.ppE_Nmod = kwargs["ppE_Nmod"] 
        if "bppe" in kwargs:
            self.bppe = kwargs["bppe"] 
        if "gIMR_Nmod_phi" in kwargs:
            self.gIMR_Nmod_phi = kwargs["gIMR_Nmod_phi"] 
        if "gIMR_phii" in kwargs:
            self.gIMR_phii = kwargs["gIMR_phii"] 
        if "gIMR_Nmod_sigma" in kwargs:
            self.gIMR_Nmod_sigma = kwargs["gIMR_Nmod_sigma"] 
        if "gIMR_sigmai" in kwargs:
            self.gIMR_sigmai = kwargs["gIMR_sigmai"] 
        if "gIMR_Nmod_beta" in kwargs:
            self.gIMR_Nmod_beta = kwargs["gIMR_Nmod_beta"] 
        if "gIMR_betai" in kwargs:
            self.gIMR_betai = kwargs["gIMR_betai"] 
        if "gIMR_Nmod_alpha" in kwargs:
            self.gIMR_Nmod_alpha = kwargs["gIMR_Nmod_alpha"] 
        if "gIMR_alphai" in kwargs:
            self.gIMR_alphai = kwargs["gIMR_alphai"] 
        if "NSflag1" in kwargs:
            self.NSflag1 = kwargs["NSflag1"] 
        if "NSflag2" in kwargs:
            self.NSflag2 = kwargs["NSflag2"] 
        bppe_type = self.ppE_Nmod*ctypes.c_double
        phi_type = self.gIMR_Nmod_phi*ctypes.c_int
        sigma_type = self.gIMR_Nmod_sigma*ctypes.c_int
        beta_type = self.gIMR_Nmod_beta*ctypes.c_int
        alpha_type = self.gIMR_Nmod_alpha*ctypes.c_int
        self.obj = rlib.MCMC_modification_struct_py(
            self.ppE_Nmod,
            bppe_type(*(np.ascontiguousarray(self.bppe,dtype=ctypes.c_double))),
            self.gIMR_Nmod_phi,
            phi_type(*(np.ascontiguousarray(self.gIMR_phii,dtype=ctypes.c_int))),
            self.gIMR_Nmod_sigma,
            sigma_type(*(np.ascontiguousarray(self.gIMR_sigmai,dtype=ctypes.c_int))),
            self.gIMR_Nmod_beta,
            beta_type(*(np.ascontiguousarray(self.gIMR_betai,dtype=ctypes.c_int))),
            self.gIMR_Nmod_alpha,
            alpha_type(*(np.ascontiguousarray(self.gIMR_alphai,dtype=ctypes.c_int))),
            self.NSflag1,
            self.NSflag2
        )
    def __del__(self):
        rlib.MCMC_modification_struct_py_destructor(self.obj)

##########################################################
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
def plot_bayesogram(filename, psd_file_in,detector, generation_method_base, generation_method_extended=None,min_dim= 0, threads=1 ,**mod_struct_kwargs):

    psd_in = np.loadtxt(psd_file_in,skiprows=1)
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
        dim = len(data[0])

        data_sub_indices = np.random.choice(np.linspace(0,len(data)-1,dtype=np.int),N)
        data_sub = data[data_sub_indices,:] 
    
        status = np.ones(dim)

        df = freqs[1]-freqs[0]
        T = 1./(df)
        dt = T / len(freqs)
        times = np.linspace(0,T, len(freqs))
        responses = np.zeros( (N,len(freqs)))
        mod_struct = MCMC_modification_struct_py(**mod_struct_kwargs)
    
        waveform_reduced = partial(create_waveform_from_MCMC_output, psd=psd, freqs=freqs, dim=dim, generation_method = generation_method_base, detector=detector, mod_struct=mod_struct)
        pool = mp.Pool(processes=threads)
        responses = pool.map(waveform_reduced, data_sub)

        for i in np.arange(N):
            ax.plot(times,responses[i],alpha=.1,color='blue' )
            #ax.hexbin(times,responses[i],mincnt=1 )
    return fig

def create_waveform_from_MCMC_output(parameters,psd,freqs, dim,generation_method, detector,mod_struct):
    df = freqs[1]-freqs[0]
    T = 1./(df)
    gparam = gpu.gen_params()
    
    temp_gen,temp_param = MCMC_prep_params_py(parameters,gparam,dim, 
        generation_method,mod_struct )
    
    #THIS IS A TEMPORARY FIX
    time_shift = temp_param[5]
    temp_param[5] = 0 
    repack_parameters_py(temp_param, gparam, "MCMC_"+temp_gen, dim)
    
    response = gwg.response_generator(freqs,detector,generation_method, gparam)
    response*= np.exp( 1j*(T-time_shift )*2*np.pi*freqs)
    
    window = tukey(len(response),4*2/T)
    response_t = np.fft.ifft(response*window/psd**.5)*df
    
    return response_t
    
