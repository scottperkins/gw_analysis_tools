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

from scipy.stats import  dirichlet
from scipy.optimize import NonlinearConstraint,differential_evolution
from scipy.special import hyp2f1, beta, gamma, loggamma,betaln

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
    model_status = []
    if "MCMC_OUTPUT/MODEL_STATUS" in f.keys():
        model_status = f["MCMC_OUTPUT/MODEL_STATUS"][chains[0]]
        for x in range(chains_N-1):
            if( "CHAIN" in chains[x+1]):
                model_status = np.insert(model_status,-1, f["MCMC_OUTPUT/MODEL_STATUS"][chains[x+1]],axis=0)
    return data, status,model_status

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
    
def dirichlet_wrapper_full(p,*bincts):
    ptemp = np.insert(p, len(p), 1-np.sum(p))

    if np.sum(ptemp<=0):
        return 10**20
    #returnval = 0
    #for x in bincts:
    #    returnval +=dirichlet.logpdf(ptemp, x)  
    #return  -(returnval)
    returnval =dirichlet.logpdf(ptemp, np.prod(bincts,axis=0)  )
    return  -(returnval)


    

def dirichlet_wrapper(p,bincts):
    ptemp = np.insert(np.asarray(p), len(p), 1-np.sum(p))
    #print("ptemp: ",ptemp)
    returnval = 0
    if np.sum(ptemp<=0):
        return 10**20
    returnval +=dirichlet.logpdf(ptemp, bincts)  
    #return  -np.exp(returnval)
    return  -(returnval)

    #return -np.exp(dirichlet.logpdf(ptemp, mvec)+dirichlet.logpdf(ptemp,mvec2))


def emcee_wrapper(p, bincts):
    if np.sum(p<0) or  (np.sum(p) > 1):
        return -10**20
    ptemp = np.insert(np.asarray(p), len(p), 1-np.sum(p))
    #print("ptemp: ",ptemp)
    returnval = 0
    for x in bincts:
        returnval +=dirichlet.logpdf(ptemp, x)  
    #return  -np.exp(returnval)
    return  (returnval)

import emcee
import corner
import scipy
from multiprocessing import Pool

def beta_int_fn(ind_counts, ind_total_counts,dim, y):
    result = 1
    prior_val = 1 #uniform
    #prior_val = .5 #Jeffreys
    for x in np.arange(len(ind_counts)):
        result *= scipy.stats.beta.ppf(y, ind_counts[x] + prior_val, ind_total_counts[x]+dim*prior_val - ind_counts[x] - prior_val)
    return y*result

def p_mean_analytic(ind_counts, ind_total_counts):
    prior_val = 1 #uniform
    #prior_val = .5 #Jeffreys
    new_ind_counts = ind_counts + prior_val
    new_ind_total_counts = ind_total_counts + prior_val
    #prior_val = .5 #Jeffreys
    nvals = new_ind_total_counts - new_ind_counts
    N = len(new_ind_counts)

    return_val= loggamma(2 - N + np.sum(new_ind_counts)) + loggamma( 2-2*N + np.sum(nvals)+np.sum(new_ind_counts))
    return_val -= loggamma(1 - N + np.sum(new_ind_counts) )
    return_val -= loggamma(3 - 2*N + np.sum(new_ind_counts)+np.sum(nvals) )
    return np.exp(return_val)


#Returns an array of probabilities pvec at corresponding values of the parameter xvec
#data has shape (M,N) where M is the number of independent datasets and N is the length of each data set (can be different for each set)
def combine_discrete_marginalized_posteriors(datasets, bins=None):
    if bins is None:
        minval = np.amin(np.hstack(datasets))
        maxval = np.amax(np.hstack(datasets))
        bins = np.linspace(minval,maxval, 20)
    elif isinstance(bins, int):
        minval = np.amin(np.hstack(datasets))
        maxval = np.amax(np.hstack(datasets))
        bins = np.linspace(minval,maxval, bins)
    
    pvec = np.ones(len(bins)-1) 
    #pvec = np.ones(len(bins)-2) 
    pvec_arr = [] 
    bincts = np.zeros( (len(datasets), len(bins)-1))
    for x in np.arange(len(datasets)):
        binctstemp,binstemp = np.histogram(datasets[x], bins=bins)
        bincts[x] = binctstemp
    bincts = np.asarray(bincts)
    #bincts = bincts+ 1e-10
    #bincts = bincts+ 1
    
    total_counts = np.sum(bincts,axis=1)


    #with Pool(5) as p :
    #    pvec = p.map(lambda x: scipy.integrate.quad(lambda y: beta_int_fn(bincts[:,x],total_counts, len(pvec), y), 0,1)[0], np.arange(len(pvec)))
    for x in np.arange(len(pvec)):
        #pvec[x],err = scipy.integrate.quad(lambda y: beta_int_fn(bincts[:,x],total_counts, len(pvec) ,y), 0,1)
        pvec[x]=p_mean_analytic(bincts[:,x], total_counts)
    pvec /= np.sum(pvec)

    #constr = NonlinearConstraint(lambda x : 1-np.sum(x), 0, 1)
    #bnds = []
    #for x in np.arange(len(bins)-2):
    #    bnds.append((0,1))

    #res = differential_evolution(dirichlet_wrapper_full,args=bincts, bounds=bnds,constraints=(constr), maxiter=10000)
    #print(res)
    #pvec = np.insert(res.x,len(res.x), 1-np.sum(res.x))

    #pvec = np.ones(len(bins)-1) 
    #for x in bincts: 
    #    print(x)
    #    res = differential_evolution(dirichlet_wrapper,args=[x], bounds=bnds,constraints=(constr),maxiter=10000)
    #    print(res)
    #    pvectemp = np.insert(res.x,len(res.x), 1-np.sum(res.x))
    #    pvec_arr.append(pvectemp/(bins[1]-bins[0]))
    #    pvec*=pvectemp
    #    pvec = pvec/(np.sum(pvec) )

    #nwalkers = 250
    #ndim = len(pvec)
    #p0=[]
    #for x in np.arange(nwalkers):
    #    p0.append(bincts[0,:-1]/np.sum(bincts[0,:])*np.random.rand(ndim))
    #sampler = emcee.EnsembleSampler(nwalkers, ndim, emcee_wrapper,threads=10, args=[bincts])

    #pos, prob, state = sampler.run_mcmc(p0,10000,progress=True)
    #
    ##print("Ac time",sampler.get_autocorr_time())
    #temp = sampler.flatchain
    ##fig = corner.corner(temp)
    ##plt.savefig("test_corner.pdf")
    ##plt.close()
    #for x in np.arange(len(temp[0])):
    #    pvec[x] = np.median(temp[:,x])
    #pvec = np.insert(pvec,len(pvec), 1-np.sum(pvec))

    pvec = pvec/(bins[1]-bins[0])
    bin_midpoints = (bins[:-1] + bins[1:])/2
    
    return pvec, bin_midpoints,bins
