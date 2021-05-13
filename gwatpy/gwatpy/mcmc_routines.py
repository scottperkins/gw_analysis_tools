import multiprocessing as mp
import ctypes
import gwatpy.config as cf
import h5py
import numpy as np
import matplotlib.pyplot as plt
import emcee
import gwatpy.waveform_generator as gwg
import gwatpy.util as gpu
from scipy.signal.windows import tukey
from functools import partial

from scipy.stats import  dirichlet,kde, spearmanr
from scipy.optimize import NonlinearConstraint,differential_evolution
from scipy.special import hyp2f1, beta, gamma, loggamma,betaln
import math

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
    ctypes.c_void_p,
    ctypes.c_bool
    ]
rlib.MCMC_prep_params_py.restype = ctypes.c_char_p

rlib.pack_local_mod_structure_py.argtypes = [
    ctypes.c_void_p,
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
    ctypes.c_char_p,
    ctypes.c_void_p,    
    ctypes.c_void_p
    ]
rlib.pack_local_mod_structure_py.restype = None


rlib.mcmc_data_interface_py.argtypes = [\
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_bool]
rlib.mcmc_data_interface_py.restype = ctypes.c_void_p

rlib.mcmc_data_interface_destructor_py.argtypes = [ctypes.c_void_p]
rlib.mcmc_data_interface_destructor_py.restype=None

rlib.MCMC_likelihood_extrinsic_py.argtypes = [
    ctypes.c_bool,
    ctypes.c_void_p,
    ctypes.c_char_p,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_char_p,
    ctypes.c_bool,
    ctypes.c_char_p,
    ctypes.c_int
    ]
rlib.MCMC_likelihood_extrinsic_py.restype=ctypes.c_double

rlib.MCMC_likelihood_extrinsic_pyv2.argtypes = [
    ctypes.c_bool,
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_void_p,
    ctypes.c_int,
    ctypes.c_char_p,
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_char_p,
    ctypes.c_bool,
    ctypes.c_char_p,
    ctypes.c_int
    ]
rlib.MCMC_likelihood_extrinsic_pyv2.restype=ctypes.c_double
##########################################################
PhenomD_labels = [r"$\alpha$",r"$\sin\delta$", r"$\psi$",r"$\cos \iota$",r"$\phi_{ref}$",r"$t_c$",r"$\ln D_L$",r"$\ln \mathcal{M}$",r"$\eta$",r"$\chi_1$",r"$\chi_2$"]
PhenomD_transformedv1_labels = [r"$\alpha$",r"$\sin\delta$",r"$\cos \iota$", r"$\psi$",r"$\phi_{ref}$",r"$t_c$",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_1$",r"$\chi_2$"]
PhenomD_transformedv2_labels = [r"$\alpha$",r"$\sin\delta$",r"$\cos \iota$", r"$\psi$",r"$\phi_{ref}$",r"$t_c$",r"$D_L$",r"$m_1$",r"$m_2$",r"$\chi_1$",r"$\chi_2$"]

PhenomPv2_labels = [r"$\alpha$",r"$\sin\delta$", r"$\psi$",r"$\cos \iota$",r"$\phi_{ref}$",r"$t_c$",r"$\ln D_L$",r"$\ln \mathcal{M}$",r"$\eta$",r"$a_1$",r"$a_2$",r"$\cos \theta_1$",r"$\cos \theta_2$",r"$\phi_1$",r"$\phi_2$"]
PhenomPv2_transformedv1_labels = [r"$\alpha$",r"$\sin\delta$",r"$\cos \iota$", r"$\psi$",r"$\phi_{ref}$",r"$t_c$",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$a_1$",r"$a_2$",r"$\cos \theta_1$",r"$\cos \theta_2$",r"$\phi_1$",r"$\phi_2$"]
PhenomPv2_transformedv2_labels = [r"$\alpha$",r"$\sin\delta$",r"$\cos \iota$", r"$\psi$",r"$\phi_{ref}$",r"$t_c$",r"$D_L$",r"$m_1$",r"$m_2$",r"$a_1$",r"$a_2$",r"$\cos \theta_1$",r"$\cos \theta_2$",r"$\phi_1$",r"$\phi_2$"]
def MCMC_likelihood_extrinsic_pyv2(save_waveform, parameters, mod_struct,generation_method, frequencies, data,  psd, weights, integration_method, log10F, detectors ):
    ndet = len(data)
    length = len(data[0])
    d_type = ctypes.c_double * (ndet*length) 
    c_type = ctypes.c_char_p * (ndet) 
    i_type = ctypes.c_int * (ndet) 

    f = frequencies.reshape(ndet*length)
    dR = np.real(data).reshape(ndet*length)
    dI = np.imag(data).reshape(ndet*length)
    PSD = psd.reshape(ndet*length)
    WEIGHTS = weights.reshape(ndet*length)
    data_lengths = np.asarray([len(x) for x in data])
    dim = len(parameters)
    p_type = ctypes.c_double * (dim) 
    return rlib.MCMC_likelihood_extrinsic_pyv2(
        save_waveform,
        p_type(* (np.ascontiguousarray(parameters,dtype=ctypes.c_double))),
        mod_struct.obj,
        dim,
        generation_method.encode('utf-8'), 
        i_type(* ( np.ascontiguousarray(data_lengths,dtype=ctypes.c_int))), 
        d_type(* ( np.ascontiguousarray(f,dtype=ctypes.c_double))), 
        d_type(* ( np.ascontiguousarray(dR,dtype=ctypes.c_double))), 
        d_type(* ( np.ascontiguousarray(dI,dtype=ctypes.c_double))), 
        d_type(* ( np.ascontiguousarray(PSD,dtype=ctypes.c_double))), 
        d_type(* ( np.ascontiguousarray(WEIGHTS,dtype=ctypes.c_double))), 
        integration_method.encode('utf-8'), 
        log10F, 
        detectors.encode('utf-8'),
        ndet 
        ) 

def MCMC_likelihood_extrinsic_py(save_waveform, parameters, generation_method, frequencies, data,  psd, weights, integration_method, log10F, detectors ):
    ndet = len(data)
    length = len(data[0])
    d_type = ctypes.c_double * (ndet*length) 
    c_type = ctypes.c_char_p * (ndet) 
    i_type = ctypes.c_int * (ndet) 

    f = frequencies.reshape(ndet*length)
    dR = np.real(data).reshape(ndet*length)
    dI = np.imag(data).reshape(ndet*length)
    PSD = psd.reshape(ndet*length)
    WEIGHTS = weights.reshape(ndet*length)
    data_lengths = np.asarray([len(x) for x in data])
    return rlib.MCMC_likelihood_extrinsic_py(
        save_waveform,
        parameters.obj,
        generation_method.encode('utf-8'), 
        i_type(* ( np.ascontiguousarray(data_lengths,dtype=ctypes.c_int))), 
        d_type(* ( np.ascontiguousarray(f,dtype=ctypes.c_double))), 
        d_type(* ( np.ascontiguousarray(dR,dtype=ctypes.c_double))), 
        d_type(* ( np.ascontiguousarray(dI,dtype=ctypes.c_double))), 
        d_type(* ( np.ascontiguousarray(PSD,dtype=ctypes.c_double))), 
        d_type(* ( np.ascontiguousarray(WEIGHTS,dtype=ctypes.c_double))), 
        integration_method.encode('utf-8'), 
        log10F, 
        detectors.encode('utf-8'),
        ndet 
        ) 
class mcmc_data_interface_py(object):
    min_dim = 0
    max_dim = 0
    chain_id = 0
    nested_model_number = 0
    chain_number = 0
    RJ_step_width = 0.0
    burn_phase = False
    obj=None
    def __init__(self, **kwargs):
        if "min_dim" in kwargs:
            self.min_dim = kwargs["min_dim"]
        if "max_dim" in kwargs:
            self.max_dim = kwargs["max_dim"]
        if "chain_id" in kwargs:
            self.chain_id = kwargs["chain_id"]
        if "nested_model_number" in kwargs:
            self.nested_model_number = kwargs["nested_model_number"]
        if "chain_number" in kwargs:
            self.chain_number = kwargs["chain_number"]
        if "RJ_step_width" in kwargs:
            self.RJ_step_width = kwargs["RJ_step_width"]
        if "burn_phase" in kwargs:
            self.burn_phase = kwargs["burn_phase"]
        self.obj = rlib.mcmc_data_interface_py(
                self.min_dim,   
                self.max_dim, 
                self.chain_id, 
                self.nested_model_number, 
                self.chain_number, 
                self.RJ_step_width, 
                self.burn_phase)
    def __del__(self):
        rlib.mcmc_data_interface_destructor_py(self.obj)
    

def pack_local_mod_structure_py(interface, param,status,  waveform_extended,parameters,  full_struct,local_struct):
    p_type = ctypes.c_double * (len(param)) 

    s_type = ctypes.c_int * (len(status)) 
     

    rlib.pack_local_mod_structure_py(
        interface.obj,
        p_type(* ( np.ascontiguousarray(param,dtype=ctypes.c_double))), 
        s_type(* ( np.ascontiguousarray(status,dtype=ctypes.c_int))), 
        #ctypes.c_void_p,#Right now, this isn't even used.
        waveform_extended.encode('utf-8'), 
        full_struct.obj,
        local_struct.obj
        )

def MCMC_prep_params_py(param,  gen_params, dim,generation_method, mod_struct,save_gmst=True):
    p_type = ctypes.c_double * (len(param)) 
    tp = (ctypes.c_double * (len(param)) )()
     
    return_meth = rlib.MCMC_prep_params_py(
        p_type(* ( np.ascontiguousarray(param,dtype=ctypes.c_double))), 
        tp, 
        gen_params.obj, 
        dim,
        generation_method.encode('utf-8'), 
        mod_struct.obj,
        save_gmst
        )
    temp = return_meth.decode("utf-8")
    return temp, np.asarray(tp)

def repack_parameters_py(param,  gen_params, generation_method, dim):
    p_type = ctypes.c_double * (len(param)) 
    return rlib.repack_parameters_py(
        p_type(* ( np.ascontiguousarray(param,dtype=ctypes.c_double))), gen_params.obj, 
        generation_method.encode('utf-8'), 
        dim)

def populate_dummy_vals(kwargs):
    if kwargs["Nmod"] !=0:
        kwargs["betappe"] = np.zeros(kwargs["Nmod"])
    if kwargs["Nmod_phi"] !=0:
        kwargs["delta_phi"] = np.zeros(kwargs["Nmod_phi"])
    if kwargs["Nmod_sigma"] !=0:
        kwargs["delta_sigma"] = np.zeros(kwargs["Nmod_sigma"])
    if kwargs["Nmod_beta"] !=0:
        kwargs["delta_beta"] = np.zeros(kwargs["Nmod_beta"])
    if kwargs["Nmod_alpha"] !=0:
        kwargs["delta_alpha"] = np.zeros(kwargs["Nmod_alpha"])
def transform_mcmc_mod_struct_to_gen_param(mod_struct):
    kwargs = {"bppe":mod_struct.bppe,
        "Nmod":mod_struct.ppE_Nmod,
        "Nmod_phi":mod_struct.gIMR_Nmod_phi ,
        "Nmod_sigma":mod_struct.gIMR_Nmod_sigma ,
        "Nmod_beta":mod_struct.gIMR_Nmod_beta ,
        "Nmod_alpha":mod_struct.gIMR_Nmod_alpha ,
        "phii":mod_struct.gIMR_phii ,
        "sigmai":mod_struct.gIMR_sigmai ,
        "betai":mod_struct.gIMR_betai ,
        "alphai":mod_struct.gIMR_alphai 
    }
    return kwargs
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
    chains = [chain for chain in chains if "CHAIN" in chain ]
    chains_N = len(chains)
    data = f["MCMC_OUTPUT"][chains[-1]][local_trim::local_ac]
    dim = len(data[0])

    cols = 2
    dims_per_fig = 5
    plots_per_fig = dims_per_fig * 2
    rows = dim
    if rows > dims_per_fig:
        rows = dims_per_fig
     
    fig=[]
    axes=[]
    for x in range(int(math.ceil(dim/dims_per_fig))):
        fig_t,ax_t = plt.subplots(nrows=rows,ncols=cols,figsize=(15,15))
        fig.append(fig_t)
        axes.append(ax_t)
    alpha = .6;
    iterations = 20
    step = int(len(data)/iterations)
    means = np.zeros((dim,iterations))
    variances = np.zeros((dim,iterations))
    pts = np.zeros(iterations)
    data_full = None
    time_steps_full = None
    lengths= []
    ids= [0]
    
    for j in chains:
        data = f["MCMC_OUTPUT"][str(j)][local_trim::local_ac]
        lengths.append(len(data))
        ids.append(ids[-1]+lengths[-1])
        if data_full is None:
            data_full = data
            time_steps_full = np.arange(len(data))
        else:
            data_full = np.append(data_full, data, axis=0)
            time_steps_full = np.append(time_steps_full, np.arange(len(data)))
        print(np.shape(data_full))
        print(np.shape(time_steps_full))
        for x in np.arange(len(data[0])):
            meanT = np.mean(data[:,x])
            varT = np.var(data[:,x])
            for y in range(iterations):
                mean = abs(np.mean(data[y*step: (y+1)*step,x]))/meanT
                means[x,y]+= (mean)
                var = np.var(data[y*step: (y+1)*step,x])/varT
                variances[x,y]+=(var)
                #frac_diff = abs((mean-data[:,x])/mean)
                pt = (y*step+(y+1)*step)/2
                pts[y]=(pt)
    
    for x in np.arange(dim):
        for y in range(iterations):
            means[x,y]/=chains_N
            variances[x,y]/=chains_N
    for x in np.arange(dim):
        z = int(x%dims_per_fig)
        k = int(x/(dims_per_fig))
        print("Dimension: ",x)
        axes[k][z,0].scatter(pts,means[x,:],alpha=alpha,color="black",label='mean')
        axes[k][z,0].scatter(pts,variances[x,:],alpha=alpha,color="blue",label='variance')
        axes[k][z,0].set_ylabel("Dim {}".format(x))
        #for l in np.arange(len(lengths)):
        #   axes[k][z,1].hexbin(np.arange(len(data_full[ids[l]:ids[l+1],x])),data_full[ids[l]:ids[l+1],int(x)],alpha=alpha,color="black")
        axes[k][z,1].hexbin(time_steps_full,data_full[:,int(x)],alpha=alpha,gridsize=50,mincnt=1,color="black")
        axes[k][z,0].xaxis.set_ticks([])
        axes[k][z,1].xaxis.set_ticks([])
    for x in np.arange(int(dim / (dims_per_fig))): 
        print(x)
        axes[x][0,0].legend()
    return fig 

def trim_thin_file(filename,trim=None, ac=None, recalc_ac=False,calc_correlation=False):
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
    data_correlations = []
    if(calc_correlation):
        data_correlations.append(data)
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
            newdat = f["MCMC_OUTPUT"][chains[x+1]][int(trim_local)::int(ac_local),:]
            data = np.insert(data,-1, newdat,axis=0)
            if(calc_correlation):
                data_correlations.append(newdat)

    if(calc_correlation):
        ccmat = np.zeros((len(data_correlations),len(data_correlations)))
        pvalmat = np.ones((len(data_correlations),len(data_correlations)))
        for d in np.arange(len(data_correlations)):
            for y in np.arange(d):
                ccs = np.ones(len(data_correlations[0][0]))
                pvals = np.ones(len(data_correlations[0][0]))
                for dim in np.arange(len(data_correlations[0][0])):
                    cc,pval = spearmanr(data_correlations[d][dim],data_correlations[y][dim])
                    ccs[dim] = cc
                    pvals[dim] = pval
                print(d,y,np.amax(abs(ccs)),np.amin(pvals))
                ccmat[d,y] = np.amax(abs(ccs))
                pvalmat[d,y] =np.amin(pvals)
        for d in np.arange(len(data_correlations[0][0])):
            plt.plot(data_correlations[30][:,d])
            plt.plot(data_correlations[16][:,d])
            plt.show()
            plt.close()
        print("Max correlations: ",np.amax(ccmat))
        print("Min Pvals: ",np.amin(pval))
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

########################################################################################

def plot_injection(injection,injection_status, psd_file_in,detector, generation_method_base, psd_column=0,generation_method_extended=None,min_dim= 0, max_dim=None,threads=1 ,xlim=None,ylim=None,data_stream_file=None,mod_struct_kwargs={},gpstime=1126259462.4,figsize=None,axis=None,alpha=1):

    gmst = gpu.gps_to_GMST_radian_py(gpstime)
    psd_in = np.loadtxt(psd_file_in,skiprows=1)
    freqs = psd_in[:,0]
    psd = psd_in[:,psd_column+1]
    df = freqs[1]-freqs[0]
    T = 1./(df)
    dt = T / len(freqs)
    times = np.linspace(0,T, len(freqs))

    ax = None
    fig = None
    if axis is None:
        fig,ax = plt.subplots(nrows=1,ncols=1,figsize=figsize)
    else:
        ax = axis
    waveform_reduced=None
    mod_struct = MCMC_modification_struct_py(**mod_struct_kwargs)
    if(generation_method_extended is not None):
        waveform_reduced = partial(create_waveform_from_MCMC_output, psd=psd, freqs=freqs, min_dim=min_dim,max_dim=len(data_sub[0]), generation_method_base = generation_method_base,generation_method_extended = generation_method_extended, detector=detector, mod_struct=mod_struct,gmst=gmst)
    else:
        dim = len(injection)
        waveform_reduced = partial(create_waveform_from_MCMC_output, psd=psd, freqs=freqs, min_dim=dim,max_dim=dim, generation_method_base = generation_method_base,generation_method_extended = generation_method_base, detector=detector, mod_struct=mod_struct,gmst=gmst)

    data_packed = [injection,np.asarray(injection_status,dtype=np.int32)]

    response = waveform_reduced( data_packed)
    ax.plot(times,np.real(response),alpha=alpha,color='red' ,linewidth=1)

    if(xlim is not None):
        ax.set_xlim(xlim)
    if(ylim is not None):
        ax.set_ylim(ylim)
    return fig

#Thanks to Neil for the term 'Bayesogram..'
def plot_bayesogram(filename, psd_file_in,detector, generation_method_base, psd_column=0,generation_method_extended=None,min_dim= 0, max_dim=None,threads=1 ,xlim=None,ylim=None,data_stream_file=None,mod_struct_kwargs={},injection=None,injection_status=None,gpstime=1126259462.4,figsize=None,axis=None,color='#254159',alpha=.8,symmetric_percentile=.90,N_waveforms=1000):

    psd_in = np.loadtxt(psd_file_in,skiprows=1)
    freqs = psd_in[:,0]
    psd = psd_in[:,psd_column+1]
    df = freqs[1]-freqs[0]
    T = 1./(df)
    dt = T / len(freqs)
    times = np.linspace(0,T, len(freqs))

    gmst = gpu.gps_to_GMST_radian_py(gpstime)

    N = N_waveforms

    ax = None
    fig = None
    if axis is None:
        fig,ax = plt.subplots(nrows=1,ncols=1,figsize=figsize)
    else:
        ax = axis
    data_sub, status_sub = [],[]
    waveform_reduced=None
    mod_struct = MCMC_modification_struct_py(**mod_struct_kwargs)
    if(generation_method_extended is not None):
        data,status,model_status = RJPTMCMC_unpack_file(filename)     
        data_sub_indices = np.random.choice(np.linspace(0,len(data)-1,dtype=np.int),N)
        data_sub = data[data_sub_indices,:] 
        status_sub = status[data_sub_indices,:]
        waveform_reduced = partial(create_waveform_from_MCMC_output, psd=psd, freqs=freqs, min_dim=min_dim,max_dim=len(data_sub[0]), generation_method_base = generation_method_base,generation_method_extended = generation_method_extended, detector=detector, mod_struct=mod_struct,gmst=gmst)

    else:
        data = trim_thin_file(filename)     
        dim = len(data[0])
        data_sub_indices = np.random.choice(np.linspace(0,len(data)-1,dtype=np.int),N)
        data_sub = data[data_sub_indices,:] 
        status_sub = np.ones((N,dim))
        waveform_reduced = partial(create_waveform_from_MCMC_output, psd=psd, freqs=freqs, min_dim=dim,max_dim=dim, generation_method_base = generation_method_base,generation_method_extended = generation_method_base, detector=detector, mod_struct=mod_struct,gmst=gmst)

    data_sub_packed = [[data_sub[x],status_sub[x]] for x in np.arange(len(data_sub))]
    responses = np.zeros( (N,len(freqs)))
    #pool = mp.Pool(processes=threads)

    #parallel
    #responses = pool.map(waveform_reduced, data_sub_packed)
    #Serial
    for x in np.arange(len(data_sub_packed)):
        responses[x] = np.real(waveform_reduced( data_sub_packed[x]))
    print("DONE")

    final_response_max = np.quantile(responses,.5 + symmetric_percentile/2,axis=0 )
    final_response_min = np.quantile(responses,.5-symmetric_percentile/2,axis=0 )
    ax.fill_between(times,final_response_max,final_response_min,alpha=alpha,color=color ,linewidth=.5)
    #################################################################################
    #################################################################################
    #for i in np.arange(N):
    #    #ax.plot(times,np.real(responses[i]),alpha=.05,color='blue' ,linewidth=.1)
    #    ax.plot(times,np.real(responses[i]),alpha=.05,color=color ,linewidth=.5)
    #################################################################################
    #################################################################################

    #ninetyp = np.quantile(np.real(responses),.90,axis=0)
    #tenthp = np.quantile(np.real(responses),.1,axis=0)
    #ax.plot(times,ninetyp)
    #ax.plot(times,tenthp)

    #time_arr,responses_arr = times,np.real(responses[0])
    #for i in np.arange(N-1):
    #    time_arr = np.insert(time_arr, len(time_arr), times,axis=0)
    #    responses_arr = np.insert(responses_arr, len(responses_arr), np.real(responses[i+1]),axis=0)
    #h, xsides, ysides,im = plt.hist2d(time_arr, responses_arr,bins=100,range=[xlim,ylim],cmin=1)

    #time_arr,responses_arr = times,np.real(responses[0])
    #for i in np.arange(N-1):
    #    time_arr = np.insert(time_arr, len(time_arr), times,axis=0)
    #    responses_arr = np.insert(responses_arr, len(responses_arr), np.real(responses[i+1]),axis=0)
    #nbins = 30
    #k = kde.gaussian_kde([time_arr, responses_arr])
    #xi,yi = np.mgrid[time_arr.min():time_arr.max():nbins*1j, responses_arr.min():responses_arr.max():nbins*1j]
    #zi = k(np.vstack([xi.flatten(),yi.flatten()]))
    #ax.pcolormesh(xi,yi,zi.reshape(xi.shape))

    if injection is not None:
        data_packed = [injection,np.asarray(injection_status,dtype=np.int32)]

        response = np.real(waveform_reduced( data_packed))
        ax.plot(times,np.real(response),alpha=1,color='red' ,linewidth=1)

    if data_stream_file is not None:
        datastream = np.loadtxt(data_stream_file, skiprows=3)
        data_ft = np.fft.fft(datastream) * df
        dataplot = np.fft.ifft(data_ft / psd**.5)*dt
        ax.plot(times, dataplot,color='red')

    if(xlim is not None):
        ax.set_xlim(xlim)
    if(ylim is not None):
        ax.set_ylim(ylim)
    return fig

def create_waveform_from_MCMC_output(parameters_status,psd,freqs, min_dim,max_dim,generation_method_base,generation_method_extended, detector,mod_struct,gmst):
    parameters = parameters_status[0]
    status = parameters_status[1]
    df = freqs[1]-freqs[0]
    T = 1./(df)
    ###############################
    kwargs = transform_mcmc_mod_struct_to_gen_param(mod_struct)
    kwargs["gmst"]=gmst
    populate_dummy_vals(kwargs)
    gparam = gpu.gen_params(**kwargs)
    ###############################
    response = None

    dimct =  np.sum(status)
    
    temp_params = parameters[ status == 1]
    time_shift=0
    if(len(parameters)==min_dim):
        temp_gen,temp_temp_param = MCMC_prep_params_py(temp_params,gparam,min_dim, 
            generation_method_base,mod_struct ,save_gmst=True)

        #THIS IS A TEMPORARY FIX
        time_shift = temp_temp_param[5]
        temp_temp_param[5] = 0 
        repack_parameters_py(temp_temp_param, gparam, "MCMC_"+temp_gen, min_dim)
        response = gwg.response_generator(freqs,detector,generation_method_base, gparam)

    else:
        interface_kwargs = {"min_dim":min_dim,"max_dim":max_dim}
        interface = mcmc_data_interface_py(**interface_kwargs) 
        mod_struct_local=MCMC_modification_struct_py() 

        pack_local_mod_structure_py( interface, parameters, status, generation_method_extended, ctypes.c_void_p, mod_struct, mod_struct_local)

        temp_gen,temp_temp_param = MCMC_prep_params_py(temp_params,gparam,dimct, 
            generation_method_extended,mod_struct_local ,save_gmst=True)
        #THIS IS A TEMPORARY FIX
        time_shift = temp_temp_param[5]
        temp_temp_param[5] = 0 
        repack_parameters_py(temp_temp_param, gparam, "MCMC_"+temp_gen, dimct)
        response = gwg.response_generator(freqs,detector,generation_method_extended, gparam)
    

    response*= np.exp( 1j*(T-time_shift )*2*np.pi*freqs)
    
    window = tukey(len(response),.4*2/T)
    #window = np.ones(len(response))
    response_t = np.fft.ifft(response*window/psd**.5)*df
    #response_t = np.fft.ifft(response*window)*df
    return response_t
    
#######################################################################################
#######################################################################################

def RJcorner(data,status,figsize=None,marginal_bins=20,cov_bins=20,show_quantiles=False,titles=None,marginal_color='black',cov_color='gray',alpha=1):
    data_shape = np.shape(data)
    dim = data_shape[1]
    
    fig, axes = plt.subplots(nrows=dim,ncols=dim,figsize=figsize,sharex='col')
    
    for x in np.arange(dim):
        for y in np.arange(x+1):
            axes[x,y].grid(False)
            if y == x:
                #marginalized
                mask = status[:,x] == 1
                if(np.sum(mask) !=0):
                    axes[x,y].hist(data[mask,x],bins=marginal_bins,density=True,edgecolor=marginal_color,histtype='step', color=marginal_color,alpha=alpha)
                axes[x,y].get_yaxis().set_ticks([])
                if titles is not None:
                    if show_quantiles and np.sum(mask) != 0:
                        fifty = np.quantile(data[mask,x],.50)
                        upper = np.quantile(data[mask,x],.84)#1 sigma, for symmetric, gaussian like distribution
                        lower = np.quantile(data[mask,x],.16)#1 sigma, for symmetric, gaussian like distribution
                        #axes[x,y].set_title(str(titles[x])+r"$={0:1.2e}^{{+ {1:1.2e} }}_{{- {2:1.2e} }}$".format(fifty,upper-fifty,fifty-lower))
                        axes[x,y].set_title(r"${0:1.2e}^{{+ {1:1.2e} }}_{{- {2:1.2e} }}$".format(fifty,upper-fifty,fifty-lower))
                    else:
                        axes[x,y].set_title(str(titles[x]))
                else:
                    if show_quantiles and np.sum(mask) != 0:
                        fifty = np.quantile(data[mask,x],.50)
                        upper = np.quantile(data[mask,x],.84)#1 sigma, for symmetric, gaussian like distribution
                        lower = np.quantile(data[mask,x],.16)#1 sigma, for symmetric, gaussian like distribution
                        axes[x,y].set_title(r"${0:1.2e}^{{+ {1:1.2e} }}_{{- {2:1.2e} }}$".format(fifty,upper-fifty,fifty-lower))
            else:
                #covariance
                mask = (status[:,x] == 1) & (status[:,y] == 1)
                if(np.sum(mask) !=0):
                    axes[x,y].hexbin(data[mask,y],data[mask,x],bins=cov_bins,mincnt=1,cmap=cov_color,alpha=alpha)
                axes[x,y].get_shared_y_axes().join(axes[0,y])
                if y != 0 :
                    axes[x,y].get_yaxis().set_ticks([])
                else:
                    if titles is not None:
                        axes[x,y].set_ylabel(titles[x])
    for x in np.arange(dim-1):
        for y in np.arange(x+1,dim):
            axes[x,y].axis('off')
    if titles is not None:
        for x in np.arange(dim):
            axes[dim-1,x].set_xlabel(titles[x])
    fig.subplots_adjust(hspace=0.03,wspace=0.03)
    return fig


#######################################################################################
#######################################################################################
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

def check_initial_distribution(method,dim, distribution):
    output_distribution = []
    if(method == "IMRPhenomD" and dim == 4):
        for x in distribution:
            output_distribution.append(x)
            if output_distribution[-1][0]<0:
                output_distribution[-1][0] = 1
            if output_distribution[-1][1]<.01 or output_distribution[-1][1]>.245:
                output_distribution[-1][1] = .245
            if output_distribution[-1][2]<-.9 or output_distribution[-1][2]>.9:
                output_distribution[-1][2] = 0
            if output_distribution[-1][3]<-.9 or output_distribution[-1][3]>.9:
                output_distribution[-1][3] = 0
    if(method == "IMRPhenomD" and dim == 11):
        for x in distribution:
            output_distribution.append(x)
            if output_distribution[-1][0]<0:
                output_distribution[-1][0] = 1
            if output_distribution[-1][7]<0:
                output_distribution[-1][7] = 1
            if output_distribution[-1][8]<.01 or output_distribution[-1][8]>.245:
                output_distribution[-1][8] = .245
            if output_distribution[-1][9]<-.9 or output_distribution[-1][9]>.9:
                output_distribution[-1][9] = 0
            if output_distribution[-1][10]<-.9 or output_distribution[-1][10]>.9:
                output_distribution[-1][10] = 0

