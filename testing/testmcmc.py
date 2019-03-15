import gw_analysis_tools_py.mcmc_routines_ext as mcmc
import gw_analysis_tools_py.waveform_generator_ext as wg
import numpy as np
from phenompy.utilities import calculate_mass1, calculate_mass2
from phenompy.gr import IMRPhenomD as imr

chirpmass = 20
eta = .24
mass1 = calculate_mass1(chirpmass, eta)
mass2 = calculate_mass2(chirpmass, eta)
DL = 100
spin1 = np.ascontiguousarray([0,.2,.5],dtype=np.float64)
spin2 = np.ascontiguousarray([0,.1,.2],dtype=np.float64)
phic = 0
tc = 0
bppe = 0
betappe = 0
theta = 0 
phi = 0
incl_angle = 0
freqs = np.ascontiguousarray(np.linspace(10,1000,100),dtype = np.float64)
noiseroot, noisefunc, temp = imr.populate_noise("Hanford_O1", int_scheme="quad")
psd  = noisefunc(freqs)**2
params = wg.gen_params_py(mass1,mass2,DL,spin1,spin2,phic,tc,bppe,betappe,theta,phi,incl_angle,False)

waveform = wg.fourier_waveform_py(freqs, b"IMRPhenomD",params)

for i in [.1,.2,.5,.6,.7,1.1,1.2]:
    chirpmass =10*i
    mass1 = calculate_mass1(chirpmass, eta)
    mass2 = calculate_mass2(chirpmass, eta)
    params = wg.gen_params_py(mass1,mass2,DL,spin1,spin2,phic,tc,bppe,betappe,theta,phi,incl_angle,False)
    
    data_real = np.ascontiguousarray(waveform.real,dtype=np.float64)
    data_imag = np.ascontiguousarray(waveform.imag,dtype=np.float64)
    
    plan = mcmc.fftw_outline_py(len(freqs))
    mcmc.initiate_likelihood_function_py(plan, len(freqs))
    oldll = mcmc.maximized_coal_log_likelihood_IMRPhenomD_Full_Param_py(freqs,data_real, data_imag, psd,chirpmass, eta, spin1[2],spin2[2],DL, theta, phi, incl_angle,False,plan)
    
    newll = mcmc.maximized_coal_Log_Likelihood_py(data_real,data_imag,psd,freqs,params,b"Hanford",b"IMRPhenomD",plan)
    
    print("OLD",oldll)
    print("NEW",newll)
    print("Percent diff: ",(newll-oldll)/oldll)
    mcmc.deactivate_likelihood_function_py(plan)
