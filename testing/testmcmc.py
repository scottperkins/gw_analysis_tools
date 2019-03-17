import gw_analysis_tools_py.mcmc_routines_ext as mcmc
import gw_analysis_tools_py.waveform_generator_ext as wg
import numpy as np
import matplotlib.pyplot as plt
from phenompy.utilities import calculate_mass1, calculate_mass2, mpc, s_solm 
from phenompy.gr import IMRPhenomD_detector_frame as imr
from phenompy.analysis_utilities import log_likelihood_maximized_coal_Full_Param as ll 

chirpmass = 20
length =2.4e4
eta = .20
mass1 = calculate_mass1(chirpmass, eta)
mass2 = calculate_mass2(chirpmass, eta)
print(mass1,mass2)
DL = 100.2
spin1 = np.ascontiguousarray([0,0,.3],dtype=np.float64)
spin2 = np.ascontiguousarray([0,0,.5],dtype=np.float64)
phic = 0.0
tc = 0.0
bppe = -0
betappe = 0.
theta = .1 
phi = 3 
incl_angle =2.0
freqs = np.ascontiguousarray(np.linspace(10,1000,length),dtype = np.float64)
noiseroot, noisefunc, temp = imr.populate_noise("Hanford_O1", int_scheme="quad")
psd  = np.ascontiguousarray(noisefunc(freqs)**2, dtype=np.float64)
params = wg.gen_params_py(mass1,mass2,DL,spin1,spin2,phic,tc,bppe,betappe,theta,phi,incl_angle,False)

#waveform, hcross= wg.fourier_waveform_polarizations_py(freqs, b"ppE_IMRPhenomD_IMR",params)
#plt.plot(freqs,waveform.real)

waveform, hcross = wg.fourier_waveform_polarizations_py(freqs, b"IMRPhenomD",params)
#plt.plot(freqs,waveform.real)

#model = imr(mass1=mass1*s_solm, mass2=mass2*s_solm,spin1=spin1[2],spin2=spin2[2],collision_phase=phic, collision_time=tc, Luminosity_Distance=DL*mpc)
#a,p,h = model.calculate_waveform_vector(freqs)
#plt.plot(freqs,h)
#plt.show()
#plt.close()

plan = mcmc.fftw_outline_py(len(freqs))
mcmc.initiate_likelihood_function_py(plan, len(freqs))

frac_diff = [[],[]]
llvec = [[],[],[]]

chirpmasses = np.linspace(5,50,10,dtype=np.float64)
#for i in [.1,.2,.5,.6,.7,1.,2,3,5,100]:
counter =0
data_real = np.ascontiguousarray(waveform.real,dtype=np.float64)
data_imag = np.ascontiguousarray(waveform.imag,dtype=np.float64)
#print(data_real[0],data_real[-1])
#print(data_imag[0],data_imag[-1])
#print(psd[0],psd[-1])
for i in chirpmasses:
    print(counter)
    counter+=1
    mass1 = np.asarray([calculate_mass1(i, eta)],dtype=np.float64)[0]
    mass2 = np.asarray([calculate_mass2(i, eta)],dtype=np.float64)[0]
    #print("Test File: ",mass1,mass2,spin1[2],spin2[2])
    params = wg.gen_params_py(mass1,mass2,DL,spin1,spin2,phic,tc,bppe,betappe,theta,phi,incl_angle,False)
    
    oldll = mcmc.maximized_coal_log_likelihood_IMRPhenomD_Full_Param_py(freqs,data_real, data_imag, psd,i, eta, spin1[2],spin2[2],DL, theta, phi, incl_angle,False,plan)  
    
    oldoldll = ll(waveform, freqs, psd, i*s_solm, eta, spin1[2], spin2[2],DL*mpc, theta,phi, incl_angle, 0.0,bppe, False)
    
    newll = mcmc.maximized_coal_Log_Likelihood_py(data_real,data_imag,psd,freqs,params,b"Hanford",b"IMRPhenomD",plan) 
    print(newll)
    print(oldll)
    print("OLDOLD",oldoldll)
    
    newllppe = mcmc.maximized_coal_Log_Likelihood_py(data_real,data_imag,psd,freqs,params,b"Hanford",b"ppE_IMRPhenomD_IMR",plan)
    #print("OLD",oldll)
    #print("NEW",newll)
    llvec[0].append(oldll)
    llvec[1].append(newll)
    llvec[2].append(oldoldll)
    #frac_diff[0].append((newll-oldll)/oldll)
    frac_diff[0].append((newll-oldll)/(newll))
    #print(mass1,mass2,i,frac_diff[0][-1])
    
    #print("outside pyx file")
    #print(mass1)
    #print(mass2)
    #print(DL)
    #print(phic)
    #print(tc)
    #print(bppe)
    #print(betappe)
    #print(incl_angle)
    #print(theta)
    #print(phi)
    frac_diff[1].append((newllppe-oldll)/oldll)
print(frac_diff[0])
plt.plot(chirpmasses,frac_diff[0])
plt.show()
plt.close()
plt.plot(chirpmasses,llvec[0],label='old')
plt.plot(chirpmasses,llvec[1],label='new')
plt.plot(chirpmasses,llvec[2],label='oldold')
#plt.plot(chirpmasses,frac_diff[1])
plt.legend()
plt.show()
plt.close()





mcmc.deactivate_likelihood_function_py(plan)
