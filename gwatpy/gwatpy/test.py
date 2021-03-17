import numpy as np
import gwatpy.mcmc_routines as gmcmc
import gwatpy.plot as gplot
import gwatpy.util as gutil
import gwatpy.waveform_generator as gwg
import gwatpy.detector_util as gdutil

kwargs = {}
kwargs['mass1'] = 7
kwargs['mass2'] = 4
kwargs['spin1'] = [0,0,.1]
kwargs['spin2'] = [0,0,.1]
kwargs['phiRef'] = 0
kwargs['tc'] = 14
kwargs['RA'] = 3
kwargs['DEC'] = 1
kwargs['psi'] = 3
kwargs['incl_angle'] = 3
kwargs['gmst'] = 3
kwargs['Luminosity_Distance'] = 10
kwargs['f_ref'] = 20
kwargs['shift_time'] = True
kwargs['shift_phase'] = True
kwargs['equatorial_orientation'] = False
kwargs['horizon_coord'] = False
kwargs['cosmology'] = "PLANCK15"
parameters=gutil.gen_params(**kwargs)

detectors=["Hanford","Livingston","Virgo"]
freqs = np.asarray([np.arange(10,1000,1/16.) for x in np.arange(3)])
psd = np.asarray([ np.asarray(gdutil.populate_noise_py(freqs[x],"AdLIGOAPlus",48))**2 for x in np.arange(3)])
weights =np.asarray( [np.arange(len(freqs[0])) for x in np.arange(3)])
data = np.asarray([np.asarray(gwg.response_generator(freqs[x],detectors[x],"IMRPhenomD",parameters)) for x in np.arange(3)])

ll = gmcmc.MCMC_likelihood_extrinsic_py(True, parameters, "IMRPhenomD", freqs,data, psd, weights, "SIMPSONS",False,"HLV")
print(ll)
