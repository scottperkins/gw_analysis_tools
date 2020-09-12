import numpy as np
import matplotlib.pyplot as plt
import gwatpy.util as gpu
import gwatpy.detector_util as gpdu
import gwatpy.mcmc_routines  
import gwatpy.waveform_generator as gwg
from time import time


#LAT, LONG, LOC, D = gpdu.get_detector_parameters_py("Hanford")
#print(LAT,LONG,LOC,D)
for x in np.arange(100):
    kwargs = {"mass1":11,"mass2":9,"Luminosity_Distance":500,"spin1":[.1,.8,.1],"Nmod":1,"bppe":[-1],"betappe":[10]}
    gp = gpu.gen_params(**kwargs)
    freq = np.linspace(1,100,1000)
    gen_meth = "ppE_IMRPhenomPv2_IMR"
    start = time()
    wfp = gwg.response_generator(freq,"Hanford", gen_meth, gp)
    
    kwargs = {"mass1":11,"mass2":9,"Luminosity_Distance":500,"spin1":[.1,.8,.1],"Nmod":1,"bppe":[-1],"betappe":[0]}
    gp2 = gpu.gen_params(**kwargs)
    wfp2 = gwg.response_generator(freq,"Hanford", gen_meth, gp2)
    end = time()
    print(end - start)


plt.semilogx(np.imag(wfp))
plt.semilogx(np.imag(wfp2),linestyle=':')
plt.show()
plt.close()


