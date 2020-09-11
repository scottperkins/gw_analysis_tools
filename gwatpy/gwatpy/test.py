import numpy as np
import gwatpy.util as gpu
import gwatpy.detector_util as gpdu
import gwatpy.mcmc_routines  
import gwatpy.waveform_generator as gwg


t = gwg.Tester(6)
t.get_s1()

#LAT, LONG, LOC, D = gpdu.get_detector_parameters_py("Hanford")
#print(LAT,LONG,LOC,D)
#kwargs = {"mass1":11,"mass2":9}
#gp = gwg.gen_params(**kwargs)
#print(gp.obj)
#gp.print_params_py()
#freq = np.linspace(1,100,100)
#gen_meth = "IMRPhenomD"
#wfp, wfc = gwg.waveform_generator(freq, gen_meth, gp)



#gp = gwg.gen_params(mass1=11)
#gp = gwg.gen_params()
