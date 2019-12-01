import numpy as np
import matplotlib.pyplot as plt
import gwatpy.mcmc_routines_ext as mcmc

#data =np.genfromtxt('data/mcmc_pv2_chain_data.csv')
#data = data[60000:]
#data = np.ascontiguousarray(data)
#mcmc.write_auto_corr_file_from_data_py(b'mcmc_ac2.csv',data,4e4,11, 5, 0.01, 10, False);
mcmc.write_auto_corr_file_from_data_file_py(b'mcmc_ac2.csv',b'data/mcmc_pv2_chain_data.csv',1.5e5,11, 10, 0.01, 10, False);
