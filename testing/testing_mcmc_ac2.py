import numpy as np
import matplotlib.pyplot as plt
import gwatpy.mcmc_routines_ext as mcmc

mcmc.write_auto_corr_file_from_data_file_py(b'mcmc_ac2.csv',b'data/mcmc_pv2_chain.csv',1e5,13, 5, 0.01, 10, False);
