#from pycbc.waveform import get_fd_waveform_sequence, get_fd_waveform, get_td_waveform
#from pycbc.waveform import fd_approximants, utils
#from pycbc.types import Array
#import pycbc
#from pycbc.types import frequencyseries
#from time import time
#from gw_analysis_tools_py import waveform_generator_ext 
from gw_analysis_tools_py import mcmc_routines_ext as mcmc
import numpy as np

acfile = b"data/py_ac.csv"
dat = b"data/mcmc_output_dCS.csv"
datstring = "data/mcmc_output_dCS.csv"
data = np.ascontiguousarray(np.loadtxt(datstring, delimiter=','))
num_segments = 50
dim = 9
length = 20000
threads = 8
target_corr = .01

#mcmc.write_auto_corr_file_from_data_file_py(acfile, dat, length, dim, num_segments, target_corr, threads)
mcmc.write_auto_corr_file_from_data_py(acfile, data, length, dim, num_segments, target_corr, threads)


print(data[0][0])
print(data[len(data)-1][dim-1])
