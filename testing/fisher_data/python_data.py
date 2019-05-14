import numpy as np
from phenompy.gr import IMRPhenomD_detector_frame as imrdf
from phenompy.utilities import s_solm, mpc
from time import time

parameters = np.loadtxt("parameters.csv",delimiter=',')
num_samples = len(parameters)
#num_samples = 10
dimension = len(parameters[0])
fishers = []
times = []
freqs = np.linspace(15,1000,3000)
for i in np.arange(num_samples):
    start = time()
    m1 = parameters[i][0]*s_solm
    m2 = parameters[i][1]*s_solm
    s1 = parameters[i][2]
    s2 = parameters[i][3]
    DL = parameters[i][4]*mpc
    tc = parameters[i][5]
    phic = parameters[i][6]
    model = imrdf(mass1=m1,mass2=m2,Luminosity_Distance=DL, spin1=s1, spin2=s2, collision_time = tc, collision_phase = phic)
    #fish, invfish = model.calculate_fisher_matrix_vector("Hanford_O1")
    fish,ifish = model.calculate_fisher_matrix_vector("Hanford_O1",lower_freq=freqs[0], upper_freq=freqs[-1], stepsize=freqs[1]-freqs[0])
    fishers.append(fish)
    end = time()
    times.append(end-start)
    print(i)
fishers_trans = np.zeros((num_samples, dimension*dimension)) 
for i in np.arange(num_samples):
    for j in np.arange(dimension):
        for k in np.arange(dimension):
            fishers_trans[i][j*dimension+k] = fishers[i][j][k]
np.savetxt("python.csv",fishers_trans) 
np.savetxt("py_times.csv",np.asarray(times))
    
