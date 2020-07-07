import numpy as np

data = np.loadtxt("Voyager.dat")
np.savetxt("Voyager.csv",data, delimiter=',')
