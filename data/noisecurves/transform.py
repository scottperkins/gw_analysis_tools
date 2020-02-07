import numpy as np

data = np.loadtxt("CE2narrow.dat")
np.savetxt("CE2narrow.csv",data, delimiter=',')
