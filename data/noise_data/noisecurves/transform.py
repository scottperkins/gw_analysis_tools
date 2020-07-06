import numpy as np

data = np.loadtxt("ETDXylophoneDwyer.dat")
np.savetxt("ETDXylophoneDwyer.csv",data, delimiter=',')
