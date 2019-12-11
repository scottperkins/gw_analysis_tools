import numpy as np

data = np.loadtxt("AdLIGOMidHigh.dat")
np.savetxt("AdLIGOMidHigh.csv",data, delimiter=',')
