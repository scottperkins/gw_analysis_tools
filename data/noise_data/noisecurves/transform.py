import numpy as np

data = np.loadtxt("hn_AdvLIGO.txt",skiprows=1)
np.savetxt("hn_AdvLIGO.csv",data, delimiter=',')
