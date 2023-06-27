import numpy as np

#data = np.genfromtxt("data/Klein16_PopIII.dat")
#np.savetxt("data/Klein16_PopIII.csv",data,delimiter=',')

data = np.genfromtxt("./AtomClockGW_v2.txt")
outdata = data[1:]
print(outdata.shape)
np.savetxt("./currently_supported/AtomClockGW_v2.csv",outdata,delimiter=',')
