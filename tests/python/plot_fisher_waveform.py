import numpy as np
import matplotlib.pyplot as plt

dataN = np.loadtxt("data/N_wf.csv",delimiter=',',unpack=True)
dataA = np.loadtxt("data/AD_wf.csv",delimiter=',',unpack=True)
plt.plot(dataN[0],dataN[1])
plt.plot(dataA[0],dataA[1])
plt.plot(dataN[0],dataN[2])
plt.plot(dataA[0],dataA[2])
plt.show()
plt.close()


plt.plot(dataN[0],(dataN[1]-dataA[1])*2/(abs(dataN[1])+abs(dataA[1])))
plt.plot(dataN[0],(dataN[2]-dataA[2])*2/(abs(dataN[2])+abs(dataA[2])))
plt.show()
plt.close()

plt.plot(dataN[0],(dataN[1]-dataA[1]))
plt.plot(dataN[0],(dataN[2]-dataA[2]))
plt.show()
plt.close()

AA = np.sqrt( dataA[1]*dataA[1] + dataA[2]*dataA[2])
AN = np.sqrt( dataA[1]*dataA[1] + dataA[2]*dataA[2])
plt.plot(dataN[0],AA)
plt.plot(dataN[0],AN)
plt.show()
plt.close()


plt.plot(dataN[0],(AA-AN)*2/(abs(AN) + abs(AA)))
plt.show()
plt.close()
