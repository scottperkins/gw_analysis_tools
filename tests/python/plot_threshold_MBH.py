import numpy as np
import matplotlib.pyplot as plt

T = 365*24*3600
data = np.loadtxt("data/threshold_output_MBH.csv",delimiter=',',unpack=True)
plt.semilogx(data[0]/T,data[1])
plt.semilogx(np.asarray([data[0][0],data[0][-1]])/T,[8,8])
#plt.loglog(data[0]/T,data[1])
plt.xlabel("Time (YRS)")
plt.ylabel("SNR")
plt.savefig("plots/threshold_MBH.pdf")
plt.close()
