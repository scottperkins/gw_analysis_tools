import numpy as np
import matplotlib.pyplot as plt
#from gwatpy.util import T_year_py as T

T = 24*3600*365
data1 = np.loadtxt("data/times_N_MBH.csv",delimiter=',',unpack=True)
data2 = np.loadtxt("data/times_AD_MBH.csv",delimiter=',',unpack=True)
data3 = np.loadtxt("data/times_0PN_MBH.csv",delimiter=',',unpack=True)
plt.loglog(data1[0],data1[1]/T*365+400,label="N")
plt.loglog(data2[0],data2[1]/T*365+400,label="AD")
plt.loglog(data3[0],data3[1]/T*365+400,label="0PN")
plt.legend()
plt.xlabel("freqs")
plt.ylabel("time (days)")
plt.savefig("plots/times_MBH.pdf")
plt.close()


plt.loglog(data1[1]/T,abs(data3[1]-data2[1])*2/(abs(data3[1])+abs(data2[1])),label="0pnAD")
plt.loglog(data1[1]/T,abs(data3[1]-data1[1])*2/(abs(data3[1])+abs(data1[1])),label="0pnN")
#plt.plot(data1[1]/T,(data3[1]-data2[1])/T*365,label="0pnAD")
#plt.plot(data1[1]/T,(data1[1]-data2[1])/T*365,label="NAD")
plt.xlabel("Time (yrs) -- merger at 0")
#plt.ylabel("error in days")
plt.ylabel("Fractional Error")
plt.legend()
plt.savefig("plots/times_fractional_diff_MBH.pdf")
plt.close()
