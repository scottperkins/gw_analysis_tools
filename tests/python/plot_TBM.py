import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

T = 365*24*3600
data_base = np.loadtxt("data/Tbm_times_base.csv",delimiter=',',unpack=True)
data_inv_N = np.loadtxt("data/Tbm_times_inv_N.csv",delimiter=',',unpack=True)
data_inv_AD = np.loadtxt("data/Tbm_times_inv_AD.csv",delimiter=',',unpack=True)
f = interp1d(data_base[1],data_base[0])
#print(data_inv_N[0])

plt.plot(data_base[0],data_base[1]/T+.1,label="True time")
plt.plot(data_inv_N[0],data_inv_N[1]/T+.1,label="N")
plt.plot(data_inv_AD[0],data_inv_AD[1]/T+.1,linestyle='-.',label="AD")
plt.yscale('log')
plt.xscale('log')
plt.xlabel("Freqs")
plt.ylabel("Times")
plt.legend()
plt.savefig("plots/Tbm.pdf")
plt.close()

#resid_N = abs(f(data_inv_N[1]) - data_inv_N[0])
#resid_AD = abs(f(data_inv_AD[1]) - data_inv_AD[0])
resid_N = abs(data_base[0] - data_inv_N[0])
resid_AD = abs(data_base[0] - data_inv_AD[0])
plt.loglog(data_inv_N[1]/T,resid_N,label="N")
plt.loglog(data_inv_AD[1]/T,resid_AD,label="AD")
plt.xlabel("Times (years)")
plt.ylabel("Error (Hz)")
plt.legend()
plt.savefig("plots/Tbm_resid.pdf")
plt.close()

plt.loglog(data_inv_N[1]/T,resid_N/data_base[0],label="N")
plt.loglog(data_inv_AD[1]/T,resid_AD/data_base[0],label="AD")
plt.legend()
plt.savefig("plots/Tbm_resid_fractional.pdf")
plt.close()


#plt.loglog(f(data_base[1]),data_base[1]/T,label="True")
#plt.loglog(data_inv_N[0],data_inv_N[1]/T,label="N")
#plt.loglog(data_inv_AD[0],data_inv_AD[1]/T,label="AD")
#plt.legend()
#plt.savefig("plots/Tbm_interp.pdf")
#plt.close()
