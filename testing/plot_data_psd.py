import matplotlib.pyplot as plt
import numpy as np

strain = np.loadtxt("data/gw150914_data_H.csv",delimiter=',',unpack=True)
psd = np.loadtxt("data/gw150914_psd_H.csv",delimiter=',',unpack=True)
freq = np.loadtxt("data/gw150914_freq_H.csv",delimiter=',',unpack=True)

plt.loglog(freq,psd,label='neil')
#plt.savefig("neil_psd.pdf")
plt.show()
plt.close()
plt.plot(freq,strain[0], label='neil dat')
#plt.savefig("neil_strain.pdf")
plt.show()
plt.close()

strain = np.loadtxt("data/gw_150914_data.csv",delimiter=',',unpack=True)
psd = np.loadtxt("data/gw_150914_psd.csv",delimiter=',',unpack=True)
freq = np.loadtxt("data/gw_150914_freq.csv",delimiter=',',unpack=True)

plt.loglog(freq,psd,label='pycbc')
plt.legend()
plt.show()
#plt.savefig("pycbc_psd.pdf")
plt.close()
plt.plot(freq,strain[0],label='pycbc')
plt.legend()
plt.show()
#plt.savefig("pycbc_strain.pdf")
plt.close()
