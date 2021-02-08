import numpy as np
import matplotlib.pyplot as plt

data_ppe_in = np.loadtxt("data/PNSeries_Raw_ppE.csv",delimiter=',')
data_ppe = data_ppe_in[0] + 1j*data_ppe_in[1]
data_PN_in = np.loadtxt("data/PNSeries_PNS_ppE.csv",delimiter=',')
data_PN = data_PN_in[0] + 1j*data_PN_in[1]
plt.plot(np.real(data_ppe))
plt.plot(np.real(data_PN))
plt.show()
plt.close()
plt.plot(np.imag(data_ppe))
plt.plot(np.imag(data_PN))
plt.show()
plt.close()

print(abs(np.imag(data_ppe)-np.imag(data_PN))/(abs(np.imag(data_ppe))+abs(np.imag(data_PN))))
print(abs(np.real(data_ppe)-np.real(data_PN))/(abs(np.real(data_ppe))+abs(np.real(data_PN))))

plt.plot(abs(np.imag(data_ppe)-np.imag(data_PN))/(abs(np.imag(data_ppe))+abs(np.imag(data_PN))))
plt.plot(abs(np.real(data_ppe)-np.real(data_PN))/(abs(np.real(data_ppe))+abs(np.real(data_PN))))
plt.show()
plt.close()
