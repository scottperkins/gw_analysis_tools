import numpy as np
import matplotlib.pyplot as plt
from scipy.signal.windows import tukey

data = np.loadtxt("data/data_output.csv",delimiter=',',unpack =True)

plt.loglog(data[0],np.sqrt(data[1]),label="Hanford PSD")
plt.loglog(data[0],np.sqrt(data[2]),label="Livingston PSD")
plt.loglog(data[0],data[3],label="Hanford data Re")
plt.loglog(data[0],data[4],label="Hanford data Im")
plt.loglog(data[0],data[5],label="Livingston data Re")
plt.loglog(data[0],data[6],label="Livingston data Im")
#plt.xlim([40,1000])
plt.legend()
plt.show()
plt.close()



plt.loglog(data[0],data[4]**2 + data[3]**2,label="Hanford data strain")
plt.loglog(data[0],data[5]**2 + data[6]**2,label="Livingston data strain")
#plt.xlim([40,1000])
plt.legend()
plt.show()
plt.close()

window = tukey(len(data[0]),alpha=.4)
ts = np.linspace(0,1./(data[0][1]-data[0][0]),len(data[0]))

inter = ( data[3] + 1j*data[4])/np.sqrt(data[1])*window
#print(np.where(data[0] <500, data[0], np.zeros(len(data[0]))))
inter2 = np.where(data[0] <200, inter, 0) 
whitenedH = np.fft.ifft( inter)
whitenedH2 = np.fft.ifft( inter2)
#h_in = int((3./4. + .05)*len(ts))
#l_in = int((3./4. - .05)*len(ts))
h_in = len(ts)
l_in = 0
plt.plot(ts[l_in:h_in],whitenedH[l_in:h_in])
plt.plot(ts[l_in:h_in],whitenedH2[l_in:h_in])
plt.legend()
plt.show()
plt.close()
inter = ( data[5] + 1j*data[6])/np.sqrt(data[2])*window
inter2 = np.where(data[0] <200, inter, 0) 
whitenedL = np.fft.ifft( inter)
whitenedL2 = np.fft.ifft( inter2)
plt.plot(ts[l_in:h_in],whitenedL[l_in:h_in])
plt.plot(ts[l_in:h_in],whitenedL2[l_in:h_in])
plt.legend()
plt.show()
plt.close()
