import numpy as np
import matplotlib.pyplot as plt
from scipy.signal.windows import tukey

data = np.loadtxt("data/whitened_data.csv",delimiter=',',unpack =True)
#dataN = np.loadtxt("data/testing_data.dat",delimiter=',',unpack =True)
#dataN2 = np.loadtxt("data/testing_PSDs.dat",delimiter=',',unpack =True)
lw=.5
plt.loglog(data[0],(data[1]),label="Hanford PSD",linewidth=lw)
plt.loglog(data[0],(data[2]),label="Livingston PSD",linewidth=lw)
plt.loglog(data[0],(data[3]**2+data[4]**2),label="Hanford data",linewidth=lw)
plt.loglog(data[0],(data[5]**2 +data[6]**2),label="Livingston data ",linewidth=lw)
plt.xlim([10,350])
plt.legend()
plt.savefig("plots/raw_whitened.pdf")
#plt.show()
plt.close()



plt.loglog(data[0],data[4]**2 + data[3]**2,label="Hanford data strain")
plt.loglog(data[0],data[5]**2 + data[6]**2,label="Livingston data strain")
#plt.xlim([40,1000])
plt.legend()
plt.savefig("plots/strains.pdf")
#plt.show()
plt.close()

window = tukey(len(data[0]),alpha=2.*.4/8)
ts = np.linspace(0,1./(data[0][1]-data[0][0]),len(data[0]))
df = data[0][1]-data[0][0]
dt = ts[1]-ts[0]
factor=len(data[3])*df*np.sqrt(dt)

inter = ( data[3] + 1j*data[4])/np.sqrt(data[1])
#print(np.where(data[0] <500, data[0], np.zeros(len(data[0]))))
inter2 = np.where( ( data[0] <350) & (data[0]> 15) , inter, 0) 
whitenedH = np.fft.ifft( inter)*factor
whitenedH2 = np.fft.ifft( inter2)*factor
#whitenedH = np.fft.ifft( inter)
#whitenedH2 = np.fft.ifft( inter2)
#h_in = int((3./4. + .05/8)*len(ts))
#l_in = int((3./4. - .1/8)*len(ts))
h_in = len(ts)
l_in = 0
#plt.plot(ts[l_in:h_in],whitenedH[l_in:h_in],linewidth=lw)
plt.plot(ts[l_in:h_in],whitenedH2[l_in:h_in],linewidth=lw,linestyle='-.',label="D1")
#plt.legend()
#plt.savefig("plots/whitened_plot1.pdf")
#plt.show()
#plt.close()
inter = ( data[5] + 1j*data[6])/np.sqrt(data[2])
inter2 = np.where( ( data[0] <350) & (data[0]> 15) , inter, 0) 
whitenedL = np.fft.ifft( inter)*factor
whitenedL2 = np.fft.ifft( inter2)*factor
#whitenedL = np.fft.ifft( inter)
#whitenedL2 = np.fft.ifft( inter2)
#plt.plot(ts[l_in:h_in],whitenedL[l_in:h_in],linewidth=lw)
plt.plot(ts[l_in:h_in],whitenedL2[l_in:h_in],linewidth=lw,linestyle='-.',label="D2")
plt.legend()
plt.savefig("plots/whitened_plot.pdf")
#plt.show()
plt.close()
