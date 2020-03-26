import numpy as np
import matplotlib.pyplot as plt
from scipy.signal.windows import tukey

T = 365*24*3600
T_day = 24*3600

data0 = np.loadtxt("data/response_0tc.csv",delimiter=',',unpack=True)
data = np.loadtxt("data/response_nonzero_tc.csv",delimiter=',',unpack=True)
datac = np.loadtxt("data/response_cyclic_tc.csv",delimiter=',',unpack=True)

fig,ax = plt.subplots(nrows=2,ncols=1)
ax[0].loglog(data0[0],data0[1]*data0[1] + data0[2]*data0[2],label="0")
ax[0].loglog(data0[0],data[1]*data[1] + data[2]*data[2],linestyle='-.',label="3/4 T")
ax[0].loglog(datac[0],datac[1]*datac[1] + datac[2]*datac[2],linestyle='-.',label="cyclic")
ax[0].legend()
ax[1].semilogx(data0[0],data0[3]/T,label="0")
ax[1].semilogx(data0[0],data[3]/T,linestyle='-.',label="3/4 T")
ax[1].semilogx(datac[0],datac[3]/T,linestyle='-.',label="cyclic")
ax[1].legend()
ax[1].set_xlabel("freq (Hz)")
ax[1].set_ylabel("Time (yrs)")
ax[0].set_xlabel("freq (Hz)")
ax[0].set_ylabel("strain")
#plt.show()
plt.tight_layout()
plt.savefig("plots/strain_time_plot.pdf")
plt.close()


window = tukey(len(data[0]), alpha=2*T_day*7/T)
fig,ax = plt.subplots(nrows=2,ncols=1)
data0ft = np.real(np.fft.ifft(window*(data0[1] + 1j*data0[2])))
dataft = np.real(np.fft.ifft(window*(data[1] + 1j*data[2])))
datacft = np.real(np.fft.ifft(window*(datac[1] + 1j*datac[2])))
times = np.linspace(0,1./(data[0][1] - data[0][0]),len(data[0]))
ax[0].plot(times/T, data0ft,label="0")
ax[0].plot(times/T, dataft,linestyle='-.',label="3/4 T")
ax[0].plot(times/T, datacft,linestyle='-.',label="cyclic")
ax[1].plot(-data0[3]/T, data0ft,label="0")
ax[1].plot(-data[3]/T, dataft,linestyle='-.',label="3/4 T")
ax[1].plot(-datac[3]/T, datacft,linestyle='-.',label="cyclic")
ax[0].legend()
ax[1].legend()
ax[1].set_xlabel("Time (yrs)")
ax[1].set_ylabel("strain ")
ax[0].set_xlabel("time (yrs)")
ax[0].set_ylabel("strain")
plt.tight_layout()
plt.savefig("plots/TD_plot.pdf")
#plt.show()
plt.close()
