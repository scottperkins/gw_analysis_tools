import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/bayesline_raw_data.csv",delimiter=',').T
estimate = np.loadtxt("data/bayesline_estimate_psd.csv",delimiter=',')
print(estimate)
plt.loglog(data[0], data[1]*data[1] + data[2]*data[2])
plt.loglog(data[0], estimate)
plt.savefig("plots/bayesline_raw_data.pdf")
plt.close()



#data = np.loadtxt("data/initial_position.csv",delimiter=',')
#plot_data = data[3:200*2:2]
#freqs = data[2:200*2-1:2]
#plt.semilogy(freqs,np.exp(plot_data))
#plt.show()
#plt.close()
