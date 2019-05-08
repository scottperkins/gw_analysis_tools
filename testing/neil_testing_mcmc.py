import numpy as np
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer

data = np.loadtxt("data/mcmc_output.csv",delimiter=',',unpack=True)
burnin =1000
#data = data[burnin:]
print(data)
burnin =0
plt.plot(data[0][burnin:])
plt.plot(data[1][burnin:])
plt.show()
plt.close()
#plt.hist(data,bins=100,density=True)
#x = np.linspace(-3,3)
plt.hist2d(data[0][burnin:],data[1][burnin:],bins=200)
#plt.plot(x,np.exp(-x**2/4.)/np.sqrt(4*np.pi))
#plt.plot(x,np.exp(-x**2/10.)/np.sqrt(10*np.pi))
plt.show()
plt.close()
plt.hist(data[0][burnin:],bins=100,density=True)
plt.show()
plt.close()
plt.hist(data[1][burnin:],bins=100,density=True)
plt.show()
plt.close()
