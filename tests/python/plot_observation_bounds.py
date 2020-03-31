import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/observation_bounds.csv",delimiter=',',unpack=True)
T = 365*24*3600
fig, ax = plt.subplots(nrows = 2, ncols = 1)
ax[0].loglog(data[0],data[1],label="integrand")
ax[0].loglog([data[0][0],data[0][-1]],[.1,.1])
ax[1].loglog(data[0],data[3]/T)
ax[1].loglog([data[0][0],data[0][-1]],[4,4])
#plt.loglog(data[0],data[2],label="Noise")
plt.savefig("plots/integrand.pdf")
