import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/integration_bounds.csv",delimiter=',',unpack=True)

plt.loglog(data[0],data[1],label="integrand")
plt.loglog([data[0][0],data[0][-1]],[.1,.1])
#plt.loglog(data[0],data[2],label="Noise")
plt.savefig("plots/integrand.pdf")
