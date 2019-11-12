import numpy as np
import corner
import matplotlib.pyplot as plt

data = np.genfromtxt("data/mcr_sampling.dat",delimiter=',')
#plt.hist(data.T[0], bins=100, density=True)
#plt.show()
#plt.close()
#plt.hist(data.T[1], bins=100, density=True)
#plt.show()
#plt.close()
#plt.hist(data.T[2], bins=100, density=True)
#plt.show()
#plt.close()

labels=["x1","x2","x3"]
figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.show()
plt.close()
