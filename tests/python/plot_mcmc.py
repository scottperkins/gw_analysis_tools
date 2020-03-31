import numpy as np 
import matplotlib.pyplot as plt
import corner 


data = np.loadtxt("data/mcmc_output.csv",delimiter=',')
labels=["x","Y"]
fig= corner.corner(data, labels=labels, quantiles=[.16,.5,.84],show_titles=True)
plt.show()
plt.close()
