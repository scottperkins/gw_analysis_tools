import matplotlib.pyplot as plt
import numpy  as np 

#data = np.loadtxt("data/mcmc_output.csv",delimiter=',',unpack=True)
data = np.loadtxt("data/mcmc_output_RB.csv",delimiter=',',unpack=True)
for x in data:
    plt.plot(x)
    plt.show()
    plt.close()
