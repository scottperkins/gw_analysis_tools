import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/bayesline_raw_data.csv",delimiter=',').T
print(data)
plt.loglog(data[0], data[1]*data[1] + data[2]*data[2])
plt.savefig("plots/bayesline_raw_data.pdf")
plt.close()
