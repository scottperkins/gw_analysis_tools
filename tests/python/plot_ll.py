import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("data/injection_ll.csv",delimiter=',',unpack=True)

plt.plot(data[0],label="ll")
plt.plot(data[1],label="lp")
plt.legend()
plt.show()
plt.close()
