import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data = np.loadtxt("data/test2_log.csv",delimiter=',',unpack=True)
data2 = np.loadtxt("data/test_output.csv",delimiter=',',unpack=True)
#plt.scatter(data2[0],data[0])
#plt.scatter(data2[1],data[0])
fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")
ax.scatter(data2[0],data2[1],data[0])
ax.scatter(data2[0],data2[1],data[1])
plt.show()
plt.close()
