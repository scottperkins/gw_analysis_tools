import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns; sns.set()

data = np.loadtxt("data/test2_log.csv",delimiter=',',unpack=True)
data2 = np.loadtxt("data/test_output.csv",delimiter=',',unpack=True)
plt.plot(data2[0])
plt.plot(data2[1])
#plt.show()
plt.close()
plt.plot(data[0])
#plt.show()
plt.close()
data = data[:,len(data[0])-10000:]
data2 = data2[:,len(data2[0])-10000:]
print(np.shape(data),np.shape(data2))
print(data[0,0],data2[:,0])
#plt.scatter(data2[0],data[0])
#plt.scatter(data2[1],data[0])
fig = plt.figure()
ax = fig.add_subplot(111,projection="3d")
#ax.scatter(data2[0],data2[1],data[0],label = "Likelihood")
#ax.scatter(data2[0],data2[1],data[1],label="Prior")
ax.scatter(data2[0],data2[1],data[0]+data[1],label="LL+P")
plt.legend()
plt.xlabel("RA")
plt.ylabel("DEC")
#plt.show()
plt.close()

data_grid = np.loadtxt("data/grid_search.csv",delimiter=',',unpack=True)
sns.heatmap(data_grid)
plt.show()
plt.close()
