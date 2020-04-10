import numpy as np
import matplotlib.pyplot as plt
import sys
for x in sys.argv[1:]:
    data = np.loadtxt(x,delimiter=',',unpack=True)
    plt.loglog(data[0],data[1])
    plt.title(x)
    plt.show()
    plt.close()
