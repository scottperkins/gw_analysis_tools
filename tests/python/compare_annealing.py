import numpy as np
import matplotlib.pyplot as plt

data_ann = np.loadtxt("data/post_anneal.csv",delimiter=',',unpack=True)
data_expl = np.loadtxt("data/post_explore.csv",delimiter=',',unpack=True)
for x in range(len(data_ann)):
    plt.plot(data_ann[x],label="Ann")
    plt.plot(data_expl[x],label="Expl")
    plt.legend()
    plt.show()
    plt.close()
