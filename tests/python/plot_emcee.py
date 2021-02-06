import numpy  as np
import matplotlib.pyplot  as plt

dat = np.loadtxt("data/emcee_samples_multi_gaussian.csv",delimiter=',')
for x in np.arange(len(dat[0])):
    plt.plot(dat[:,x])
    plt.show()
    plt.close()
