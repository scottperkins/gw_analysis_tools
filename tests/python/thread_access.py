import matplotlib.pyplot as plt
import numpy as np

dat = np.loadtxt('data/thread_access_test.csv',delimiter=',')
times = np.arange(len(dat))
plt.hexbin(times,dat,gridsize=(100,12),mincnt=1)
plt.savefig('plots/thread_access.pdf')
plt.close()
