import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("data/EA_fully_restricted_v1_waveform.csv",delimiter=',')
plt.loglog(data[0]*data[0]+data[1]*data[1])
plt.savefig("plots/EA_wf.pdf")
plt.close()
