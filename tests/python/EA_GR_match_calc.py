import numpy as np
import matplotlib.pyplot as plt


ID = 0

data = np.loadtxt("data/EA_GR_COMP_{}.csv".format(ID), delimiter=',').T

EAwf = data[0] + 1j*data[1]
GRwf = data[2] + 1j*data[3]
EAphase = data[4]
GRphase = data[5]

plt.loglog(abs(EAwf))
plt.loglog(abs(GRwf))
plt.savefig("plots/EA_GR_COMP_{}.pdf".format(ID))
plt.close()



plt.plot(EAphase)
plt.plot(GRphase)
plt.savefig("plots/EA_GR_COMP_phase_{}.pdf".format(ID))
plt.close()
