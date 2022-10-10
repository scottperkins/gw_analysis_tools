import numpy as np
import matplotlib.pyplot as plt
import corner

iterations = 12

for i in np.arange(iterations):
    data = np.loadtxt("data/EA_parameter_MC.csv".format(i),delimiter=',',unpack=True)

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 18}

plt.rc('font', **font)
# Plotting compactness as a function of lambda for body 1 and then body 2
# The relation is the same for both bodies so this basically doubles our data points.
# Compare to figure 4 of arXiv:1903.03909
plt.plot(data[6], data[8], 'ko', label='NS 1')
plt.plot(data[7], data[9], 'cx', label='NS 2')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$\lambda$', fontsize=15)
plt.ylabel('C', fontsize=15)
plt.title(r'$C(\lambda)$', fontsize=20)
plt.savefig("plots/C_Love", bbox_inches="tight")
plt.close()

# Plotting sensitivity as a function of compactness for body 1 and then body 2.
# The relation is the same for both bodies so this basically doubles our data points.
# Compare to figure 3 of arXiv:2104.04596v1
plt.plot(data[8], data[10], 'ko', label='NS 1')
plt.plot(data[9], data[11], 'cx', label='NS 2')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.xlabel('C')
plt.ylabel('s')
plt.title(r'$s(C)$')
plt.savefig("plots/sensitivity", bbox_inches="tight")
plt.close()
