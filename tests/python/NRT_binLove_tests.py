import numpy as np
import matplotlib.pyplot as plt
import corner

iterations = 8

# Loading in my data
for i in np.arange(iterations):
    q50 = np.loadtxt("data/binLove/binLove_50.csv".format(i),delimiter=',',unpack=True)
    q75 = np.loadtxt("data/binLove/binLove_75.csv".format(i),delimiter=',',unpack=True)
    q90 = np.loadtxt("data/binLove/binLove_90.csv".format(i),delimiter=',',unpack=True)

# Loading in data from arXiv:1903.03909v6
for i in np.arange(2):
    q50_comp = np.loadtxt("data/binLove/q05_fit.csv".format(i),delimiter=',',unpack=True)
    q75_comp = np.loadtxt("data/binLove/q075_fit.csv".format(i),delimiter=',',unpack=True)
    q90_comp = np.loadtxt("data/binLove/q09_fit.csv".format(i),delimiter=',',unpack=True)

# These data files are all of the form m1, m2, lambda_s, lambda_a
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 15}

plt.rc('font', **font)
# Making a plot to compare what we got from the binary Love relations with what
# arXiv:1903.03909v6 had. We used their relations, so these should match. 
plt.plot(q50[6], q50[7], 'ro', label=(r'$q_{code} = 0.50$'), markersize=5)
plt.plot(q75[6], q75[7], 'cx', label=(r'$q_{code} = 0.75$'), markersize=5)
plt.plot(q90[6], q90[7], 'g+', label=(r'$q_{code} = 0.90$'), markersize=5)
plt.plot(q50_comp[0], q50_comp[1], 'ko', label=(r'$q_{paper} = 0.50$'), markersize=5)
plt.plot(q75_comp[0], q75_comp[1], 'kx', label=r'$q_{paper} = 0.75$')
plt.plot(q90_comp[0], q90_comp[1], 'k+', label=r'$q_{paper} = 0.90$')
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize='small')
plt.xlabel(r'$\lambda_s$', fontsize=15)
plt.ylabel(r'$\lambda_a$', fontsize=15)
plt.title(r'$\lambda_a(\lambda_s)$', fontsize=20)
plt.savefig("plots/binLove", bbox_inches="tight")
plt.close()
