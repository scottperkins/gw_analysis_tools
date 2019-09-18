import numpy as np
import matplotlib.pyplot as plt
import corner

data_import = np.loadtxt("data/mcmc_output_reduceddim.csv",delimiter=',')
for x in data_import:
    x[0] = np.exp(x[0])
labels = ["chirp","eta","chi1","chi2"]

figure = corner.corner(data_import, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("mcmc_testing_reduceddim_50k.pdf")
plt.close()
