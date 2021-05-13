import numpy as np
import gwatpy.mcmc_routines as gmcmc
from corner import corner
import matplotlib.pyplot as plt

source ="GW151226"
theory ="EdGB"
version="v1"
data = gmcmc.trim_thin_file("data/data_runs/{}/{}/{}/output_{}.hdf5".format(source,theory,version,source),ac=None,trim=None)

compute_confidence=True

data[:,6]=np.exp(data[:,6])
data[:,7]=np.exp(data[:,7])


labels = gmcmc.PhenomPv2_transformedv1_labels
if theory != "GR":
    for ct in np.arange(len(data[0]) - 15):
        labels.append(r'$\beta_{}$'.format(ct))

if compute_confidence:
    print("5% confidence, median, 95% confidence, 90% confidence")
    ct = 0
    for d in data.T:
        print(labels[ct])
        print(np.quantile(d,.05),np.median(d),np.quantile(d,.95), np.quantile(d,.9))
        ct+=1

fig = corner(data,labels=labels,bins=30,show_titles=True)
plt.savefig("plots/mcmc_data_{}_{}_{}.pdf".format(source,theory,version))
plt.close()
