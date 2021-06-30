import numpy as np
import matplotlib.pyplot as plt
from corner import corner


data = np.loadtxt("data/test_mvn.csv",delimiter=',')
cov = np.loadtxt("data/test_mvn_cov.csv",delimiter=',')
mean = np.loadtxt("data/test_mvn_mean.csv",delimiter=',')
fake_data = np.random.multivariate_normal(mean=mean, cov=cov, size=len(data))
fig = corner(data)
corner(fake_data,fig=fig, color='blue')
plt.savefig("plots/mvn_test.pdf")
plt.close()
