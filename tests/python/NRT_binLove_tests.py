import numpy as np
import matplotlib.pyplot as plt
import corner

iterations = 8
iterations_residuals = 9

# Loading in my data
for i in np.arange(iterations):
    q50 = np.loadtxt("data/binLove/binLove_50.csv".format(i),delimiter=',',unpack=True)
    q75 = np.loadtxt("data/binLove/binLove_75.csv".format(i),delimiter=',',unpack=True)
    q90 = np.loadtxt("data/binLove/binLove_90.csv".format(i),delimiter=',',unpack=True)

# Loading in my data with error marginalization implemented
for i in np.arange(iterations):
    q50_error = np.loadtxt("data/binLove/binLove_error_50.csv".format(i),delimiter=',',unpack=True)
    q75_error = np.loadtxt("data/binLove/binLove_error_75.csv".format(i),delimiter=',',unpack=True)
    q90_error = np.loadtxt("data/binLove/binLove_error_90.csv".format(i),delimiter=',',unpack=True)

# Loading in my data with the residual listed explicitly
for i in np.arange(iterations_residuals):
    q50_res = np.loadtxt("data/binLove/binLove_residuals_50.csv".format(i),delimiter=',',unpack=True)
    q75_res = np.loadtxt("data/binLove/binLove_residuals_75.csv".format(i),delimiter=',',unpack=True)
    q90_res = np.loadtxt("data/binLove/binLove_residuals_90.csv".format(i),delimiter=',',unpack=True)


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

# Making a plot to compare what we have from the binary Love relations (including
# error marginalization) with what arXiv:1903.03909v6 had.
plt.plot(q50_error[6], q50_error[7], 'ro', label=(r'$q_{code} = 0.50$'), markersize=5)
plt.plot(q75_error[6], q75_error[7], 'cx', label=(r'$q_{code} = 0.75$'), markersize=5)
plt.plot(q90_error[6], q90_error[7], 'g+', label=(r'$q_{code} = 0.90$'), markersize=5)
plt.plot(q50_comp[0], q50_comp[1], 'ko', label=(r'$q_{paper} = 0.50$'), markersize=5)
plt.plot(q75_comp[0], q75_comp[1], 'kx', label=r'$q_{paper} = 0.75$')
plt.plot(q90_comp[0], q90_comp[1], 'k+', label=r'$q_{paper} = 0.90$')
plt.xscale('log')
plt.yscale('log')
plt.legend(fontsize='small')
plt.xlabel(r'$\lambda_s$', fontsize=15)
plt.ylabel(r'$\lambda_a$', fontsize=15)
plt.title(r'$\lambda_a(\lambda_s)$', fontsize=20)
plt.savefig("plots/binLove_error", bbox_inches="tight")
plt.close()

# Making a plot for comparison with figure 3 of arXiv:1804.03221v2
plt.plot(q50_error[7], q50_error[6], 'ro', label=(r'$q = 0.50$'), markersize=5)
plt.plot(q75_error[7], q75_error[6], 'cx', label=(r'$q = 0.75$'), markersize=5)
plt.plot(q90_error[7], q90_error[6], 'g+', label=(r'$q= 0.90$'), markersize=5)
plt.xlim([0,5000])
plt.ylim([0,5000])
plt.legend(fontsize='small')
plt.ylabel(r'$\lambda_s$', fontsize=15)
plt.xlabel(r'$\lambda_a$', fontsize=15)
plt.title(r'$\lambda_s(\lambda_a)$ linear scale', fontsize=20)
plt.savefig("plots/binLove_error_fig3", bbox_inches="tight")
plt.close()

#print(q50_res[8]/q50_res[7])
#print("residual min = ",np.amin((q50_res[8]-q50_res[7])/q50_res[7]))
#print("residual max = ", np.amax((q50_res[8]-q50_res[7])/q50_res[7]))
n_bins = 25
# Making plots to compare Lambda_a with the error. Note that
# q50_res[8]= new_tidal_a and q50_res[7]=tidal_a from binary love (before error addition)
plt.hist((q50_res[8]-q50_res[7])/q50_res[7], n_bins)
plt.title(r'$error/\lambda_a$, q=.5')
plt.savefig("plots/binLove_residuals_hist_50")
plt.close()

plt.hist((q75_res[8]-q75_res[7])/q75_res[7], n_bins)
plt.title(r'$error/\lambda_a$, q=.75')
plt.savefig("plots/binLove_residuals_hist_75")
plt.close()

plt.hist((q90_res[8]-q90_res[7])/q90_res[7], n_bins)
plt.title(r'$error/\lambda_a$, q=.90')
plt.savefig("plots/binLove_residuals_hist_90")
plt.close()

# Making plots to compare with figure 5 of arXiv:1903.03909v7
plt.plot( q50_res[6], abs(q50_res[8]-q50_res[7]), 'ro')
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$|\Lambda_a - \Lambda_a^{fit}|$')
plt.xlabel(r'$\Lambda_s$')
plt.title(r'$q=.50$')
plt.savefig("plots/binLove_residuals_50", bbox_inches="tight")

plt.plot(q75_res[6], abs(q75_res[8]-q75_res[7]), 'ro')
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$|\Lambda_a - \Lambda_a^{fit}|$')
plt.xlabel(r'$\Lambda_s$')
plt.title(r'$q=.75$')
plt.savefig("plots/binLove_residuals_75", bbox_inches="tight")

plt.plot(q90_res[6], abs(q90_res[8]-q90_res[7]), 'ro')
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$|\Lambda_a - \Lambda_a^{fit}|$')
plt.xlabel(r'$\Lambda_s$')
plt.title(r'$q=.90$')
plt.savefig("plots/binLove_residuals_90", bbox_inches="tight")
