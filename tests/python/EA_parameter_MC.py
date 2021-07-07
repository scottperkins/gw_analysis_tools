import numpy as np
import matplotlib.pyplot as plt
import corner

iterations = 13
n_bins = 50
half_bins = 25
for i in np.arange(iterations):
    #data = np.loadtxt("data/EA_parameter_MC.csv".format(i),delimiter=',',unpack=True)
    data = np.loadtxt("plots/EA_uniformdist/EA_parameter_MC.csv".format(i),delimiter=',',unpack=True)
#print(np.amax(data[11]))

# Regular histograms.
titles = [r'$c_a$', r'$c_\theta$', r'$c_\omega$', r'$c_\sigma$',
    r'$c_T$', r'$c_S$', r'$c_V$', r'$\alpha_{ppE,0}^{(2,T)}$',
    r'$g_{b1}$', r'$a_{bL}$', r'$g_{X1}$', r'$\alpha_1$',
    r'$\alpha_2$']
plot_names = ["ca", "ctheta", "cw", "csigma", "cT", "cS", "cV", "alphaT", "gb1",
    "abL", "gX1", "alpha1", "alpha2"]
for i in np.arange(iterations):
    if(i==4):
        bins = np.linspace(.5,1.5,n_bins)
        plt.hist(data[i], bins)
    else:
        plt.hist(data[i],n_bins)
    plt.title("{}".format(titles[i]))
    plt.savefig("plots/hist_{}".format(plot_names[i]))
    plt.close()

# Corner plots.
#cdata = np.hstack([data[0], data[1], data[2], data[3], data[4], data[5],
    #data[6], data[7], data[8], data[9], data[10]])
cdata = np.array([data[0], data[1], data[2], data[3], data[11], data[12], data[4],
    np.log10(data[5]), np.log10(data[6]), np.log10(abs(data[7])), data[8], data[9]])

figure = corner.corner(cdata.T, labels=[r'$c_a$', r'$c_\theta$', r'$c_\omega$',
    r'$c_\sigma$', r'$\alpha_1$', r'$\alpha_2$', r'$c_T$', r'log $c_S$', r'log $c_V$', r'log(|$\alpha_{ppE,0}^{(2,T)}$|)',
    r'$g_{b1}$', r'$a_{bL}$'],
    quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})
figure.savefig("plots/corner")
#figure.close()

"""
# Histograms with a log scale.
def plot_loghist(x, bins):
  hist, bins = np.histogram(x, bins=bins)
  logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
  plt.hist(x, bins=logbins)
  plt.xscale('log')

def plot_neg_loghist(x, bins):
    hist, bins = np.histogram(x, bins=bins)
    logbins = np.flip(-np.logspace(-24., np.log10(abs(bins[0])), len(bins)))
    #logbins = np.flip(-np.logspace(np.log10(abs(bins[-1])), np.log10(abs(bins[0])), len(bins)))
    #print("logbins: ", logbins)
    plt.hist(x, bins=logbins)
    plt.xscale('symlog', base=10, linthresh=10**-25)

def plot_symloghist(x, bins):
    hist, bins = np.histogram(x, bins=bins)
    pos_logbins = np.logspace(-24., np.log10(bins[-1]), len(bins))
    neg_logbins = np.flip(- np.logspace(-24., np.log10(abs(bins[0])), len(bins)))
    #print(pos_logbins)
    #print("negative: ", neg_logbins)
    logbins = np.concatenate((neg_logbins, pos_logbins), axis=None)
    #print("logbins: ", logbins)
    #logbins = np.logspace(np.log10(abs(bins[0])), np.log10(bins[-1]),len(bins))
    plt.hist(x, bins=logbins)
    plt.xscale('symlog', base=10, linthresh=10**-25, linscale=2)

#plt.hist(np.log10(abs(data[3])), n_bins)
#plt.savefig("plots/hist_logtest")
#plt.close()

plot_loghist(data[0], n_bins)
plt.title(r'$c_a$')
plt.savefig("plots/hist_ca")
plt.close()

plot_symloghist(data[1], half_bins)
plt.title(r'$c_\theta$')
plt.savefig("plots/hist_ctheta")
plt.close()

plot_symloghist(data[2], half_bins)
plt.title(r'$c_\omega$')
plt.savefig("plots/hist_cw")
plt.close()

plot_symloghist(data[3], half_bins)
plt.title(r'$c_\sigma$')
plt.savefig("plots/hist_csigma")
plt.close()


bins = np.linspace(.5,1.5,n_bins)
#print(bins)
plt.hist(data[4], bins=bins)
#plot_loghist(data[4], n_bins)
plt.title(r'$c_T$')
plt.savefig("plots/hist_cT")
plt.close()

plot_loghist(data[5], n_bins)
plt.title(r'$c_S$')
plt.savefig("plots/hist_cS")
plt.close()

plot_loghist(data[6], n_bins)
plt.title(r'$c_V$')
plt.savefig("plots/hist_cV")
plt.close()

plot_neg_loghist(data[7], n_bins)
plt.title(r'$\alpha_{ppE,0}^{(2,T)}$')
plt.savefig("plots/hist_alphaT")
plt.close()

plot_loghist(data[8], n_bins)
plt.title(r'$g_{b1}$')
plt.savefig("plots/hist_gb1")
plt.close()

#plt.hist(data[9], n_bins)
plot_symloghist(data[9], half_bins)
plt.title(r'$a_{bL}$')
plt.savefig("plots/hist_abL")
plt.close()

plot_symloghist(data[10], half_bins)
plt.title(r'$g_{X1}$')
plt.savefig("plots/hist_gX1")
plt.close()

plot_symloghist(data[11], half_bins)
plt.title(r'$\alpha_1$')
plt.savefig("plots/hist_alpha1")
plt.close()

plot_symloghist(data[12], half_bins)
plt.title(r'$\alpha_2$')
plt.savefig("plots/hist_alpha2")
plt.close()
"""
