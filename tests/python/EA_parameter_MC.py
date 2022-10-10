import numpy as np
import matplotlib.pyplot as plt
import corner

iterations = 13
n_bins = 50
half_bins = 25
for i in np.arange(iterations):
    data = np.loadtxt("data/EA_parameter_MC.csv".format(i),delimiter=',',unpack=True)
    #data = np.loadtxt("plots/EA_logdist/EA_parameter_MC.csv".format(i),delimiter=',',unpack=True)
#print("ctheta min = ",np.amin(data[1]))
#print("comega min = ", np.amin(data[2]))

#print("alpha1 max = ",np.amax(data[11]))
#print("alpha2 max = ",np.amax(data[12]))
#print("gX1 max = ",np.amax(data[10]))
#print("alpha2T max = ",np.amax(data[7]))
#print("abL min = ", np.amin(data[9]))
"""
cwneg = []
cwnegctheta = []
cwnegalpha1 = []
cwnegalpha2 = []
#cwneg = {} # A dictionary to hold points when cw < 0
for i in np.arange(len(data[2])):
    if data[2][i] < 0:
        cwneg.append(data[2][i])
        cwnegctheta.append(data[1][i])
        cwnegalpha1.append(data[11][i])
        cwnegalpha2.append(data[12][i])

print(len(cwneg))

#cwnegcorner = np.array([cwneg, cwnegctheta, cwnegalpha1, cwnegalpha2])
cwnegcorner = np.array([data[2], data[1], data[11], data[12]])

#print(cwnegcorner)
#cdata = np.array([data[0], data[1], data[2], data[3], data[11], data[12], data[4],
#    np.log10(data[5]), np.log10(data[6]), np.log10(abs(data[7])), data[8], data[9]])

figure = corner.corner(cwnegcorner.T, labels=[ r'$c_\omega$', r'$c_\theta$',
     r'$\alpha_1$', r'$\alpha_2$'],
    quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})
figure.savefig("plots/corner_cwneg")
"""

def plot_loghist(x, bins):
    hist, bins = np.histogram(x, bins=bins)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(x, bins=logbins)#, density=True)#, log=True)
    plt.xscale('log')

def plot_neg_loghist(x, bins, upper_bound):
    hist, bins = np.histogram(x, bins=bins)
    logbins = np.flip(-np.logspace(upper_bound, np.log10(abs(bins[0])), len(bins)))
    plt.hist(x, bins=logbins)#,density=True, log=True)
    plt.xscale('symlog', base=10, linthresh=10**upper_bound)

def plot_symloghist(x, bins, posbins, negbins, bound):
    hist, bins = np.histogram(x, bins=bins)
    pos_logbins = np.logspace(bound, np.log10(bins[-1]), posbins)
    neg_logbins = np.flip(- np.logspace(bound, np.log10(abs(bins[0])), negbins))
    logbins = np.concatenate((neg_logbins, pos_logbins), axis=None)
    plt.hist(x, bins=logbins)#,density=True, log=True)
    plt.xscale('symlog', base=10, linthresh=10**bound, linscale=2)


# Plots for data drawn from a uniform distribution.
plt.hist(data[0], n_bins)
plt.title(r'$c_a$')
plt.savefig("plots/hist_ca")
plt.close()

plt.hist(data[1], n_bins)
plt.title(r'$c_\theta$')
plt.savefig("plots/hist_ctheta")
plt.close()

plt.hist(data[2], n_bins)
plt.title(r'$c_\omega$')
plt.savefig("plots/hist_cw")
plt.close()

plt.hist(data[3], n_bins)
plt.title(r'$c_\sigma$')
plt.savefig("plots/hist_csigma")
plt.close()

bins = np.linspace(.5,1.5,n_bins)
plt.hist(data[4], bins=bins)
plt.title(r'$c_T$')
plt.savefig("plots/hist_cT")
plt.close()

plot_loghist(data[5], n_bins)
#plt.hist(data[5], n_bins)
plt.title(r'$c_S$')
plt.savefig("plots/hist_cS")
plt.close()

plot_loghist(data[6], n_bins)
#plt.hist(data[6], n_bins)
plt.title(r'$c_V$')
plt.savefig("plots/hist_cV")
plt.close()

plot_neg_loghist(data[7], n_bins, -8.)
plt.title(r'$\alpha_{ppE,0}^{(2,T)}$')
plt.savefig("plots/hist_alphaT")
plt.close()

plt.hist(data[8], n_bins)
plt.title(r'$g_{b1}$')
plt.savefig("plots/hist_gb1")
plt.close()

plt.hist(data[9], n_bins)
plt.title(r'$a_{bL}$')
plt.savefig("plots/hist_abL")
plt.close()

plot_neg_loghist(data[10], n_bins, -30.)
#bins = np.linspace(-10**-15,10**-15,n_bins)
#plt.hist(data[10], bins=bins)
plt.title(r'$g_{X1}$')
plt.savefig("plots/hist_gX1")
plt.close()

#plot_neg_loghist(data[11], half_bins)
plt.hist(data[11], n_bins)
plt.title(r'$\alpha_1$')
plt.savefig("plots/hist_alpha1")
plt.close()

plt.hist(data[12], n_bins)
plt.title(r'$\alpha_2$')
plt.savefig("plots/hist_alpha2")
plt.close()

"""
# Corner plots.
#cdata = np.hstack([data[0], data[1], data[2], data[3], data[4], data[5],
    #data[6], data[7], data[8], data[9], data[10]])
cdata = np.array([data[0], data[1], data[2], data[3], data[11], data[12],
    np.log10(data[5]), np.log10(data[6]), np.log10(abs(data[7])), data[8], data[9]])

figure = corner.corner(cdata.T, labels=[r'$c_a$', r'$c_\theta$', r'$c_\omega$',
    r'$c_\sigma$', r'$\alpha_1$', r'$\alpha_2$', r'log $c_S$', r'log $c_V$', r'log(|$\alpha_{ppE,0}^{(2,T)}$|)',
    r'$g_{b1}$', r'$a_{bL}$'], axes_scale = 'symlog',
    quantiles=[0.16, 0.5, 0.84], show_titles=True, title_kwargs={"fontsize": 12})
figure.savefig("plots/corner")
#figure.close()
"""
"""
# Plots for data drawn from a log distribution.
plot_loghist(data[0], n_bins)
plt.title(r'$c_a$')
plt.savefig("plots/hist_ca")
plt.close()

plot_symloghist(data[1], n_bins, 25, 25, -20.)
#plot_loghist(data[1], n_bins)
plt.title(r'$c_\theta$')
plt.savefig("plots/hist_ctheta")
plt.close()

plot_symloghist(data[2], n_bins, 40, 10, -20.)
#plot_symloghist(data[2], n_bins, 40, n_bins-40, -20.)
plt.title(r'$c_\omega$')
plt.savefig("plots/hist_cw")
plt.close()

plot_symloghist(data[3], n_bins, 25, 25, -20.)
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

plot_neg_loghist(data[7], n_bins, -8.)
plt.title(r'$\alpha_{ppE,0}^{(2,T)}$')
plt.savefig("plots/hist_alphaT")
plt.close()

plot_loghist(data[8], n_bins)
plt.title(r'$g_{b1}$')
plt.savefig("plots/hist_gb1")
plt.close()

#plt.hist(data[9], n_bins)
plot_symloghist(data[9], n_bins, 25, 25, -2.)
plt.title(r'$a_{bL}$')
plt.savefig("plots/hist_abL")
plt.close()

plot_symloghist(data[10], n_bins, 25, 25, -26.)
plt.title(r'$g_{X1}$')
plt.savefig("plots/hist_gX1")
plt.close()

plot_symloghist(data[11], n_bins, 15, 35, -22.)
plt.title(r'$\alpha_1$')
plt.savefig("plots/hist_alpha1")
plt.close()

plot_symloghist(data[12], n_bins, 15, 35, -22.)
plt.title(r'$\alpha_2$')
plt.savefig("plots/hist_alpha2")
plt.close()
"""

"""
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
"""
