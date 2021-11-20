import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

def credible_region(data_1D,CR_interval=.68,npoints=1e4):
    #Create KDE of data for easy evaluation of PDF
    kde = gaussian_kde(data_1D)

    #Calculate the grid on which we will search
    grid = np.linspace(np.amin(data_1D),np.amax(data_1D),num =int( npoints) )

    #Create array of probablitities at each grid point
    #(actually, ``p(x)dx'' so its properly normalized for the integral)
    diff = (grid[1:]-grid[:-1])
    new_grid = grid[1:]

    probs = kde.pdf(new_grid)*diff#*norm

    #Sort the array -- this creats an array of INDICES that is sorted.
    #We need to sort the grid in the exact same way as the probabilities
    #-- flip reverses the array because argsort only does ascending, we want descending
    indices_prob_sorted = np.flip(probs.argsort())

    #create sorted arrays for probability and grid points, so probs_sorted[index]
    #corresponds to grid_sorted[index]
    probs_sorted = probs[indices_prob_sorted]
    #print("First ten probability: ",probs_sorted[:10])
    grid_sorted = new_grid[indices_prob_sorted]

    #calculate the cumulative sum --
    #cumulative_sum[index] = Sum_{i}^{index} probs_sorted[i]
    cumulative_sum = np.cumsum(probs_sorted)
    #print("cumulative_sum: ", cumulative_sum)

    #Determine where on the grid the interval crosses the desired region
    #(say 90% or 68%)
    CR_prob_index = np.argwhere(cumulative_sum >= CR_interval) [0,0]

    #for the array of values that contributed to the integral, calculate the
    #min/max x values
    min_val = np.amin(grid_sorted[:CR_prob_index])
    max_val = np.amax(grid_sorted[:CR_prob_index])

    kde_approx =kde.integrate_box_1d(min_val,max_val)
    #print("Relative error (relative to KDE integration from min_val to max_val -- INCORRECT FOR DOUBLY PEAKED DISTRIBUTIONS)", (kde_approx - cumulative_sum[CR_prob_index])/kde_approx)
    return min_val,max_val

def plot_loghist(x, bins):
    hist, bins = np.histogram(x, bins=bins)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(x, bins=logbins)#,density=True)
    #plt.hist(x, weights=np.ones(len(x))/len(x),bins=logbins,density=False)
    plt.xscale('log')

def plot_neg_loghist(x, bins, upper_bound):
    hist, bins = np.histogram(x, bins=bins)
    logbins = np.flip(-np.logspace(upper_bound, np.log10(abs(bins[0])), len(bins)))
    plt.hist(x, bins=logbins)#,density=True)
    plt.xscale('symlog', base=10, linthresh=10**upper_bound)

def plot_symloghist(x, bins, posbins, negbins, bound):
    hist, bins = np.histogram(x, bins=bins)
    pos_logbins = np.logspace(bound, np.log10(bins[-1]), posbins)
    neg_logbins = np.flip(- np.logspace(bound, np.log10(abs(bins[0])), negbins))
    logbins = np.concatenate((neg_logbins, pos_logbins), axis=None)
    plt.hist(x, bins=logbins)#,density=True)
    plt.xscale('symlog', base=10, linthresh=10**bound, linscale=2)

# Loading in actual data.
iterations = 13
n_bins = 50
for i in np.arange(iterations):
    #data = np.loadtxt("data/EA_parameter_MC.csv".format(i),delimiter=',',unpack=True)
    data = np.loadtxt("data/uniform/case1/EA_parameter_MC3.csv".format(i), delimiter=',', unpack=True)
    #data = np.loadtxt("data/uniform/case2/EA_parameter_MC.csv".format(i), delimiter=',', unpack=True)
    #data = np.loadtxt("data/log/case1/EA_parameter_MC2.csv".format(i), delimiter=',', unpack=True)
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 20}
# Calculate credible regions
minca,maxca = credible_region(data[0])
print("Credible region for ca: ({}, {})\n".format(minca, maxca))
plt.hist(data[0], n_bins, density=True)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax = plt.gca()
ax.yaxis.offsetText.set_fontsize(15)
ax.xaxis.offsetText.set_fontsize(15)
#plot_loghist(data[0], n_bins)
plt.axvline(minca,color='r', label="68% credible region \n ({:0.2},{:0.2})".format(minca, maxca))
plt.axvline(maxca,color='r')

plt.legend(loc="upper right", fontsize=12)
plt.title(r'$c_a$', fontsize=20)
fig = plt.gcf()
fig.set_size_inches(5, 3)
plt.savefig("plots/credible_region_ca", bbox_inches="tight")
plt.close()


minctheta,maxctheta = credible_region(data[1])
print("Credible region for ctheta: ({}, {})\n".format(minctheta, maxctheta))
plt.hist(data[1], n_bins, density=True)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax = plt.gca()
ax.yaxis.offsetText.set_fontsize(15)
ax.xaxis.offsetText.set_fontsize(15)
#plot_loghist(data[1], n_bins)
plt.axvline(minctheta,color='r', label="68% credible region \n ({:0.2},{:0.2})".format(minctheta, maxctheta))
plt.axvline(maxctheta,color='r')
plt.legend(loc="upper right", fontsize=12)
plt.title(r'$c_\theta$', fontsize=20)
fig = plt.gcf()
fig.set_size_inches(5, 3)
plt.savefig("plots/credible_region_ctheta", bbox_inches="tight")
plt.close()


mincw,maxcw = credible_region(data[2])
print("Credible region for cw: ({}, {})\n".format(mincw, maxcw))
plt.hist(data[2], n_bins, density=True)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax = plt.gca()
ax.yaxis.offsetText.set_fontsize(15)
ax.xaxis.offsetText.set_fontsize(15)
#plt.axvline(mincw,color='r')
#plt.axvline(maxcw,color='r')
plt.title(r'$c_\omega$', fontsize=20)
fig = plt.gcf()
fig.set_size_inches(5, 3)
plt.savefig("plots/credible_region_cw", bbox_inches="tight")
plt.close()


mincsigma, maxcsigma = credible_region(data[3])
print("Credible region for csigma:({},{})\n".format(mincsigma, maxcsigma))
plt.hist(data[3], n_bins, density=True)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax = plt.gca()
ax.yaxis.offsetText.set_fontsize(15)
ax.xaxis.set_major_locator(plt.MaxNLocator(5))
ax.xaxis.offsetText.set_fontsize(15)
#plot_symloghist(data[3], n_bins, 25, 25, -20.)
#plt.axvline(mincsigma,color='r')
#plt.axvline(maxcsigma,color='r')
plt.title(r'$c_\sigma$', fontsize=20)
fig = plt.gcf()
fig.set_size_inches(5, 3)
plt.savefig("plots/credible_region_csigma", bbox_inches="tight")
plt.close()
