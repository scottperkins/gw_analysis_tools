import numpy as np
import emcee 
import corner 
import matplotlib.pyplot as plt
from multiprocessing import Pool
from time import time
import os 
os.environ["OMP_NUM_THREADS"] = "10"

distribution_data = np.loadtxt("data/rosenbock_parameters.csv")
a = distribution_data[0]
mu = distribution_data[1]
n = int(distribution_data[2])
n1 = int(distribution_data[3])
n2 = int(distribution_data[4])
b = distribution_data[5:]
print(distribution_data)

def log_l(x):
    LL = - a*(x[0] - mu)**2
    for j in np.arange(n2):
        for i in np.arange(1,n1):
            if(i == 1):
                LL -= b[j*(n1-1) +i] *(x[j*(n1-1) + i] - x[0]*x[0])**2 
            else:
                LL -= b[j*(n1-1) +i] *(x[j*(n1-1) + i] - x[j*(n1-1) + i-1]*x[j*(n1-1) + i-1])**2 
    return LL

#def log_l(x, mu, cov):
#    diff = x - mu
#    return -0.5 * np.dot(diff, np.linalg.solve(cov, diff))

start = time()
with Pool() as pool:
    #ndim, nwalkers = 5, 100
    ndim, nwalkers = n, 300
    ivar = None
    p0 = np.random.rand(nwalkers,ndim)
    #mu = np.random.rand(ndim)
    #cov = 0.5 - np.random.rand(ndim ** 2).reshape((ndim, ndim))
    #cov = np.triu(cov)
    #cov += cov.T - np.diag(cov.diagonal())
    #cov = np.dot(cov, cov)
    
    #sampler = emcee.EnsembleSampler(nwalkers, ndim, log_l, args=[mu,cov],pool=pool)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_l, args=[],pool=pool)
    state = sampler.run_mcmc(p0, 10000,progress=True)
    sampler.reset()
    print("Burned")
    sampler.run_mcmc(state,50000,progress=True)
    ac = sampler.get_autocorr_time(tol=0)
    print(ac)
print("finished sampling",time()-start)
samplesnonflat = sampler.get_chain(flat=False)
meantau=np.zeros(ndim)
print(samplesnonflat.shape)
for y in range(nwalkers):
    for x in range(ndim):
        #data_thinned = []
        #for y in np.arange(len(samples)):
        #        if y%1 == 0:
        #            data_thinned.append(samples[y][x])
        #plt.plot(data_thinned)
        #plt.show()
        #plt.close()
        #print(emcee.autocorr.integrated_time(data_thinned, tol=10))
        meantau[x]+=emcee.autocorr.integrated_time(samplesnonflat[:,y,x], tol=0)[0]
meantau/=nwalkers
max_meantau = int(np.amax(meantau))
print(max_meantau)
samples = sampler.get_chain(flat=True,thin=max_meantau)
np.savetxt("data/emcee_output.dat",samples)
fig = corner.corner(samples)
plt.savefig("plots/emcee_comparison.pdf")
plt.close()

