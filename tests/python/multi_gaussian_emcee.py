import emcee
import numpy as np
import matplotlib.pyplot as plt
import os
from corner import corner

os.environ["OMP_NUM_THREADS"] = "10"



def log_like(x, fisher, mean,scale):
    return1 = 0
    diff = (mean - x)
    return1 -= np.tensordot(np.tensordot(diff,fisher,axes=1),diff,axes=1)/2
    diff = (-mean - x)
    return2 = 0
    return2 -= np.tensordot(np.tensordot(diff,fisher,axes=1),diff,axes=1)/2
    #for l in np.arange(len(x)):
    #    for j in np.arange(len(x)):
    #        return1-=(mean[l] - x[l])*(mean[j] - x[j])/2 * fisher[l][j]
    #return2 = 0
    #for l in np.arange(len(x)):
    #    for j in np.arange(len(x)):
    #        return2-=(-mean[l] - x[l])*(-mean[j] - x[j])/2 * fisher[l][j]
    return np.log(np.exp(return1) + np.exp(return2) ) / scale
def log_p(x, fisher, mean,scale):
    return1 = 0
    #for l in np.arange(len(x)):
    #    for j in np.arange(len(x)):
    #        return1-=(mean[l] - x[l])*(mean[j] - x[j])/2 * fisher[l][j]
    #return return1/scale
    diff = (mean - x)
    return1 -= np.tensordot(np.tensordot(diff,fisher,axes=1),diff,axes=1)/2
    return return1/scale


def prob(x, fisherL,fisherP,meanL,meanP,scale):
    return log_like(x,fisherL,meanL,scale) + log_p(x,fisherP,meanP,scale)

from multiprocessing import Pool

with Pool() as pool:

    ndim = 5
    prior_mean = np.loadtxt("data/multi_gaussian_prior_mean.csv",delimiter=',')
    like_mean = np.loadtxt("data/multi_gaussian_like_mean.csv",delimiter=',')
    prior_cov = np.loadtxt("data/multi_gaussian_prior_cov.csv",delimiter=',')
    like_cov = np.loadtxt("data/multi_gaussian_like_cov.csv",delimiter=',')
    like_fish = np.linalg.inv(like_cov)
    prior_fish = np.linalg.inv(prior_cov)
    
    nwalkers = 200
    #p0 = np.random.rand(nwalkers, ndim)
    p0 = [-(1) * like_mean*np.random.rand(ndim) for x in np.arange(nwalkers)]
    #p0 = [  10*np.diagonal(like_cov)*np.random.rand(ndim) for x in np.arange(nwalkers)]
    scale = 1
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, prob, args=[like_fish,prior_fish, like_mean,prior_mean,scale],pool=pool)
    
    state = sampler.run_mcmc(p0, 10000)
    sampler.reset()
    sampler.run_mcmc(state, 40000,progress=True)
    tau = np.max(sampler.get_autocorr_time(tol=0))
    print("AC: ",tau)
    samples = sampler.get_chain(flat=True,thin=int(tau))
    np.savetxt("data/emcee_samples_multi_gaussian.csv",samples,delimiter=',')
    print("Samples: ",len(samples))
    fig = corner(samples)
    plt.savefig("plots/gaussian_emcee.pdf")
    plt.close()
    for x in np.arange(len(samples[0])):
        print(emcee.autocorr.integrated_time(samples[:,x]))
    
