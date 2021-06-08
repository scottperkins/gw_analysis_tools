import bilby 
import numpy as np

from bilby.core.likelihood import AnalyticalMultidimensionalCovariantGaussian,\
AnalyticalMultidimensionalBimodalCovariantGaussian

label = 'multivariate_gaussian_prior'
outdir = "data/bilby/"
bilby.utils.check_directory_exists_and_if_not_mkdir(outdir)

mean1 = np.loadtxt("data/multi_gaussian_like_mean.csv",delimiter=',')
mean2 = -1*mean1
covL = np.loadtxt("data/multi_gaussian_like_cov.csv",delimiter=',')
likelihood = AnalyticalMultidimensionalBimodalCovariantGaussian(mean1,mean2,covL)



mus = np.loadtxt("data/multi_gaussian_prior_mean.csv",delimiter=',')
names = ["x{}".format(i) for i in np.arange(len(mus))]
bounds = [(-10,10) for i in np.arange(len(mus))]
cov = np.loadtxt("data/multi_gaussian_prior_cov.csv",delimiter=',')
mvg = bilby.core.prior.MultivariateGaussianDist(names, nmodes=1, mus=mus, covs=cov,bounds=bounds)

priors = dict()
for i in np.arange(len(names)):
    priors[names[i]] = bilby.core.prior.MultivariateGaussian(mvg,names[i]) 

#result = bilby.run_sampler(
#    likelihood=likelihood, priors=priors, sampler='dynesty', nlive=4000,
#    outdir=outdir, label=label)

#result = bilby.run_sampler(
#    likelihood=likelihood, priors=priors, sampler='ptemcee', nsamples = 1e3,nwalkers=200,threads=8,pos0='prior',ntemps=10,
#    outdir=outdir, label=label)

if __name__ == "__main__":
    from multiprocessing import Pool
    #with Pool(8) as p:
    #    result = bilby.run_sampler(
    #    likelihood=likelihood, priors=priors, sampler='dynesty', nlive = 1e3,
    #    outdir=outdir, label=label,pool=p)

    result = bilby.run_sampler(
        likelihood=likelihood, priors=priors, sampler='ptemcee', nsamples = 1e3,nwalkers=200,threads=8,pos0='prior',ntemps=10,
        outdir=outdir, label=label)
