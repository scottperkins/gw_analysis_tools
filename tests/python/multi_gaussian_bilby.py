import bilby 
import numpy as np

from bilby.core.likelihood import AnalyticalMultidimensionalCovariantGaussian,\
AnalyticalMultidimensionalBimodalCovariantGaussian

label = 'multivariate_gaussian_prior'
outdir = "data/bilby/"
bilby.utils.check_directory_exists_and_if_not_mkdir(outdir)

mean1 = np.loadtxt("data/multi_gaussian_like_mean.csv",delimiter=',')
dim = len(mean1)
mean2 = -1*mean1
covL = np.loadtxt("data/multi_gaussian_like_cov.csv",delimiter=',')
likelihood = AnalyticalMultidimensionalBimodalCovariantGaussian(mean1,mean2,covL)


#################################################
## PRIOR ##
#################################################

#mus = np.loadtxt("data/multi_gaussian_prior_mean.csv",delimiter=',')
#names = ["x{}".format(i) for i in np.arange(dim)]
#bounds = [(-10,10) for i in np.arange(dim)]
#cov = np.loadtxt("data/multi_gaussian_prior_cov.csv",delimiter=',')
#mvg = bilby.core.prior.MultivariateGaussianDist(names, nmodes=1, mus=mus, covs=cov,bounds=bounds)
#
#priors = dict()
#for i in np.arange(dim):
#    priors[names[i]] = bilby.core.prior.MultivariateGaussian(mvg,names[i]) 


priors = bilby.core.prior.PriorDict()
priors.update({"x{0}".format(i): bilby.core.prior.Uniform(-10, 10, "x{0}".format(i)) for i in range(dim)})
#################################################

if __name__ == "__main__":
    from multiprocessing import Pool
    #with Pool(8) as p:
    #    result = bilby.run_sampler(
    #    likelihood=likelihood, priors=priors, sampler='dynesty', nlive = 1e3,
    #    outdir=outdir, label=label,pool=p)

    result = bilby.run_sampler(
        likelihood=likelihood, priors=priors, sampler='ptemcee', nsamples = 1e3,nwalkers=20,threads=6,pos0='prior',ntemps=5,
        outdir=outdir, label=label)
