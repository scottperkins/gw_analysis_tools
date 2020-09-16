import numpy as np
import matplotlib.pyplot as plt
import gwatpy.util as gpu
import gwatpy.detector_util as gpdu
import gwatpy.mcmc_routines  as gmcmc
import gwatpy.waveform_generator as gwg
from time import time


from scipy.stats import binom, dirichlet
from scipy.optimize import curve_fit,minimize, NonlinearConstraint,differential_evolution
import scipy
print(scipy.__version__)

#pvec = [.3,.3,.4]
means = np.random.uniform(size=10)*0
sigmas = np.random.uniform(size=10)*10
data = [np.random.normal(means[x],sigmas[x], 10000) for x in np.arange(len(means))]

bins =100
#bins = np.linspace(minval, maxval, 7)
#mvec, bins = np.histogram(data1, bins=bins)
#mvec2, bins = np.histogram(data2, bins=bins)
#
##mvec = [400,100,140,100,30,10]
##mvec2 = [500,400,240,20,300,10]
#
##mvec /= np.sum(mvec)
#mvecnorm = mvec/np.sum(mvec)
#mvecnorm2 = mvec2/np.sum(mvec2)
#mvecprod = mvecnorm*mvecnorm2/np.sum(mvecnorm*mvecnorm2)
#print(mvec)
#print("vec 1 norm",mvecnorm)
#print("vec 2 norm",mvecnorm2)
#print("product",mvecprod)
#xvals= np.linspace(0,10,len(mvec))
#
##print(dirichlet.pdf(pvec, mvec))
#
#
#def dirichlet_wrapper_test(p):
#    ptemp = np.insert(np.asarray(p), len(p), 1-np.sum(p))
#    #print("ptemp: ",ptemp)
#    return -np.exp(dirichlet.logpdf(ptemp, mvec)+dirichlet.logpdf(ptemp,mvec2))
#
#
##constr = ({'type':'ineq','fun':lambda x: 1-np.sum(x)})
#constr = NonlinearConstraint(lambda x : 1-np.sum(x), 0, 1)
##bnds = ((0,1),(0,1))
#bnds = ((0,1),(0,1),(0,1),(0,1),(0,1))
##bnds = [(0,0),(1,1)]
##print(dirichlet_wrapper(np.ones(len(mvec))/len(mvec)))
##res = minimize(dirichlet_wrapper, np.ones(len(mvec))/(len(mvec)), method='SLSQP',bounds=bnds,constraints=constr)
##res = minimize(dirichlet_wrapper, np.ones(len(mvec)-1)/(len(mvec)), method='SLSQP',bounds=bnds,constraints=constr,options={'maxiter':1000,'disp':True})
#
#res = differential_evolution(dirichlet_wrapper_test, bounds=bnds,constraints=(constr))
#pvec_fit = np.insert(res.x,len(res.x), 1-np.sum(res.x))
#print(res)
#print(pvec_fit)
#plt.plot(bins[1:],pvec_fit,label='dirichlet')
#plt.plot(bins[1:],mvecnorm,label='1')
#plt.plot(bins[1:],mvecnorm2,label='2')
#plt.plot(bins[1:],mvecprod,label='prod')
#plt.hist(data1,bins=bins,density=True,alpha=.5)
#plt.hist(data2,bins=bins,density=True,alpha=.5)


pvec_test, pvec_arr,bins_test = gmcmc.combine_discrete_marginalized_posteriors(data,bins=bins)
plt.plot(bins_test,pvec_test,label='dirichlet test')
#bins_test = np.linspace(np.amin(data),np.amax(data), bins)

for i in np.arange(len(means)):
    plt.hist(data[i],bins=bins_test,alpha=.5,density=True,label='data {}'.format(i))

#prod = np.ones(bins)
#for i in np.arange(len(means)):
#    hist, b = np.histogram(data[i],bins=bins_test)
#    prod *= hist
#prod /=  (np.sum(prod)*(b[1]-b[0]))
#print(np.sum(prod))
#plt.plot(b[1:],prod,label='prod')

print(sigmas)
sigma_tot = np.sqrt( 1./ np.sum( 1/sigmas**2))
print(sigma_tot)
#xvals = np.linspace(np.amin(data),np.amax(data),1000)
xvals = bins_test
vals = 1./np.sqrt( sigma_tot**2 * 2 *np.pi) * np.exp( - ( xvals)**2/(2.*sigma_tot**2))
print(vals)
plt.plot(xvals,vals,label='true')

#true_dist = 1./np.sqrt( sigma1**2 * 2 *np.pi) * np.exp( - ( bins_test[1:]-mean1)**2/(2.*sigma1**2))*\
#    1./np.sqrt( sigma2**2 * 2 *np.pi) * np.exp( - ( mean2-bins_test[1:])**2/(2.*sigma2**2))
#plt.plot(bins_test[1:],true_dist,label='true')


plt.legend()
plt.show()

#print( dirichlet.pdf(pvec_fit, mvec))
#print( multinom.pdf(pvec_fit, mvec))
#print(np.insert(res, len(res),1-np.sum(res)))


#xcts = [100,120,100]
#ps = [.3,.4,.3]
#rv = binom(xcts,np.sum(xcts), ps)
#plt.plot(x,rv.pmf(




#LAT, LONG, LOC, D = gpdu.get_detector_parameters_py("Hanford")
#print(LAT,LONG,LOC,D)
#for x in np.arange(100):
#    kwargs = {"mass1":11,"mass2":9,"Luminosity_Distance":500,"spin1":[.1,.8,.1],"Nmod":1,"bppe":[-1],"betappe":[10]}
#    gp = gpu.gen_params(**kwargs)
#    freq = np.linspace(1,100,1000)
#    gen_meth = "ppE_IMRPhenomPv2_IMR"
#    start = time()
#    wfp = gwg.response_generator(freq,"Hanford", gen_meth, gp)
#    
#    kwargs = {"mass1":11,"mass2":9,"Luminosity_Distance":500,"spin1":[.1,.8,.1],"Nmod":1,"bppe":[-1],"betappe":[0]}
#    gp2 = gpu.gen_params(**kwargs)
#    wfp2 = gwg.response_generator(freq,"Hanford", gen_meth, gp2)
#    end = time()
#    print(end - start)
#
#
#plt.semilogx(np.imag(wfp))
#plt.semilogx(np.imag(wfp2),linestyle=':')
#plt.show()
#plt.close()


