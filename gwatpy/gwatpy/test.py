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
N = 3
means = np.random.uniform(size=N)*3-1.5
sigmas = np.random.uniform(size=N)*5
cts = np.random.uniform(size=N)*1000+1000
#sigmas = np.ones(10)
data = np.asarray([np.random.normal(means[x],sigmas[x], int(cts[x])) for x in np.arange(len(means))])

bins =20
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


pvec_test,bins_mid,bins_edges = gmcmc.combine_discrete_marginalized_posteriors(data,bins=bins)
#plt.plot(bins_mid,pvec_test,label='dirichlet test')
#bins_test = np.linspace(np.amin(data),np.amax(data), bins)

for i in np.arange(len(means)):
    plt.hist(data[i],bins=bins_edges,alpha=.5,density=True,label='data {}'.format(i))

prod = np.ones(bins-1)
dists = []
for i in np.arange(len(means)):
    hist, b = np.histogram(data[i],bins=bins_edges)
    prod *= hist
    dists.append(scipy.stats.rv_histogram((hist,b)))
prod /=  (np.sum(prod)*(b[1]-b[0]))
print("PVEC",pvec_test)
print("PROD",prod)
print(np.sum(prod))
#plt.plot(bins_mid,prod,label='prod')
fits = np.product([fn.pdf(bins_mid) for fn in dists], axis=0)
fits /= (np.sum(fits)*(bins_mid[1]-bins_mid[0]))
#plt.plot(bins_mid, fits,label='fits')

print(sigmas)
sigma_tot = np.sqrt( 1./ np.sum( 1/sigmas**2))
mean_tot = np.sum(means/(2*sigmas**2)) / (np.sum(1./(2*sigmas**2)))
print(sigma_tot)
#xvals = np.linspace(np.amin(data),np.amax(data),1000)
xvals = bins_mid
def true_dist(x,sigmatot,meantot):
    return 1./np.sqrt( sigmatot**2 * 2 *np.pi) * np.exp( - (meantot- x)**2/(2.*sigmatot**2))
vals = true_dist(xvals, sigma_tot,mean_tot)
#plt.plot(xvals,vals,label='true')
xvals = np.linspace(np.amin(bins_mid),np.amax(bins_mid),1000)
vals = true_dist(xvals, sigma_tot,mean_tot)
plt.plot(xvals,vals,label='true-full')
print("True sigma",sigma_tot)
print("True mean",mean_tot)


dirichlet_fn = scipy.interpolate.interp1d(bins_mid, pvec_test,kind='linear')
prod_fn = scipy.interpolate.interp1d(bins_mid, prod,kind='linear')
per_dirichlet_upper = scipy.optimize.fsolve(lambda x: 0.05-scipy.integrate.quad(dirichlet_fn,x,xvals[-1])[0], x0=0)
per_prod_upper = scipy.optimize.fsolve(lambda x: 0.05-scipy.integrate.quad(prod_fn,x,xvals[-1])[0], x0=0)
per_true_upper = scipy.optimize.fsolve(lambda x: 0.05-scipy.integrate.quad(lambda y: true_dist(y,sigma_tot,mean_tot),x,xvals[-1])[0], x0=0)
per_dirichlet_lower = scipy.optimize.fsolve(lambda x: 0.05-scipy.integrate.quad(dirichlet_fn,xvals[0],x)[0], x0=0)
per_prod_lower = scipy.optimize.fsolve(lambda x: 0.05-scipy.integrate.quad(prod_fn,xvals[0],x)[0], x0=0)
per_true_lower = scipy.optimize.fsolve(lambda x: 0.05-scipy.integrate.quad(lambda y: true_dist(y,sigma_tot,mean_tot),xvals[0],x)[0], x0=0)
print("Dirichlet 90%",per_dirichlet_upper)
print("Prod 90%",per_prod_upper)
print("True 90%",per_true_upper)
print("Dirichlet 90%",per_dirichlet_lower)
print("Prod 90%",per_prod_lower)
print("True 90%",per_true_lower)

plt.plot(xvals, dirichlet_fn(xvals),label='dirichlet')
plt.plot(xvals, prod_fn(xvals),label='prod')

#print("95 dirichlet",per_dirichlet)
#print("95 dirichlet",per_dirichlet)


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


