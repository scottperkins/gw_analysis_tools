import gwatpy.mcmc_routines_ext as mcmc
#import emcee
import numpy as np
import matplotlib
#matplotlib.use('macosx')
import matplotlib.pyplot as plt
import emcee

#data =np.loadtxt("data/mcmc_output.csv",delimiter=',',unpack=True)
#data = np.sin(np.linspace(1,10,50))
#datamid = np.asarray([[x] for x in data])
#dim = 1
#datain = np.ascontiguousarray(datamid)
#mcmc.write_auto_corr_file_from_data_py(b'data/ac.csv',datain,len(datain),dim,1,0.01,10,True)
def next_pow_two(n):
    i = 1
    while i < n:
        i = i << 1
    return i


def autocorr_func_1d(x, norm=True):
    x = np.atleast_1d(x)
    if len(x.shape) != 1:
        raise ValueError("invalid dimensions for 1D autocorrelation function")
    n = next_pow_two(len(x))
    # Compute the FFT and then (from that) the auto-correlation function
    f = np.fft.fft(x - np.mean(x), n=2 * n)
    acf = np.fft.ifft(f * np.conjugate(f))[: len(x)].real
    acf /= 4 * n

    # Optionally normalize
    if norm:
        acf /= acf[0]

    return acf
def auto_window(taus, c):
    m = np.arange(len(taus)) < c * taus
    if np.any(m):
        return np.argmin(m)
    return len(taus) - 1


# Following the suggestion from Goodman & Weare (2010)
def autocorr_gw2010(y, c=5.0):
    f = autocorr_func_1d(np.mean(y, axis=0))
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]


def autocorr_new(y, c=5.0):
    f = np.zeros(y.shape[1])
    for yy in y:
        f += autocorr_func_1d(yy)
    f /= len(y)
    taus = 2.0 * np.cumsum(f) - 1.0
    window = auto_window(taus, c)
    return taus[window]


data =np.loadtxt("data/mcmc_output.csv",delimiter=',',unpack=True)
dataTest = np.asarray([data[0]])
#dataTest = np.asarray([data.T])
c = 5
#print( autocorr_new(dataTest))
#print(autocorr_func_1d(dataTest[0]))
taus = 2.0 * np.cumsum ( autocorr_func_1d(dataTest[0])) -1.0
window = auto_window(taus,c)
#print(window)
#print( taus[window])

#exit()
data =np.loadtxt("data/mcmc_output.csv",delimiter=',')
dim = len(data[0])
for x in range(dim):
    data_thinned = []
    for y in np.arange(len(data)):
            if y%1 == 0:
                data_thinned.append(data[y][x])
    print(emcee.autocorr.integrated_time(data_thinned, tol=10))
#exit()
data_thinned = []
for x in np.arange(len(data)):
        if x%1 == 0:
            data_thinned.append(data[x])
data_thinned = np.ascontiguousarray(data_thinned)
mcmc.write_auto_corr_file_from_data_py(b'data/ac.csv',data_thinned,len(data_thinned),dim,5,0.01,10,True)
