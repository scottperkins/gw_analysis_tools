import numpy as np
import matplotlib.pyplot as plt

def signal_model( t, alphas,T):
    yti = np.sum([alphas[x] * (t/T)**x for x in np.arange(len(alphas))])
    return yti 

beta = 5
N = 100
dt = 1
M = 6

time = np.linspace(0,N*dt,N)
noise = np.random.normal(0, .1, N)
alphas = np.random.normal(0, beta,M)
signal = np.asarray([signal_model(t, alphas,N*dt) for t in time])
plt.plot(time,signal)
plt.plot(time,noise+signal)
plt.show()
plt.close()
np.savetxt("data/clean_data_transdimensional_{}_{}_{}_{}.csv".format(beta,int(M),dt,int(N)),signal,delimiter=',')
np.savetxt("data/full_data_transdimensional_{}_{}_{}_{}.csv".format(beta,int(M),dt,int(N)),signal+noise,delimiter=',')


