import numpy as np
import matplotlib.pyplot as plt
import gwatpy.mcmc_routines as gmcmc

raw_data,raw_status, raw_model = gmcmc.RJPTMCMC_unpack_file("data/bayesline_output.hdf5")
#plt.hist(raw_data[:,0])
#plt.hist(raw_data[:,1])
knots = raw_data[:,2:41*2:2]
SNs = raw_data[:,3:41*2+1:2]
#knots = raw_data[:,2:4*2:2]
#SNs = raw_data[:,3:4*2+1:2]
knotsflat = np.array(knots.flat)
SNsflat = np.array(SNs.flat)
plt.hist(knotsflat,bins=500)
plt.savefig("plots/temp_hist.pdf")
plt.close()
plt.hist(SNsflat,bins=500)
plt.savefig("plots/temp_hist2.pdf")
plt.close()

#print(knotsflat.shape)
#knotsflat = np.array([ [x] for x in knotsflat])
#print(knotsflat.shape)
#plot_data = np.insert(knotsflat, 0, SNsflat, axis=1)
#print(plot_data.shape)

import seaborn as sns

#sns.histplot(plot_data)
sns.histplot(x=knotsflat, y=SNsflat)
plt.savefig("plots/temp_heatmap.pdf")
plt.close()


#from scipy.interpolate import interp1d
#freqs = np.linspace(100,300,1000)
##for x in np.arange(len(SNs[:,0])):
#for x in np.arange(10):
#    fn = interp1d(knots[x],SNs[x],kind='cubic')
#    plt.scatter(knots[x],SNs[x])
#    try:
#        SN = fn(freqs)
#        plt.plot(freqs,SN)
#    except:
#        continue
#plt.savefig("plots/reinterp.pdf")
#plt.close() 


#exit()



data = np.loadtxt("data/bayesline_raw_data.csv",delimiter=',').T
estimate = np.loadtxt("data/bayesline_estimate_psd.csv",delimiter=',')
print(estimate)
plt.semilogy(data[0], data[1]*data[1] + data[2]*data[2])
plt.semilogy(data[0], estimate)
#plt.savefig("plots/bayesline_raw_data.pdf")
#plt.close()
data = np.loadtxt("data/initial_position.csv",delimiter=',')
plot_data = data[3:4*2:2]
freqs = data[2:4*2-1:2]
#plt.semilogy(freqs,np.exp(plot_data))
plt.savefig("plots/bayesline_raw_data.pdf")
plt.close()
