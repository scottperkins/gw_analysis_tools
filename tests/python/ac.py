import gwatpy.mcmc_routines_ext as mcmc
#import emcee
import numpy as np
import matplotlib
#matplotlib.use('macosx')
import matplotlib.pyplot as plt
import emcee


#data =np.loadtxt("data/test_output.csv",delimiter=',')
data =np.loadtxt("data/injection_output.csv",delimiter=',')
#data =np.loadtxt("data/mcmc_output.csv",delimiter=',')
#data2 = data[:int(len(data)/2)]
#data = data[int(len(data)/2):]
dim = len(data[0])

#dataT =np.loadtxt("data/injection_output.csv",delimiter=',', unpack=True)
dataT =np.loadtxt("data/mcmc_output.csv",delimiter=',', unpack=True)
#dataT =np.loadtxt("data/test_output.csv",delimiter=',', unpack=True)
i = 0 
#for x in dataT:
#    #plt.plot(x[2000:])
#    plt.plot(x)
#    plt.savefig("plots/temp{}".format(i))
#    plt.close()
#    i+=1

for x in range(dim):
    data_thinned = []
    for y in np.arange(len(data)):
            if y%1 == 0:
                data_thinned.append(data[y][x])
    #plt.plot(data_thinned)
    #plt.show()
    #plt.close()
    print(emcee.autocorr.integrated_time(data_thinned, tol=50))
#for x in range(dim):
#    data_thinned = []
#    for y in np.arange(len(data2)):
#            if y%1 == 0:
#                data_thinned.append(data2[y][x])
#    #plt.plot(data_thinned)
#    #plt.show()
#    #plt.close()
#    print(emcee.autocorr.integrated_time(data_thinned, tol=1000))

data_thinned = []
for x in np.arange(len(data)):
        if x%1 == 0:
            data_thinned.append(data[x])
data_thinned = np.ascontiguousarray(data_thinned)
#data_thinned = data_thinned[100:]
mcmc.write_auto_corr_file_from_data_py(b'data/ac.csv',data_thinned,len(data_thinned),dim,5,0.01,10,True)
#data_thinned = []
#for x in np.arange(len(data2)):
#        if x%1 == 0:
#            data_thinned.append(data2[x])
#data_thinned = np.ascontiguousarray(data_thinned)
#mcmc.write_auto_corr_file_from_data_py(b'data/ac2.csv',data_thinned,len(data_thinned),dim,5,0.01,10,True)
