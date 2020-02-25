import gwatpy.mcmc_routines_ext as mcmc
#import emcee
import numpy as np
import matplotlib
#matplotlib.use('macosx')
import matplotlib.pyplot as plt

#data =np.loadtxt("data/mcmc_pv2_chain.csv",unpack=True,delimiter=',')
#for i in np.arange(len(data)):
#    print(emcee.autocorr.integrated_time(x=data[i],c=5,tol=50,quiet=False))

#data =np.loadtxt("data/mcmc_output_uncorr_P.csv",delimiter=',')
#data =np.loadtxt("data/mcmc_output_uncorr_D.csv",delimiter=',')
#data =np.loadtxt("data/test.csv",delimiter=',')
data =np.loadtxt("data/test2.csv",delimiter=',')
#data =np.loadtxt("data/mcmc_output_uncorr_Pv2_in.csv",delimiter=',')
dim = len(data[0])
#data =np.loadtxt("data/neil_mcmc_output_uncorr.csv",delimiter=',')
#data=data[2000:]

#dataT =np.loadtxt("data/mcmc_output_uncorr_P.csv",delimiter=',', unpack=True)
#dataT =np.loadtxt("data/mcmc_output_uncorr_D.csv",delimiter=',', unpack=True)
dataT =np.loadtxt("data/test2.csv",delimiter=',', unpack=True)
#dataT =np.loadtxt("data/mcmc_output_uncorr_Pv2_in.csv",delimiter=',', unpack=True)
#dataT =np.loadtxt("data/test.csv",delimiter=',', unpack=True)
i = 0 
for x in dataT:
    #plt.plot(x[2000:])
    plt.plot(x)
    plt.savefig("temp{}".format(i))
    plt.close()
    i+=1

#data = data[10000:]
data_thinned = []
for x in np.arange(len(data)):
        if x%1 == 0:
            data_thinned.append(data[x])
data_thinned = np.ascontiguousarray(data_thinned)
#data_thinned = data_thinned[100:]
mcmc.write_auto_corr_file_from_data_py(b'data/ac.csv',data_thinned,len(data_thinned),dim,10,0.01,10,True)
#mcmc.write_auto_corr_file_from_data_file_py(b'data/ac.csv',b'data/mcmc_pv2_chain.csv',200000,12,10,0.01,10,True)


#data_thinned2 = []
#for x in np.arange(len(data_thinned)):
#        if x%10 == 0:
#            data_thinned2.append(data_thinned[x])
#data_thinned2 = np.ascontiguousarray(data_thinned2)
##data = data[10000:]
#mcmc.write_auto_corr_file_from_data_py(b'data/ac2.csv',data_thinned2,len(data_thinned2),dim,10,0.01,10,True)
##mcmc.write_auto_corr_file_from_data_file_py(b'data/ac.csv',b'data/mcmc_pv2_chain.csv',200000,12,10,0.01,10,True)
#
#data_thinned3 = []
#for x in np.arange(len(data_thinned2)):
#        if x%10 == 0:
#            data_thinned3.append(data_thinned2[x])
#data_thinned3 = np.ascontiguousarray(data_thinned3)
##data = data[10000:]
#mcmc.write_auto_corr_file_from_data_py(b'data/ac3.csv',data_thinned3,len(data_thinned3),dim,10,0.01,10,True)
##mcmc.write_auto_corr_file_from_data_file_py(b'data/ac.csv',b'data/mcmc_pv2_chain.csv',200000,12,10,0.01,10,True)
#
