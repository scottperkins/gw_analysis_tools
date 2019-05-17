import numpy as np
import matplotlib.pyplot as plt
import corner 

#data = np.loadtxt("data/neil_mcmc_output.csv",delimiter=',',unpack=True)
#burnin =1000
#burnin =0
#data = data[burnin:]
#plt.plot(data[0][burnin:])
#plt.plot(data[1][burnin:])
#plt.show()
#plt.close()
#plt.hist(data,bins=100,density=True)
#x = np.linspace(-3,3)
data = np.loadtxt("data/neil_mcmc_output.csv",delimiter=',')
ndim, nsamples = 2, len(data) 
labels = [r"x",r"y"]

#for x in data:
#    x[0] = abs(x[0])
#    x[1] = abs(x[1])
figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("neil_mcmc_testing.pdf")
plt.close()


autocorr = np.loadtxt("data/neil_auto_corr_mcmc.csv",delimiter=',')
lengths = autocorr[0]
autocorr = autocorr[1:]
for i in autocorr:
    plt.plot(lengths,i)
plt.savefig("neil_autocorr_testing.pdf")
plt.close()




#plt.hist(data,bins=100,density=True)
#x = np.linspace(-3,3)
data = np.loadtxt("data/neil_mcmc_output2.csv",delimiter=',')
newdata = []
for x in data:
    if x[1]!=0:
        newdata.append([abs(x[0]*x[1]), abs(x[0]/x[1])])
    #temp = x
    #x[0] = (temp[1]*temp[0])
    #x[1] = (temp[0]/temp[1])
print(len(newdata))
newdata = np.asarray(newdata)
ndim, nsamples = 2, len(newdata) 
labels = [r"x",r"y"]

figure = corner.corner(newdata, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("neil_mcmc_testing2.pdf")
plt.close()


