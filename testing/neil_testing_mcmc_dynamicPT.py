import numpy as np
import matplotlib.pyplot as plt
import corner 
burn = False
burnlength = 200000
sets =8

#data = np.loadtxt("data/neil_mcmc_output.csv",delimiter=',',unpack=True)
#burnin =1000
#burnin =0
#data = data[burnin:]
#plt.hist(data,bins=100,density=True)
#x = np.linspace(-3,3)
datasets = []
for i in np.arange(sets):
    datasets.append( np.loadtxt("data/neil_mcmc_output{}_dynamicPT.csv".format(i+1),delimiter=','))
data = [[],[]]
for i in datasets:
    for j in i:
        data[0].append(j[0])
        data[1].append(j[1])
     

plt.plot(data[0])
plt.plot(data[1])
plt.show()
plt.close()
datasets=[]
data = []
for i in np.arange(sets):
    datasets.append( np.loadtxt("data/neil_mcmc_output{}_dynamicPT.csv".format(i+1),delimiter=','))
for x in datasets:
    for y in x:
        data.append(y)
ndim, nsamples = 2, len(data) 
labels = [r"x",r"y"]
if burn:
    data = data[burnlength:]
#for x in data:
#    x[0] = abs(x[0])
#    x[1] = abs(x[1])
figure = corner.corner(data, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
plt.savefig("neil_mcmc_testing_dynamicPT.pdf")
plt.close()


#autocorr = np.loadtxt("data/neil_auto_corr_mcmc1_dynamicPT.csv",delimiter=',')
#lengths = autocorr[0]
#autocorr = autocorr[1:]
#for i in autocorr:
#    plt.plot(lengths,i)
#plt.savefig("neil_autocorr_testing_dynamicPT.pdf")
#plt.close()




#plt.hist(data,bins=100,density=True)
#x = np.linspace(-3,3)
#data = np.loadtxt("data/neil_mcmc_output2.csv",delimiter=',')
#newdata = []
#for x in data:
#    if x[1]!=0:
#        newdata.append([abs(x[0]*x[1]), abs(x[0]/x[1])])
#    #temp = x
#    #x[0] = (temp[1]*temp[0])
#    #x[1] = (temp[0]/temp[1])
#print(len(newdata))
#newdata = np.asarray(newdata)
#ndim, nsamples = 2, len(newdata) 
#labels = [r"x",r"y"]
#
#figure = corner.corner(newdata, labels=labels,quantiles=[.16,.5,.84], show_titles=True)
#plt.savefig("neil_mcmc_testing2.pdf")
#plt.close()


