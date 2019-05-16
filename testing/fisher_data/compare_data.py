import numpy as np
import matplotlib.pyplot as plt
from phenompy.gr import IMRPhenomD_detector_frame as imrdf
import seaborn as sns

data_num_in = np.loadtxt("numerical.csv",delimiter=',')
data_auto_in = np.loadtxt("autodiff.csv",delimiter=',')
data_py_in = np.loadtxt("python.csv")
param = np.loadtxt("parameters.csv",delimiter=',',unpack=True)
times = np.loadtxt("timing.csv",delimiter=',',unpack=True)

num_samples= len(data_num_in)
#num_samples=10
dimension = int(len(data_num_in[0])**(1/2))
dimsquared = len(data_num_in[0])

data_num = []
#for i in np.arange(len(num_samples):
#    data_num.append([])
#    for j in np.arange(dimension):
#        data_num[-1].append([])
#        for k in np.arange(dimension)):
#            data_num[-1][-1].append([])
#            data_num[-1][-1][-1] = data_num_in
#data_num =

c_diff = np.zeros((num_samples,dimsquared))
py_auto = np.zeros((num_samples,dimsquared))
py_num = np.zeros((num_samples,dimsquared))

for i in np.arange(num_samples): 
    for j in np.arange(dimsquared):
        c_diff[i][j] = ( data_num_in[i][j] - data_auto_in[i][j]) /data_auto_in[i][j]
        py_auto[i][j] = ( data_py_in[i][j] - data_auto_in[i][j]) /data_py_in[i][j]
        py_num[i][j] = ( data_num_in[i][j] - data_py_in[i][j]) /data_py_in[i][j]
#for i in np.arange(num_samples): 
c_diff_t = c_diff.transpose()
py_auto_t = py_auto.transpose()
py_num_t = py_num.transpose()
#i = 3
#print(py_auto_t[i*dimension +i][:])
#i = 4
#print(py_auto_t[i*dimension +i][:])
#i = 2
#print(py_auto_t[i*dimension +i][:])
bins = 20
names = [r'$ln A_0$',r'$t_c$',r'$\phi_c$',r'$ln \mathcal{M}$',r'$ln \eta$',r'$\chi_s$',r'$\chi_a$']
filename_base="fisher_"
filename_ext = ".pdf"
filedirectory = "plots/"
filenames = [filename_base+"{}{}".format(i,i)+filename_ext for i in np.arange(dimension)]
labels = ['C methods', 'py - auto', 'py - num']
for i in np.arange(dimension):
    #fig, axarr = plt.subplots(3)
    fig, axarr = plt.subplots(1)
    plt.suptitle(names[i])
    cplot = np.abs(c_diff_t[i*dimension +i][~np.isnan(c_diff_t[i*dimension +i])])
    py_auto_plot = np.abs(py_auto_t[i*dimension +i][~np.isnan(py_auto_t[i*dimension +i])])
    py_num_plot = np.abs(py_num_t[i*dimension +i][~np.isnan(py_num_t[i*dimension +i])])
    datarr = [cplot, py_auto_plot, py_num_plot]
    
    for j in np.arange(3):
        MIN, MAX = np.min(datarr[j]), np.max(datarr[j])
        if MIN==0:
            MIN= 1e-15
        #axarr[j].hist(datarr[j], label = labels[j],bins=np.exp(np.linspace(np.log(MIN),np.log(MAX),bins)))
        #axarr[j].set_xscale('log')
        #axarr[j].legend()
        #axarr.hist(datarr[j], label = labels[j],bins=np.exp(np.linspace(np.log(MIN),np.log(MAX),bins)))
        sns.distplot(datarr[j], label = labels[j],bins=np.exp(np.linspace(np.log(MIN),np.log(MAX),bins)), ax=axarr, hist=True,kde=False)
        axarr.set_xscale('log')
        axarr.legend()

    
    #sns.distplot(cplot, label = "auto num",ax=axarr[0],kde=False)
    #sns.distplot(py_auto_t[i*dimension +i], label = "py auto",ax=axarr[1],kde=False)
    #sns.distplot(py_num_t[i*dimension +i], label = "py num",ax=axarr[2],kde=False)
    
    #axarr[0].hist(cplot, label = "auto num",bins=np.exp(np.linspace(np.log(1e-15),np.log(1e-13),bins)))
    #axarr[0].set_xscale('log')
    #axarr[1].hist(py_auto_t[i*dimension +i], label = "py auto")
    #axarr[2].hist(py_num_t[i*dimension +i], label = "py num")
    #plt.show()
    plt.savefig(filedirectory+filenames[i])
    plt.close() 
