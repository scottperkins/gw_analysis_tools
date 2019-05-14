import numpy as np
import matplotlib.pyplot as plt
from phenompy.gr import IMRPhenomD_detector_frame as imrdf

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
i = 3
print(py_auto_t[i*dimension +i][:])
i = 4
print(py_auto_t[i*dimension +i][:])
i = 2
print(py_auto_t[i*dimension +i][:])
for i in np.arange(dimension):
    #for j in np.arange(dimension):
    plt.title(i)
    #plt.hist(c_diff_t[i*dimension +i][:],bins=20,density=True, label = "auto num")
    plt.hist(py_auto_t[i*dimension +i][:],bins=20,density=True, label = "py auto")
    #plt.hist(py_num_t[i*dimension +i][:],bins=20,density=True, label = "py num")
    plt.legend()
    plt.show()
    plt.close() 
#plt.hist(c_diff_t[0][:],bins=20,label="0")
#plt.hist(c_diff_t[2*dimension + 2][:],bins=20,label = "1")
#plt.hist(c_diff_t[3*dimension + 3][:],bins=20,normed=True)
#plt.hist(c_diff_t[int(dimension)+1][:],normed=True)
#plt.legend()
#plt.show()
#plt.close() 
