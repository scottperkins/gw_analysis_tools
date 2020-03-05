import numpy as np
import matplotlib.pyplot as plt

file1 = "currently_supported/CE2_strain.csv"
file2 = "currently_supported/CE2_strain_smoothed.csv"
data = np.loadtxt(file1, delimiter=',',unpack=True)
data_out_temp = data[1]
lag = 2
thresh_frac = .5
data_out = [ data_out_temp[i] if not ( (data_out_temp[i] - data_out_temp[i+lag]) >thresh_frac*data_out_temp[i+lag])else data_out_temp[i+lag] for i in np.arange(len(data_out_temp)-lag) ]
for x in data[1,-lag:]:
    data_out.append(x)
    
#data_out = [ data_out[i+1] for i in np.arange(len(data_out)-2) ]

plt.scatter(data[0],data_out,s=5)
plt.scatter(data[0],data[1],s=5,color="red")
plt.yscale('log')
plt.xscale('log')
plt.xlim([1,10000])
plt.ylim([1e-26,1e-18])
#plt.loglog(data[0],data[1])
plt.show()
plt.close()

output = np.asarray([data[0],data_out])
np.savetxt(file2,output.T, delimiter=',')
