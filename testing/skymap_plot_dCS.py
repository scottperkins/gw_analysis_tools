import matplotlib as mpl
mpl.use('pdf')
import astropy.coordinates as coord
import numpy as np
import matplotlib.pyplot as plt
burn=False
data = np.loadtxt("data/mcmc_output_dCS.csv",usecols=(1,2),delimiter=',',unpack=True)
if burn:
    data = data[10000:]
#mollweide projection uses RA \el [-pi, pi], not [0, 2 pi]
for x in np.arange(len(data[0])):
    if data[0][x] > np.pi:
        data[0][x]-=2*np.pi
    #data[0][x]-=np.pi
xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']

#for x in np.arange(len(data[0])):
#    if data[0][x] > np.pi:
#        data[0]-=2*np.pi
#    #data[0][x]-=np.pi
#data = [[],[]]
#data[0] = [4*np.pi/3. - 2*np.pi]
#data[1] = [np.pi/4.]
#xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']


fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111,projection="mollweide")
ax.set_xticklabels(xlab, weight=800)
ax.scatter(data[0],data[1], s=.1)
plt.grid(True)
plt.savefig("skymap_dCS.pdf")

plt.close()
