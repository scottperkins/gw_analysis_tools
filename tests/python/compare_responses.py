import numpy as np
import matplotlib.pyplot as plt
#import gwatpy.gwatpy_plot as gp ; gp.set()

for i in np.arange(1):
    data = np.loadtxt("data/response_{}.csv".format(i),delimiter=',',unpack=True)
    plt.loglog(data[0],data[1]*data[1]+data[2]*data[2],label="LAL")
    plt.loglog(data[0],data[3]*data[3]+data[4]*data[4],label="GWAT")
plt.legend()
plt.savefig("plots/comp_amp.pdf")
plt.close()
for i in np.arange(1):
    data = np.loadtxt("data/response_{}.csv".format(i),delimiter=',',unpack=True)
    plt.loglog(data[0],abs(data[1]*data[1]+data[2]*data[2]-data[3]*data[3]-data[4]*data[4])*2/(abs(data[1]*data[1]+data[2]*data[2])+abs(data[3]*data[3]+data[4]*data[4])))
#plt.show()
plt.savefig("plots/diff_amp.pdf")
plt.close()
for i in np.arange(1):
    data = np.loadtxt("data/response_{}.csv".format(i),delimiter=',',unpack=True)
    data[5] += 100
    data[6] += 100
#    #plt.loglog(data[0],abs(data[1]-data[3])*2/(abs(data[1])+abs(data[3])))
#    #plt.loglog(data[0],abs(data[2]-data[4])*2/(abs(data[2])+abs(data[4])))
#    plt.loglog(data[0],abs(data[5]-data[6])*2/(abs(data[5])+abs(data[6])))
#    #plt.loglog(data[0],abs(data[1]-data[3]))
#    #plt.loglog(data[0],abs(data[2]-data[4]))
    plt.plot(data[0],abs(data[5]-data[6])/abs(data[5]))
    #plt.plot(data[0],data[6])
#plt.show()
plt.savefig("plots/diff_phase.pdf")
plt.close()
for i in np.arange(1):
    data = np.loadtxt("data/response_{}.csv".format(i),delimiter=',',unpack=True)
#    #plt.loglog(data[0],abs(data[1]-data[3])*2/(abs(data[1])+abs(data[3])))
#    #plt.loglog(data[0],abs(data[2]-data[4])*2/(abs(data[2])+abs(data[4])))
#    plt.loglog(data[0],abs(data[5]-data[6])*2/(abs(data[5])+abs(data[6])))
#    #plt.loglog(data[0],abs(data[1]-data[3]))
#    #plt.loglog(data[0],abs(data[2]-data[4]))
    plt.plot(data[0],data[5],label="LAL")
    plt.plot(data[0],data[6],label="GWAT")
#plt.show()
plt.legend()
plt.savefig("plots/comp_phase.pdf")
plt.close()
