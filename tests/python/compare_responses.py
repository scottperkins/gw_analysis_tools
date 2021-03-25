import numpy as np
import matplotlib.pyplot as plt
#import gwatpy.gwatpy_plot as gp ; gp.set()

for i in np.arange(1):
    data = np.loadtxt("data/response_{}.csv".format(i),delimiter=',',unpack=True)
    plt.loglog(data[0],data[1]*data[1]+data[2]*data[2],label="LAL")
    plt.loglog(data[0],data[3]*data[3]+data[4]*data[4],label="GWAT")
plt.legend()
#plt.xlim([50,5000])
#plt.ylim([10**-35,10**-55])
plt.axvline(x=1291.9800951082, color='r')
plt.axvline(x=5882.29916274512, color='k')
plt.savefig("plots/comp_amp.pdf")
plt.close()
for i in np.arange(1):
    data = np.loadtxt("data/response_{}.csv".format(i),delimiter=',',unpack=True)
    plt.loglog(data[0],abs(data[1]*data[1]+data[2]*data[2]-data[3]*data[3]-data[4]*data[4])*2/(abs(data[1]*data[1]+data[2]*data[2])+abs(data[3]*data[3]+data[4]*data[4])))

#plt.show()
plt.axvline(x=1291.9800951082, color='r')
plt.axvline(x=5882.29916274512, color='k')

plt.savefig("plots/diff_amp.pdf")
plt.close()
for i in np.arange(1):
    data = np.loadtxt("data/response_{}.csv".format(i),delimiter=',',unpack=True)
    data[5] += 10
    data[6] += 10
#    #plt.loglog(data[0],abs(data[1]-data[3])*2/(abs(data[1])+abs(data[3])))
#    #plt.loglog(data[0],abs(data[2]-data[4])*2/(abs(data[2])+abs(data[4])))
#    plt.loglog(data[0],abs(data[5]-data[6])*2/(abs(data[5])+abs(data[6])))
#    #plt.loglog(data[0],abs(data[1]-data[3]))
#    #plt.loglog(data[0],abs(data[2]-data[4]))
    plt.semilogy(data[0],abs(data[5]-data[6])/abs(data[5]))

    #plt.plot(data[0],abs(data[5]-data[6])/abs(data[5]))
    print(abs(data[5]-data[6])/abs(data[5]))
    #plt.plot(data[0],data[6])
#plt.show()

plt.axvline(x=1661.11726513911, color='r')
plt.axvline(x=5882.29916274512, color='k')

#plt.xlim([50,1000])
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
plt.axvline(x=1661.11726513911, color='r')
plt.axvline(x=5882.29916274512, color='k')
plt.legend()
#plt.xlim([50,1000])
plt.savefig("plots/comp_phase.pdf")
plt.close()
