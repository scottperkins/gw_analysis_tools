import numpy as np
import matplotlib.pyplot as plt
#import gwatpy.gwatpy_plot as gp ; gp.set()

def calculate_constants(GWAT,LAL,freqs):
    G1 = GWAT[0] 
    L1 = LAL[0] 
    F1 = freqs[0]
    G2 = GWAT[100] 
    L2 = LAL[100] 
    F2 = freqs[100]
    print(G1,G2,L1,L2,F1,F2)
    tc = (L2 - G2 - (L1 - G1))/ ( 2 *np.pi *(F2-F1))
    phic = L1 - G1 - 2 * np.pi * F1*tc
    return tc, phic
    #return 0, 0

for i in np.arange(1):
    data = np.loadtxt("data/response_{}.csv".format(i),delimiter=',',unpack=True)
    plt.loglog(data[0],data[1]*data[1]+data[2]*data[2],label="LAL")
    plt.loglog(data[0],data[3]*data[3]+data[4]*data[4],label="GWAT")
plt.legend()
#plt.xlim([50,5000])
plt.ylim([10**-55,10**-49])
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
    #phic, tc = calculate_constants(data[6], data[5],data[0])
    data[5] += 10
    data[6] += 10
    #data[6] += 2*np.pi*data[0] * tc
    #data[6] += phic
#    #plt.loglog(data[0],abs(data[1]-data[3])*2/(abs(data[1])+abs(data[3])))
#    #plt.loglog(data[0],abs(data[2]-data[4])*2/(abs(data[2])+abs(data[4])))
#    plt.loglog(data[0],abs(data[5]-data[6])*2/(abs(data[5])+abs(data[6])))
#    #plt.loglog(data[0],abs(data[1]-data[3]))
#    #plt.loglog(data[0],abs(data[2]-data[4]))
    #plt.loglog(data[0],abs(data[5]-data[6])/abs(data[5]))
    plt.loglog(data[0],abs(data[5]-data[6]))

    #plt.plot(data[0],abs(data[5]-data[6])/abs(data[5]))
    print(abs(data[5]-data[6])/abs(data[5]))
    #plt.plot(data[0],data[6])
#plt.show()

#plt.axvline(x=5000, color='r')
#plt.xlim([1600,1700])
plt.savefig("plots/diff_phase.pdf")
plt.close()
for i in np.arange(1):
    data = np.loadtxt("data/response_{}.csv".format(i),delimiter=',',unpack=True)
    
    #phic, tc = calculate_constants(data[6], data[5],data[0])
    data[5] += 10
    data[6] += 10
    #data[6] += 2*np.pi*data[0] * tc
    #data[6] += phic
#    #plt.loglog(data[0],abs(data[1]-data[3])*2/(abs(data[1])+abs(data[3])))
#    #plt.loglog(data[0],abs(data[2]-data[4])*2/(abs(data[2])+abs(data[4])))
#    plt.loglog(data[0],abs(data[5]-data[6])*2/(abs(data[5])+abs(data[6])))
#    #plt.loglog(data[0],abs(data[1]-data[3]))
#    #plt.loglog(data[0],abs(data[2]-data[4]))
    plt.plot(data[0],data[5],label="LAL")
    plt.plot(data[0],data[6],label="GWAT")
#plt.show()
plt.legend()
#plt.xlim([1600,1700])
plt.savefig("plots/comp_phase.pdf")
plt.close()

