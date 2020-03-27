import numpy as np
import matplotlib.pyplot as plt
#import gwatpy.gwatpy_plot as gp ; gp.set()

data = np.loadtxt("data/noise_curves_terr.csv",delimiter=',')
dataS = np.loadtxt("data/noise_curves_space.csv",delimiter=',')
compS = np.loadtxt("../data/noise_data/NewLISATable.dat",unpack=True)
names = ["aLIGO_analytic","Hanford_O1_fitted","AdLIGOMidHigh","AdLIGOAPlus","AdLIGODesign","CE1","CE2","KAGRA_opt","KAGRA_pess","AdVIRGOPlus1","AdVIRGOPlus2_opt","AdVIRGOPlus2_pess","AdLIGODesign_smoothed","CE1_smoothed"]
#for x in np.arange(len(data)-1):
#    plt.loglog(data[0],data[x+1],label=names[x])
#plt.loglog(data[0],data[5],label=names[4])
#plt.loglog(data[0],data[-1],label=names[-1])
plt.loglog(data[0],data[11],label=names[10])
plt.loglog(data[0],data[12],label=names[11])


plt.legend()
plt.savefig("noise_curves_terr.pdf")
#plt.show()
plt.close()
for x in np.arange(len(dataS)-1):
    plt.loglog(dataS[0],dataS[x+1],label=str(x))
plt.loglog(compS[0],compS[2],linestyle='-.',label="Comp")
plt.legend()
#plt.show()
plt.legend()
plt.savefig("noise_curves_space.pdf")
plt.close()
