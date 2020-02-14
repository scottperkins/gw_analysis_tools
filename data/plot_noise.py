import numpy as np
import matplotlib.pyplot as plt
import gwatpy.gwatpy_plot as gp ; gp.set()
gp.show_minor_grid()


ELIGO = np.loadtxt("noisecurves/AdLIGO_allnoises.dat",unpack=True)
ELIGOAplus = np.loadtxt("noisecurves/Aplus.dat",unpack=True)
ECE1 = np.loadtxt("noisecurves/CE1.dat",unpack=True)
ECE2N = np.loadtxt("noisecurves/CE2narrow.dat",unpack=True)
ECE2W = np.loadtxt("noisecurves/CE2wide.dat",unpack=True)
EAdV = np.loadtxt("noisecurves/AdV_baseline_sensitivity_12May09.csv",delimiter=',',unpack=True)
EAdV2 = np.loadtxt("noisecurves/AdV_refsens_090427.csv",delimiter=',',unpack=True)
AdLIGO_design = np.loadtxt("local_noise_curves/aLIGODesign.csv",delimiter=',',unpack=True)
AdLIGO_Aplus = np.loadtxt("local_noise_curves/AplusDesign.csv",delimiter=',',unpack=True)
CE1 = np.loadtxt("local_noise_curves/CE1_strain.csv",delimiter=',',unpack=True)
CE2 = np.loadtxt("local_noise_curves/CE2_strain.csv",delimiter=',',unpack=True)
bK = np.loadtxt("local_noise_curves/kagra_DRSE.csv",delimiter=',',unpack=True)


fig, ax = plt.subplots(nrows=3,ncols=2,sharex=True,figsize=(15,15))
ax[0,0].loglog(AdLIGO_design[0],AdLIGO_design[1],label="T1800044-v5")
ax[0,0].loglog(ELIGO[0],ELIGO[4],label="Emanuele_AdLIGO_Design")

ax[0,1].loglog(AdLIGO_Aplus[0],AdLIGO_Aplus[1],label="T1800042-v4")
ax[0,1].loglog(ELIGOAplus[0],ELIGOAplus[1],label="Emanuele_AdLIGO_Aplus")

ax[2,0].loglog(EAdV[0],EAdV[1],label="Emanuele AdV phase 2")
ax[2,0].loglog(EAdV[0],EAdV[1]*1.71034,label="Emanuele AdV phase 1 analog")
#ax[2,0].loglog(EAdV2[0],EAdV2[1],label="Emanuele AdV Refnes")

ax[1,0].loglog(CE1[0],CE1[1],label="CE1 cosmicexplorer.org")
ax[1,0].loglog(ECE1[0],ECE1[1],label="Emanuele CE1")

ax[1,1].loglog(CE2[0],CE2[1],label="Scott cosmicexplorer.org")
ax[1,1].loglog(ECE2N[0],ECE2N[1],label="Emanuele CE2N")
ax[1,1].loglog(ECE2W[0],ECE2W[1],label="Emanuele CE2W")

ax[2,1].loglog(bK[0],bK[1],label="bKAGRA https://gwcenter.icrr.u-tokyo.ac.jp/ O4")
ax[2,1].loglog(bK[0],bK[1]*3.9326,label="bKAGRA O3 analog")

ax[0,0].legend()
ax[1,0].legend()
ax[0,1].legend()
ax[1,1].legend()
ax[2,1].legend()
ax[2,0].legend()
ax[2,1].set_ylim([10**-24,10**-20])
ax[0,1].set_ylim([10**-24,10**-20])
ax[0,0].set_title("AdLIGO Design")
ax[0,1].set_title("AdLIGO A+")
ax[1,0].set_title("CE1")
ax[1,1].set_title("CE2")
ax[2,0].set_title("AdVIRGO")
ax[2,1].set_title("KAGRA")
plt.savefig("Noise_curves.pdf")
#plt.show()
plt.close()
