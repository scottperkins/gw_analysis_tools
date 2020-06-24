import astropy.cosmology as cosmos
from scipy.integrate import quad
import numpy as np
import matplotlib.pyplot as plt
c = 299792458.
MPC_SEC = 3.085677581491367278913937957796471611e22/c

cosmo = cosmos.Planck15

def helper(alpha, z):
    return (1+z)**(1-alpha)*quad(lambda x: (1+x)**(alpha-2)/cosmo.H(x).to("Hz").value , 0 , z)[0]

alphas = ["0.000000","0.500000","1.000000","1.500000","2.000000","2.500000","3.000000","3.500000","4.000000"]
for alpha in alphas:
    data = np.loadtxt("data/MD_DL_from_Z_{}.csv".format(alpha),delimiter=',',unpack=True)
    fracdiff = []
    for x in data:
        DL= helper(float(alpha),x[0])/MPC_SEC
        fracdiff.append( (x[1] -DL  )/ (DL))
    plt.hist(fracdiff, bins = np.logspace(-8,0,100),label=str(alpha))
plt.legend()
plt.xscale('log')
plt.savefig("plots/alpha_MD_comp.pdf")
plt.close()
