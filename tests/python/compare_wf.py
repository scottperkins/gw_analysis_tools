import numpy as np
import matplotlib.pyplot as plt
from time import time
import matplotlib.pyplot as plt
import astropy.cosmology as cosmology
from astropy.coordinates import Distance
from astropy import units as u
import astropy.constants as consts
from pycbc.types import frequencyseries
from pycbc.waveform import get_fd_waveform_sequence, get_fd_waveform, get_td_waveform
from pycbc.waveform import fd_approximants, utils
from pycbc.types import Array

gwatp = np.loadtxt("data/hpgwat.csv",delimiter=',',unpack=True)
lalp = np.loadtxt("data/hpLAL.csv",delimiter=',',unpack=True)
freq = np.loadtxt("data/freqs.csv",delimiter=',',unpack=True)
plt.plot(freq,gwatp[0],label="gwat p")
plt.plot(freq,lalp[0],label="LAL p",linestyle='-.')
plt.legend()
plt.savefig("plots/wfd.pdf")
#plt.show()
plt.close()

ampgwatp = gwatp[0]*gwatp[0] + gwatp[1]*gwatp[1]
ampLALp = lalp[0]*lalp[0] + lalp[1]*lalp[1]
plt.semilogy(freq,ampgwatp,label="gwat p")
plt.semilogy(freq,ampLALp,label="LAL p")
plt.legend()
plt.savefig("plots/ampd.pdf")
#plt.show()
plt.close()

plt.plot(freq,(ampgwatp-ampLALp)/ampLALp,label="plus")
#plt.plot(freq,(ampgwatp-ampLALp),label="plus")
#plt.xlim([0,400])
plt.legend()
plt.savefig("plots/ampddiff.pdf")
#plt.show()
plt.close()

delta_f = freq[1]-freq[0]
gwatseriesplus = frequencyseries.FrequencySeries(gwatp[0]+ 1j*gwatp[1],delta_f=delta_f)
phasegwatp = utils.phase_from_frequencyseries(gwatseriesplus)
LALseriesplus = frequencyseries.FrequencySeries(lalp[0]+ 1j*lalp[1],delta_f=delta_f)
phaseLALp = utils.phase_from_frequencyseries(LALseriesplus)
gwatpp = np.loadtxt("data/ppgwat.csv",delimiter=',',unpack=True)
#phasegwatp = np.arctan(gwatp[1]/gwatp[0] )
#phaseLALp = np.arctan(lalp[1]/lalp[0] )
#phasegwatc = np.arctan(gwatc[1]/gwatc[0] )
#phaseLALc = np.arctan(lalc[1]/lalc[0] )
plt.plot(freq,phasegwatp,label="gwat p")
plt.plot(freq,phaseLALp,label="LAL p")
#plt.plot(freq,gwatpp+(phasegwatp[0]-gwatpp[0]),label="gwat p -- native")
plt.legend()
plt.savefig("plots/phased.pdf")
#plt.show()
plt.close()

#plt.plot(freq,(phasegwatp-phaseLALp)/phaseLALp,label="plus")
#plt.plot(freq,(phasegwatp-phaseLALp),label="plus")
#plt.plot(freq,(phasegwatc-phaseLALc),label="cross")
plt.plot(freq,(phasegwatp-phaseLALp),label="plus")
plt.legend()
plt.savefig("plots/phaseddiff.pdf")
#plt.show()
plt.close()
