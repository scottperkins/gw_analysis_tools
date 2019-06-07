import matplotlib.pyplot as plt
import csv 
import numpy as np
#from phenompy.modified_gr import Modified_IMRPhenomD_Inspiral_Freq_detector_frame as imrmod
from phenompy.modified_gr import dCS_IMRPhenomD_detector_frame as imrmod
#from phenompy.modified_gr import EdGB_IMRPhenomD_detector_frame as imrmod
from phenompy.gr import IMRPhenomD_detector_frame as imr
from phenompy.utilities import s_solm, mpc, calculate_mass1,calculate_mass2
from phenompy.analysis_utilities import log_likelihood_maximized_coal as ll
from time import time

show_plots = True

freqs = []
amp = []
phasec = []
waveform_real = []
waveform_imag = []
chirp_deriv_real = []
chirp_deriv_imag = []
derivA = []
derivp = []
for i in np.arange(8):
    derivA.append([])
    derivp.append([])
with open("data/ppeamplitude_output.csv",'r') as f:
    reader = csv.reader(f, delimiter=',')
    for line in reader:
        freqs.append(float(line[0]))
        amp.append(float(line[1]))
with open("data/ppephase_output.csv",'r') as f:
    reader = csv.reader(f, delimiter=',')
    for line in reader:
        phasec.append(float(line[1]))
with open("data/ppereal_waveform_output.csv",'r') as f:
    reader = csv.reader(f, delimiter=',')
    for line in reader:
        waveform_real.append(float(line[1]))
with open("data/ppeimag_waveform_output.csv",'r') as f:
    reader = csv.reader(f, delimiter=',')
    for line in reader:
        waveform_imag.append(float(line[1]))
#with open("data/ppederiv_amp.csv",'r') as f:
#  reader = csv.reader(f, delimiter=',')
#  for line in reader:
#      derivA[0].append(float(line[1]))
#      derivA[1].append(float(line[2]))
#      derivA[2].append(float(line[3]))
#      derivA[3].append(float(line[4]))
#      derivA[4].append(float(line[5]))
#      derivA[5].append(float(line[6]))
#      derivA[6].append(float(line[7]))
#      derivA[7].append(float(line[8]))
#with open("data/ppederiv_phase.csv",'r') as f:
#    reader = csv.reader(f, delimiter=',')
#    for line in reader:
#        derivp[0].append(float(line[1]))
#        derivp[1].append(float(line[2]))
#        derivp[2].append(float(line[3]))
#        derivp[3].append(float(line[4]))
#        derivp[4].append(float(line[5]))
#        derivp[5].append(float(line[6]))
#        derivp[6].append(float(line[7]))
#        derivp[7].append(float(line[8]))

start = time()
eta = .21
chirpm = 49.78
#model = imrmod(mass1 = calculate_mass1(chirpm,eta)*s_solm, mass2 =calculate_mass2(chirpm,eta)*s_solm, spin1=-.2, spin2=.4, collision_phase=0, collision_time = 0, Luminosity_Distance=410*mpc, phase_mod=np.exp(2), bppe= -1)
model = imrmod(mass1 = calculate_mass1(chirpm,eta)*s_solm, mass2 =calculate_mass2(chirpm,eta)*s_solm, spin1=-.2, spin2=.4, collision_phase=0, collision_time = 0., Luminosity_Distance=410*mpc, phase_mod=(1000/3.e5)**4.)
amppy,phase,h = model.calculate_waveform_vector(freqs)


#i1 = int(len(freqs)/3)
#i2 = int(2*len(freqs)/3)
#freq1 = freqs[i1 ]
#freq2 = freqs[i2]
##negative because phenomd is defined as -phi(f)i
#pycbc1 = phasec[i1]
#pycbc2 = phasec[i2]
#scott1 = phase[i1]
#scott2 = phase[i2]
#tc = (pycbc1 -scott1 - pycbc2 + scott2)/(2*np.pi*(freq1-freq2) )
#phic = -pycbc1 + scott1 + 2*np.pi * freq1*tc


#model = imrmod(mass1 = calculate_mass1(chirpm,eta)*s_solm, mass2 =calculate_mass2(chirpm,eta)*s_solm, spin1=-.2, spin2=.4, collision_phase=phic, collision_time = tc, Luminosity_Distance=410*mpc, phase_mod=np.exp(2))


#model = imr(mass1 = 200*s_solm, mass2 = 50*s_solm, spin1=-.2, spin2=.9, collision_phase=2, collision_time = 8, Luminosity_Distance=800*mpc)
print("py model time: ",time()-start)
start = time()
amppy,phase,h = model.calculate_waveform_vector(freqs)
strain = amppy*np.exp(-1j*phase)
print("py waveform time: ",time()-start)

start = time()
#for i in np.arange(100):
fish,ifish = model.calculate_fisher_matrix_vector("Hanford_O1",lower_freq=freqs[0], upper_freq=freqs[-1], stepsize=freqs[1]-freqs[0])
print("Fish time: ",time()-start)
#print(fish)

#py_deriv = []
#py_derivp = []
#start = time()
#model.calculate_derivatives()
#for i in range(1,9):
#    py_deriv.append( (model.calculate_waveform_derivative_vector(freqs,i)).real)
#    py_derivp.append(- (model.calculate_waveform_derivative_vector(freqs,i)).imag/amppy)
#print("deriv time: ",time()-start)
#py_deriv[0] = py_deriv[0]*model.A0;
#py_deriv[3] = py_deriv[3]*model.chirpm;
#py_deriv[4] = py_deriv[4]*model.symmratio;
#
#py_derivp[0] = py_derivp[0]*model.A0;
#py_derivp[3] = py_derivp[3]*model.chirpm;
#py_derivp[4] = py_derivp[4]*model.symmratio;
#
#for i in np.arange(0,8):
#    plt.plot(freqs, 2*(py_deriv[i] - derivA[i])/(derivA[i]+py_deriv[i]), label=i)
#    plt.plot(freqs, 2*(py_derivp[i] - derivp[i])/(derivp[i]+py_derivp[i]), label=i)
#plt.legend()
#plt.xscale('log')
#if show_plots:
#    plt.show()
#plt.close()





plt.plot(freqs,waveform_imag)
plt.plot(freqs,strain.imag)
if show_plots:
    plt.show()
plt.close()

plt.plot(freqs,(strain.real-waveform_real)/amppy)
plt.savefig("real_waveform_diff.pdf")
if show_plots:
    plt.show()
plt.close()

plt.plot(freqs,(strain.imag-waveform_imag)/amppy)
if show_plots:
    plt.show()
plt.close()

plt.loglog(freqs,amp)
plt.loglog(freqs, amppy)
if show_plots:
    plt.show()
plt.close()

plt.plot(freqs,phasec,label='c')
plt.plot(freqs, phase,label='py')
plt.legend()
if show_plots:
    plt.show()
plt.close()

plt.plot(freqs,(amp-amppy)/amppy)
if show_plots:
    plt.show()
plt.close()
plt.plot(freqs,(phasec-phase)/phase)
if show_plots:
    plt.show()
plt.close()
