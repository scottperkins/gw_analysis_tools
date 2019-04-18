import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation, patches
from gw_analysis_tools_py import mcmc_routines_ext as mcmc
from gw_analysis_tools_py import waveform_generator_ext as wg
from phenompy.utilities import mpc, s_solm, calculate_mass1,calculate_mass2
import seaborn as sns

sns.set(style='dark')
sns.set_palette('colorblind')
#sns.set_context('talk')

method = b"ppE_IMRPhenomD_Inspiral"
screening = False
fmin = 28
fmax = 29
#fmax = 200
steps = 1600
deltaf = (fmax-fmin)/steps
#T = 1/deltaf
T = 0
tmin = 0
tmax = T
#tc =-T/4
tc = 0

betascale = -1

framenum = 100
framenum_min = int(framenum*.30)
framenum_unscreened = int(framenum*.4)
framenum_max = int(framenum*.7)

figdim = [8,6]
# First set up the figure, the axis, and the plot element we want to animate
#fig = plt.figure()
fig, (ax1,ax2) = plt.subplots(nrows=2, ncols=1, figsize=[8,7])
#fig, (ax2) = plt.subplots(nrows=1, ncols=1,figsize=figdim)
#ax1 = plt.axes(xlim=(0, 2), ylim=(-2, 2))
ax1.set_xlim([fmin,fmax])
#ax1.set_ylim([-0,-100000])
ax2.set_xlim([fmin,fmax])
ax2.set_ylim([-2e-22,2e-22])
line2, = ax2.plot([],[],lw=2)

eta = .21
chirpmass = 49.78
spin1 = 0
spin2 =0
DL=410
dl_step = DL/(framenum_unscreened)
dl_steptot = DL/framenum
theta =122.2 * np.pi/180
phi = 150.1*np.pi/180
iota = 0*np.pi/2
beta = np.ascontiguousarray([-1], dtype=np.float64)
b = np.ascontiguousarray([-3],dtype=np.intc)
param = wg.gen_params_py(calculate_mass1(chirpmass,eta),calculate_mass2(chirpmass,eta),DL,
        np.ascontiguousarray([0,0,spin1]),np.ascontiguousarray([0,0,spin2]),T/2,0,b,beta,1,0,0,iota,100,0,False)
#plt.close()
frequencies = np.ascontiguousarray(np.linspace(fmin,fmax,steps),dtype=np.float64)
waveform = wg.fourier_waveform_py(frequencies, b"IMRPhenomD",param)
#plt.plot(np.linspace(tmin,tmax, steps), np.fft.fft(waveform).real)
#plt.show()
#plt.close() 
y= np.asarray(wg.fourier_phase_py(frequencies, b"ppE_IMRPhenomD_Inspiral",param))
x2 = frequencies 
line1, = ax1.plot(x2,y,lw=2)

plt.tight_layout()
def init():
    line2.set_data([], [])
    line1.set_data([], [])

    return line2,

# animation function.  This is called sequentially
def animate(i):

    i = i+1
    beta = np.ascontiguousarray([betascale*dl_step*(i)], dtype=np.float64)

    param = wg.gen_params_py(calculate_mass1(chirpmass,eta),calculate_mass2(chirpmass,eta),DL,
            np.ascontiguousarray([0,0,spin1]),np.ascontiguousarray([0,0,spin2]),
            0,tc,b,beta,1,0,0,iota,100,0,False)
    y2= wg.fourier_waveform_py(frequencies, method,param).real
    y= np.asarray(wg.fourier_phase_py(frequencies, method,param))
    x2 = frequencies 
    for t in np.arange(len(x2)):
        if x2[t]>28.7 and x2[t]<28.9:
            print(y[t],x2[t])
    print()
    line2.set_data(x2,y2)
    line1.set_data(x2,y)
    ax1.set_ylim([y[0]*1.1,y[-1]])
    #line2.set_data(frequencies,y2.real)
    return line2,

# call the animator.  blit=True means only re-draw the parts that have changed.

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=framenum, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('test.mp4', fps=5, extra_args=['-vcodec', 'libx264'])

#plt.show()

