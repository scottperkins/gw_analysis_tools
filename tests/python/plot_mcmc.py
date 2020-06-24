import matplotlib as mpl
#mpl.use("pdf")
import corner
import matplotlib.pyplot as plt
import numpy as np
#import gwatpy.util as gpu
#import gwatpy.gwatpy_plot as gp; gp.set()
#from phenompy.utilities import calculate_mass1, calculate_mass2

data = np.loadtxt("data/injection_output.csv",delimiter=',')
#data = np.loadtxt("data/experiment_output.csv",delimiter=',')
#data = np.loadtxt("data/test_output.csv",delimiter=',')
#injections = np.loadtxt("data/injections.csv",delimiter=',',unpack=True)
injections = np.loadtxt("data/injections.csv",delimiter=',',unpack=True)
#dim = 11
dim = 15

#data = data[100:]
for x in data:
    #x[1] = np.arcsin(x[1])
    #x[3] = np.arccos(x[3])
    #x[1] = np.arcsin(x[1])*180./np.pi
    x[6] = np.exp(x[6])
    x[7] = np.exp(x[7])
    #x[0] = np.exp(x[0])
#injections[1] = np.arcsin(injections[1])
#injections[3] = np.arccos(injections[3])
injections[6] = np.exp(injections[6])
injections[7] = np.exp(injections[7])
data_thinned = []
for x in np.arange(len(data)):
    if x%1 ==0:
        data_thinned.append(data[x]) 
#for i in np.arange(8):
#    parameter = [x[i] for x in data]
#    plt.plot(parameter)
#    plt.show()
#    plt.close()
ndim, nsamples = 11, len(data) 
labels = [r"$\alpha$",r"$\sin(\delta)$",r"$\psi$",r"$\cos(\iota)$","$\phi_{ref}$","$t_c$",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$a_{1}$",r"$a_2$",r"$\cos \theta_1$",r"$\cos \theta_2$",r"$\phi_p$",r"$\sqrt{\alpha}$"]
#labels = [r"$\alpha$",r"$\sin(\delta)$",r"$\psi$",r"$\iota$","$\phi_{ref}$","$t_c$",r"$D_L$",r"$\mathcal{M}$",r"$\eta$",r"$\chi_{1}$",r"$\chi_2$"]
data_plot=[]
for x in data_thinned:
    #chi1 = x[2]*(x[4])
    #chi2 = x[3]*(x[5])
    #chi1 = x[2]
    #chi2 = x[3]
    #m1 = gpu.calculate_mass1_py(x[0],x[1])
    #m2 = gpu.calculate_mass2_py(x[0],x[1])
    ###chi1 = x[9]*(x[11])
    ###chi2 = x[10]*(x[12])
    ###chi1 = x[9]
    ###chi2 = x[10]
    #m1 = gpu.calculate_mass1_py(x[7],x[8])
    #m2 = gpu.calculate_mass2_py(x[7],x[8])
    #x=np.append(x,(chi1*m1 + chi2*m2 ) /(m1+m2))
    data_plot.append(x)
data_plot = np.asarray(data_plot)
figure = corner.corner(data_plot, labels=labels,quantiles=[.1,.5,.9], show_titles=True)
axes = np.array(figure.axes).reshape(dim,dim)
for i in np.arange(dim):
    ax = axes[i,i]
    ax.axvline(injections[i])

for yi in np.arange(dim):
    for xi in np.arange(yi):
        ax = axes[yi,xi]
        ax.axvline(injections[xi])
        ax.axhline(injections[yi])
        ax.plot(injections[xi],injections[yi])

plt.savefig("plots/mcmc_injection.pdf")
#plt.savefig("plots/mcmc_experiment.pdf")
plt.close()
