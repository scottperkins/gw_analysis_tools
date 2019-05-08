import numpy as np
import matplotlib.pyplot as plt
from chainconsumer import ChainConsumer

data = np.loadtxt("data/mcmc_output.csv",delimiter=',')
autocorr = np.loadtxt("data/auto_corr_mcmc.csv",delimiter=',')
lengths = autocorr[0]
autocorr = autocorr[1:]
#for i in np.arange(len(data)):
#    data[i][0] = np.log(data[i][0])
#    data[i][1] = np.log(data[i][1])
for i in autocorr:
    plt.plot(lengths,i)
plt.show()
plt.close()
burnin =0
data = data[burnin:]
#burnin =0;
#plt.plot(data[0][burnin:])
#plt.plot(data[1][burnin:])
#plt.show()
#plt.close()
##plt.hist(data,bins=100,density=True)
##x = np.linspace(-3,3)
#plt.hist2d(data[0][burnin:],data[1][burnin:],bins=200)
##plt.plot(x,np.exp(-x**2/4.)/np.sqrt(4*np.pi))
##plt.plot(x,np.exp(-x**2/10.)/np.sqrt(10*np.pi))
#plt.show()
#plt.close()
#plt.hist(data[0][burnin:],bins=100,density=True)
#plt.show()
#plt.close()
#plt.hist(data[1][burnin:],bins=100,density=True)
#plt.show()
#plt.close()
c = ChainConsumer()
#c.add_chain(data, parameters=["${\cal M}$", "$M$", "$\chi_{eff}$"])
#c.add_chain(data, parameters=["${\cal M}$", "$M$", "$\chi_{eff}$","$m_g$"])
c.add_chain(data, parameters=["DL","${\cal M}$", "$\eta$", "$\chi_{1}$","$\chi_{2}$"])
#c.add_chain(data, parameters=["${\cal M}$", "$M$", "$\chi_{eff}$","$beta_{-1}$"])
c.configure(summary=True, tick_font_size=8, label_font_size=10, max_ticks=4, sigmas=[0,1,2,3])

# If we wanted to save to file, we would instead have written
fig = c.plotter.plot(filename="mcmc_testing", figsize="column")

#fig = c.plotter.plot(filename="intrinsic", figsize="column", truth=[1.0, 1.0, 3.1416], extents=[[0.5,1.5],[0.996,1.002],[2.5,4.8]])

#fig = c.plotter.plot(filename="intrinsic", figsize="column", truth=[0.7071, 1.0, 3.1416], extents=[[0,1.2],[0.5,2.5],[0,6.283]])


# If we wanted to display the plot interactively...
#fig = c.plotter.plot(display=True, figsize="column")

fig.set_size_inches(4.5 + fig.get_size_inches())
