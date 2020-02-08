import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import gwatpy.gwatpy_plot as gp

x = np.linspace(0,1)
y = np.sin(x)
y2 = np.cos(x)
y3 = np.arctan(x)
y4 = np.exp(x)
y5 = x**2
y6 = x**3
y7 = x**4
y8 = x**5
y9 = x**.5
plt.plot(x,y,label="sin")
plt.plot(x,y2,label="cos")
plt.plot(x,y3,label="arctan")
plt.plot(x,y4,label="exp")
plt.plot(x,y5,label="square")
plt.plot(x,y6,label="cube")
plt.plot(x,y7,label="quad")
plt.plot(x,y8,label="fifth")
plt.plot(x,y9,label="root")
plt.legend()
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Title")
plt.savefig("sample.pdf")
plt.close()

gp.show_minor_grid()
x = np.random.normal(size=10000)
y = np.random.normal(loc=1,size=10000)
z = np.random.normal(loc=-1,size=10000)
x = np.abs(x)
y = np.abs(y)
z = np.abs(z)
plt.hist(x,alpha=.8,bins=np.logspace(np.log10(0.001),np.log10(5),100),label='Mean 0')
plt.hist(y,alpha=.8,bins=np.logspace(np.log10(0.001),np.log10(5),100),label='Mean 1')
plt.hist(z,alpha=.8,bins=np.logspace(np.log10(0.001),np.log10(5),100),label='Mean -1')
plt.xscale('log')
plt.legend()
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Title")
plt.savefig("sample_hist.pdf")
plt.close()

gp.show_minor_grid()
x = np.linspace(0,2*np.pi)
y = np.sin(x)
plt.loglog(x,y)
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Title")
plt.savefig("sample_log.pdf")
plt.close()

gp.show_major_grid()


########################################
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10,5))
ax.plot(x,y)
plt.savefig("sample_fig.pdf")
plt.close()
########################################

