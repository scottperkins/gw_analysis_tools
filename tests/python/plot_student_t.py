import matplotlib as mpl
#mpl.use("pdf")
import corner
import matplotlib.pyplot as plt
import numpy as np
#import gwatpy.util as gpu
#import gwatpy.gwatpy_plot as gp; gp.set()
#from phenompy.utilities import calculate_mass1, calculate_mass2

data = np.loadtxt("data/mcmc_output.csv",delimiter=',')

labels = ["X","Y"]
data_plot=[]
for x in np.arange(len(data)):
    if x% 1 == 0:
        data_plot.append(data[x])
data_plot = np.asarray(data_plot)
figure = corner.corner(data_plot, labels=labels,quantiles=[.16,.5,.84], show_titles=True)

plt.savefig("plots/mcmc_student_t.pdf")
plt.close()
