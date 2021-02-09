import numpy as np
import matplotlib.pyplot as plt
import os

files = [f for f in os.listdir("data/") if "swap" in f]
print(files)

dat = np.loadtxt("data/"+files[0],delimiter=',')
for f in files[1:]:
    print(f)
    datt = np.loadtxt("data/"+f,delimiter=',')
    dat+= datt

im = plt.imshow(dat,cmap='hot',interpolation='nearest',origin=True)
ax = plt.gca()
cbar = ax.figure.colorbar(im,ax=ax)

plt.savefig("plots/heatmap_temp.pdf")
plt.close()
