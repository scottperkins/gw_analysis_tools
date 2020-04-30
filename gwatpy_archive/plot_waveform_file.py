import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_waveform_image(in_file, out_file):
    data = np.loadtxt(in_file, delimiter=',',unpack=True)
    fig, ax = plt.subplots(nrows=2,ncols=1,sharex=True)
    ax[0].loglog(data[0], np.sqrt(data[1]*data[1]+data[2]*data[2]))
    ax[1].semilogx(data[0], data[1],label="real")
    ax[1].semilogx(data[0], data[2],label="imag")
    ax[1].set_xlabel("Frequency (Hz)")
    ax[0].set_ylabel("Strain")
    ax[1].set_ylabel("Strain")
    ax[1].legend()
    plt.savefig(out_file)
    plt.close()

if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    plot_waveform_image(in_file,out_file);
