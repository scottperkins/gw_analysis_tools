import numpy as np
import matplotlib.pyplot as plt
import sys

def plot_waveform_image(in_file, out_file):
    data = np.loadtxt(in_file, delimiter=',',unpack=True)
    plt.loglog(data[0], np.sqrt(data[1]*data[1]+data[2]*data[2]))
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Strain")
    plt.savefig(out_file)
    plt.close()

if __name__ == "__main__":
    in_file = sys.argv[1]
    out_file = sys.argv[2]
    plot_waveform_image(in_file,out_file);
