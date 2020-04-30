import numpy as np
import gwatpy.util as gpu
import gwatpy.detector_util as gpdu
import matplotlib.pyplot as plt

print(gpu.calculate_chirpmass_py(20,10))

freqs = np.linspace(10,1000,1000)
detector = "AdLIGODesign"
length = 1000
integration_time = 48
psd_root = gpdu.populate_noise_py(freqs, detector, length , integration_time)
plt.loglog(freqs, psd_root)
plt.show()
plt.close()
