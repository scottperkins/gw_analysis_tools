import numpy as np
import matplotlib.pyplot as plt

def AC(data):
    h = 1
    N = len(data)
    acs = np.ones(int(N*.75))
    mean = np.mean(data)
    for l in np.arange(int(N*.75)):
        acs[l] = 1./(N-h) * np.sum( (data[h:]-mean)*(data[:N-h]-mean) )
        h +=1
    return acs

