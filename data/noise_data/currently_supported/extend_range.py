import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

range_TERR = np.asarray([2.,10000],dtype=np.float64)

def fitting_func(x, logrootS,  b,c):
    return logrootS +  b*x +c*x*x

def fitting_func2(x, logrootS,  b,c,d):
    return logrootS +  b*x +c*x*x+d*x*x*x

#files = ["aligo_O4high.csv","CE1_strain_smoothed.csv","aligo_O4high_smoothed.csv","AplusDesign.csv","AplusDesign_smoothed.csv", "AdLIGODesign.csv","AdLIGOMidHigh.csv","avirgo_O4high_NEW.csv","avirgo_O4high_NEW_smoothed.csv","avirgo_O5high_NEW.csv","avirgo_O5high_NEW_smoothed.csv","avirgo_O5low_NEW.csv","avirgo_O5low_NEW_smoothed.csv","CE1_strain.csv","CE1_strain_smoothed.csv","CE2_strain.csv","CE2_strain_smoothed.csv","Hanford_O2_Strain.csv","kagra_128Mpc.csv","kagra_25Mpc.csv","kagra_80Mpc.csv"]
files = ["aligo_O4high.csv","CE1_strain_smoothed.csv","aligo_O4high_smoothed.csv","AplusDesign.csv","AplusDesign_smoothed.csv", "AdLIGODesign.csv","AdLIGOMidHigh.csv","avirgo_O4high_NEW.csv","avirgo_O4high_NEW_smoothed.csv","avirgo_O5high_NEW.csv","avirgo_O5high_NEW_smoothed.csv","avirgo_O5low_NEW.csv","avirgo_O5low_NEW_smoothed.csv","CE1_strain.csv","CE1_strain_smoothed.csv","CE2_strain.csv","CE2_strain_smoothed.csv","kagra_128Mpc.csv","kagra_25Mpc.csv","kagra_80Mpc.csv"]
#files = ["AdLIGODesign.csv"]
#files = ["aligo_O4high.csv"]
for i in np.arange(len(files)):
    file_i = "../currently_supported_raw/" + files[i]
    dataraw = np.loadtxt(file_i, delimiter=',',unpack=True)
    upper_fit_bound = 50
    #Better results
    #if(files[i] == "AplusDesign.csv" or files[i] == "AdLIGODesign.csv" or files[i] == "AdLIGOMidHigh.csv" or files[i] == "AplusDesign_smoothed.csv" or files[i] == "AdLIGODesign_smoothed.csv" or files[i] == "AdLIGOMidHigh_smoothed.csv"):
    if(files[i] == "AdLIGODesign.csv" or files[i] == "AdLIGOMidHigh.csv" or  files[i] == "AdLIGODesign_smoothed.csv" or files[i] == "AdLIGOMidHigh_smoothed.csv"):
        upper_fit_bound = 20
        print(files[i])
    if(files[i] == "AplusDesign.csv" or files[i] == "AplusDesign_smoothed.csv"):
        upper_fit_bound = 10
        print(files[i])



    deltaf = dataraw[0][1]-dataraw[0][0]
    data_out = []
    data_out.append(dataraw[0])
    f = dataraw[0][0] - deltaf
    ct = 0
    while data_out[0][0] > range_TERR[0]:
        data_out[0] =np.insert( data_out[0],0,f)
        f-= deltaf
    deltaf = dataraw[0][-1]-dataraw[0][-2]
    f = dataraw[0][-1] + deltaf
    while data_out[0][-1] < range_TERR[1]:
        data_out[0] =np.insert( data_out[0],len(data_out[0]),f)
        f+= deltaf

    
    #popt, pcov = curve_fit(f=fitting_func, xdata=np.log(dataraw[0]),ydata=np.log(dataraw[1]))
    popt, pcov = curve_fit(f=fitting_func, xdata=np.log(dataraw[0][np.where(dataraw[0]<upper_fit_bound)]),ydata=np.log(dataraw[1][np.where(dataraw[0]<upper_fit_bound)]))
    #newrange = np.logspace(np.log10(range_TERR[0]),np.log10(range_TERR[1]),10000)
    newrange = np.logspace(np.log10(range_TERR[0]),np.log10(upper_fit_bound),10000)
    y = np.exp(fitting_func(np.log(newrange), *popt))
    #plt.loglog(newrange,y)


    data_out.append(np.array([]))
    ct = 0
    for x in data_out[0][:int(len(data_out[0])/2)]:
        if x in dataraw[0]:
            data_out[1] = np.insert(data_out[1],len(data_out[1]),dataraw[1][ct])
            ct+=1
        else:
            y =np.exp(fitting_func(np.log(x),*popt))
            if(y > dataraw[1][0]):
                data_out[1] = np.insert(data_out[1],len(data_out[1]),y)
            else:
                data_out[1] = np.insert(data_out[1],len(data_out[1]),dataraw[1][0])

    
    popt, pcov = curve_fit(f=fitting_func, xdata=np.log(dataraw[0][np.where(dataraw[0]>3000)]),ydata=np.log(dataraw[1][np.where(dataraw[0]>3000)]))
    newrange = np.logspace(np.log10(3000),np.log10(range_TERR[1]),10000)
    y = np.exp(fitting_func(np.log(newrange), *popt))
    #plt.loglog(newrange,y)

    
    for x in data_out[0][int(len(data_out[0])/2):]:
        if x in dataraw[0]:
            data_out[1] = np.insert(data_out[1],len(data_out[1]),dataraw[1][ct])
            ct+=1
        else:
            y =np.exp(fitting_func(np.log(x),*popt))
            if(y > dataraw[1][-1]):
                data_out[1] = np.insert(data_out[1],len(data_out[1]),y)
            else:
                data_out[1] = np.insert(data_out[1],len(data_out[1]),dataraw[1][-1])
    plt.loglog(dataraw[0],dataraw[1])
    plt.loglog(data_out[0],data_out[1],linestyle='--')
    plt.title(file_i)
    plt.savefig("plots/extended_comparison_{}.pdf".format(files[i]))
    #plt.show()
    plt.close()
    np.savetxt(files[i],np.transpose(data_out),delimiter=',')

