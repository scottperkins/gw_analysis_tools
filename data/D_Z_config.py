##########################################################
#This code writes a header file to store numerical fit 
#data to calculate redshift from luminosity distance 
#quicker and more efficiently. To do so, I use the 
#scipy.optimize.curve_fit function to fit coefficients 
#to the function defined below (power series in 1/2 power 
#increments) which is then written out to a header file
#which is in turn included in the util.cpp file.

#This file shouldn't need to be run once the header file is 
#created once, but its available to update cosmology information
#or to do cosmological studies
##########################################################

import numpy as np
import multiprocessing as mp
import csv
import astropy.cosmology as cosmos
import astropy.units as u
from astropy.coordinates import Distance
from scipy.integrate import quad
import matplotlib.pyplot as plt
#from IMRPhenomD import mpc
import os
from time import time
import scipy 
from scipy import optimize
c = 299792458.
MPC_SEC = 3.085677581491367278913937957796471611e22/c



deg = 11
#Function being fit to the data -- power series in steps of 1/2 power
def func(x, a, b ,c, d, e, f, g, h, i, j ,k):
    return  (a*np.sqrt(x)  + b*x + c*x**(3./2.) + d*x*x + e*x**(5./2.) + f*x*x*x  +
            g*x**(7./2.) + h*x*x*x*x + i*x**(9./2.) + j*x**(10./2.) + k*x**(11./2.) )
            
if __name__=="__main__":

    cosmologies = [cosmos.Planck15, cosmos.Planck13, cosmos.WMAP9,cosmos.WMAP7, cosmos.WMAP5]
    cosmology_names = ["PLANCK15", "PLANCK13", "WMAP9","WMAP7", "WMAP5"]

    data_filetree = os.path.dirname(os.path.realpath(__file__))
    include_filetree = data_filetree[:-4] + "include"
    config_header = open(include_filetree+'/D_Z_Config.h','w')
    config_header.write("#ifndef D_Z_CONFIG_H\n")
    config_header.write("#define D_Z_CONFIG_H\n")
    config_header.write("\n")

    for k in np.arange(len(cosmologies)):

        cosmo = cosmologies[k]
        cosmo_name = cosmology_names[k]
        def H(z):
            return cosmo.H(z).to('Hz').value
        def helper(y):
            return (1+y)*quad(lambda x: 1/H(x),0,y)[0]
        #def helperold(y):
        #    return Distance(y,unit=u.Mpc).compute_z(cosmology = cosmo)
        #def helper2(zfinal):
        #    return quad(lambda x: 1/((1+x)**2*H(x)),0,zfinal )[0]
        pool = mp.Pool(processes=mp.cpu_count())

        # dl = np.linspace(1e-3,500,1000)
        # y1 = helper(dl)
        # y2 = list(map(helperold,dl))
        # # plt.plot(dl,helper(dl),label='interpolated')
        # # plt.plot(dl,list(map(helperold,dl)),label='astropy')
        # plt.plot(dl,(y1-y2)/y2,label='astropy')
        # plt.legend()
        # plt.show()
        # plt.close()
        # H0=cosmos.Planck15.H0.to('Hz').value


        #############################################################
        #This code is used to tabulate data for the mapping of redshift to luminosity Distance
        #to be interpolated later for speed
        #File has the form: LumD , Z
        #LumD ranges from 1 to 50000 MPC
        #############################################################

        pts_seg = 1e3
        z1 = np.logspace(np.log10(1e-6),np.log10(.1),pts_seg)
        z2 = np.logspace(np.log10(.1),np.log10(1),pts_seg)
        z3 = np.logspace(np.log10(1),np.log10(20),pts_seg)
        zs = [z1,z2,z3]
        lds =[] 
        for i in np.arange(len(zs)): 
            lds.append(np.asarray(pool.map(helper, zs[i]))/MPC_SEC)

        coeffs = []
        for i in np.arange(len(zs)):
            popt, pcov = scipy.optimize.curve_fit(func,lds[i], zs[i])
            coeffs.append(popt)

        for i in np.arange(len(coeffs)):
            config_header.write("const double {}_COEFF_VEC_{}[{}] =  {{".format(cosmo_name, i, deg))
            for j in np.arange(len(coeffs[0])):
                config_header.write("{}, ".format(coeffs[i][j]))
            config_header.write("}}\n".format(cosmo_name))
        
        for i in np.arange(len(coeffs)):
            config_header.write("#define {}_DL_MIN_{} {}\n".format(cosmo_name, i, lds[i][0]))
            config_header.write("#define {}_DL_MAX_{} {}\n".format(cosmo_name, i, lds[i][-1]))

        for i in np.arange(len(coeffs)):
            interpolated = func(lds[i], *coeffs[i]) 
            print(interpolated)
            plt.plot(lds[i], (interpolated-zs[i])/zs[i])
        plt.title("{} Residuals".format(cosmo_name))
        plt.xlabel("Luminosity Distance (MPC)")
        plt.ylabel(r"Fractional Difference $(Z_{true} - Z_{interp})/Z_{true}$")
        plt.tight_layout()
        plt.savefig(data_filetree+"/{}_residuals.png".format(cosmo_name))
        plt.close()

        config_header.write("#define {}_SEG_NUM {}\n".format(cosmo_name,len(coeffs)))
        config_header.write("\n")
    config_header.write("#endif")
    config_header.close()








    #############################################################
    #This code is used to tabulate data for the mapping of redshift to cosmological distance D/(1+Z)
    # defined in Will '97 to be interpolated later for speed
    # D/(1+Z) = integral^Z_0 dz/((1+Z)**2*np.sqrt(.3(1+Z)**3 + .7))
    #using the Lambda CDM universe with Omega_M = 0.3 and Omega_Lambda =0.7
    #File has the form: Z , D
    #Z ranges from 0 to 1.5
    #############################################################
    #z = np.linspace(0,10,1e5)
    #d = np.asarray(pool.map(helper2, z))
    #with open(data_filetree+'/tabulated_Z_D.csv','w') as file:
    #    writer = csv.writer(file, delimiter=',')
    #    row = [[z[t],d[t]/mpc] for t in np.arange(len(d))]
    #    for i in row:
    #        writer.writerow(list(i))

    #TESTING the accuracy of the second set of data
    # Zcheck = np.linspace(0,.5,100)
    # start = time()
    # d1 = np.asarray(list(map(helper2,Zcheck)))
    # print(time()-start)
    # start = time()
    # dfunc = interp1d(z,d/mpc)
    # d2 = np.asarray(dfunc(Zcheck))
    # print(time()-start)
    # print((d1/mpc-d2))
    #
    # plt.plot(Zcheck,dfunc(Zcheck))
    # plt.scatter(z,np.divide(d,mpc))
    # plt.show()
    # plt.close()
