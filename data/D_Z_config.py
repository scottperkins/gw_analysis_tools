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



#Function being fit to the data -- power series in steps of 1/2 power
def func(x, l,a, b ,c, d, e, f, g, h, i, j ,k):
    return  l+(a*np.sqrt(x)  + b*x + c*x**(3./2.) + d*x*x + e*x**(5./2.) + f*x*x*x  +
            g*x**(7./2.) + h*x*x*x*x + i*x**(9./2.) + j*x**(10./2.) + k*x**(11./2.) )
            
if __name__=="__main__":

    data_filetree = os.path.dirname(os.path.realpath(__file__))
    include_filetree = data_filetree[:-4] + "include"
    #config_header = open(include_filetree+'/D_Z_Config.h','w')
    with open(include_filetree+'/D_Z_Config.h','w') as config_header:
        cosmologies = [cosmos.Planck15, cosmos.Planck13, cosmos.WMAP9,cosmos.WMAP7, cosmos.WMAP5]
        cosmology_names = ["PLANCK15", "PLANCK13", "WMAP9","WMAP7", "WMAP5"]
        deg = np.ones((len(cosmologies),))*12 #Number of parameters in the fitting function
        num_segments = np.ones((len(cosmologies),))*3 #This is how many chunks the interpolation range is broken up into 
        pts_seg = []
        for i in np.arange(len(cosmologies)):
            pts_seg.append([])
            for j in np.arange(num_segments[i]):
                pts_seg[-1] .append(1e3)#points in each segment
        boundariesZ_base = []
        for i in np.arange(len(cosmologies)):
            boundariesZ_base .append(np.logspace(np.log10(1e-6),np.log10(20), num_segments[i]+1))#points in each segment
        #boundariesZ_base = np.logspace(np.log10(1e-6),np.log10(20), num_segments+1) #Written based on chunks of Z, then boundaries of DL will be dependent on cosmology

        config_header.write("#ifndef D_Z_CONFIG_H\n")
        config_header.write("#define D_Z_CONFIG_H\n")
        config_header.write("\n")
        config_header.write("// I am aware this isn't very human-readable -- python has issues with line sizes, so the lines had to be truncated in some fashion\n")

        config_header.write("const char * cosmos[{}] =  {{".format(len(cosmology_names)))
        config_header.write("\"{}\"".format(cosmology_names[0]))
        for i in np.arange(len(cosmology_names)-1):
            config_header.write(", \"{}\"".format(cosmology_names[i+1]))
        config_header.write("};\n")
        
        boundariesZ =[]
        boundariesD =[]
        coeffsZ =[]
        coeffsD =[]
        for i in np.arange(len(cosmologies)):
            #boundariesZ.append([])
            #boundariesD.append([])
            #coeffsZ.append([])
            #coeffsD.append([])
            boundariesZ.append(boundariesZ_base[i])
            boundariesD.append(np.zeros((int(num_segments[i]),)))
            coeffsZ.append(np.zeros((int(num_segments[i]),int(deg[i]))))
            coeffsD.append(np.zeros((int(num_segments[i]),int(deg[i]))))

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

            #pts_seg = 1e3
            #z1 = np.logspace(np.log10(1e-6),np.log10(.1),pts_seg)
            #z2 = np.logspace(np.log10(.1),np.log10(1),pts_seg)
            #z3 = np.logspace(np.log10(1),np.log10(20),pts_seg)
            #zs = [z1,z2,z3]
            zs = []
            for x in np.arange(int(num_segments[k])): 
                zs.append(np.logspace(np.log10(boundariesZ[k][x]),np.log10(boundariesZ[k][x+1]),int(pts_seg[k][x])))
            lds =[] 
            for i in np.arange(len(zs)): 
                lds.append(np.asarray(pool.map(helper, zs[i]))/MPC_SEC)

            ##########################################################
            ##### redshift to luminosity distance ###############
            ##########################################################
            coeffsDZ = []
            for i in np.arange(len(zs)):
                popt, pcov = scipy.optimize.curve_fit(func,lds[i], zs[i])
                coeffsDZ.append(popt)

            #for i in np.arange(len(coeffsDZ)):
            #    #config_header.write("const double {}_COEFF_VEC_DZ_{}[{}] =  {{".format(cosmo_name, i, deg))
            #    config_header.write("const double {}_COEFF_VEC_DZ_{}[{}] =  {{".format(k, i, deg))
            #    config_header.write("{} ".format(coeffsDZ[i][0]))
            #    for j in np.arange(len(coeffsDZ[0])-1):
            #        config_header.write(", {}".format(coeffsDZ[i][j+1]))
            #    config_header.write("};\n")
            #
            #for i in np.arange(len(coeffsDZ)):
            #    #config_header.write("#define {}_DL_MIN_{} {}\n".format(cosmo_name, i, lds[i][0]))
            #    #config_header.write("#define {}_DL_MAX_{} {}\n".format(cosmo_name, i, lds[i][-1]))
            #    config_header.write("#define {}_DL_MIN_{} {}\n".format(k, i, lds[i][0]))
            #    config_header.write("#define {}_DL_MAX_{} {}\n".format(k, i, lds[i][-1]))

            #config_header.write("\n")

            ##########################################################
            #####  redshift to Luminosity distance###############
            ##########################################################
            coeffsZD = []
            for i in np.arange(int(num_segments[k])):
                popt, pcov = scipy.optimize.curve_fit(func,zs[i], lds[i])
                coeffsZD.append(popt)

            for i in np.arange(int(num_segments[k])):
                boundariesZ[k][i] = zs[i][-1]
                boundariesD[k][i] = lds[i][-1]
                for j in np.arange(int(deg[k])):
                    coeffsD[k][i][j] = coeffsZD[i][j]
                    coeffsZ[k][i][j] = coeffsDZ[i][j]
            #for i in np.arange(len(coeffsZD)):
            #    config_header.write("const double {}_COEFF_VEC_ZD_{}[{}] =  {{".format(cosmo_name, i, deg))
            #    config_header.write("{}".format(coeffsZD[i][0]))
            #    for j in np.arange(len(coeffsZD[0])-1):
            #        config_header.write(", {}".format(coeffsZD[i][j+1]))
            #    config_header.write("};\n")
            #
            #for i in np.arange(len(coeffsZD)):
            #    config_header.write("#define {}_Z_MIN_{} {}\n".format(cosmo_name, i, zs[i][0]))
            #    config_header.write("#define {}_Z_MAX_{} {}\n".format(cosmo_name, i, zs[i][-1]))

            ##########################################################
            ######## TESTING #####################
            ##########################################################
            fig, ax = plt.subplots(ncols=2, nrows=1, figsize=(8,5))
            for i in np.arange(len(coeffsDZ)):
                interpolated = func(lds[i], *coeffsDZ[i]) 
                ax[0].plot(lds[i]/1e3, (interpolated-zs[i])/zs[i])
                ax[0].set_title("Luminosity Distance to Redshift")
                ax[0].set_xlabel("Luminosity Distance (GPC)")
                ax[0].set_ylabel(r"Fractional Difference $(Z_{true} - Z_{interp})/Z_{true}$")
            for i in np.arange(len(coeffsZD)):
                interpolated = func(zs[i], *coeffsZD[i]) 
                ax[1].plot(zs[i], (interpolated-lds[i])/lds[i])
                ax[1].set_title("Redshift to Luminosity Distance")
                ax[1].set_xlabel("Redshift Z")
                ax[1].set_ylabel(r"Fractional Difference $(D_{true} - D_{interp})/D_{true}$")
            fig.suptitle(r"{} Residuals".format(cosmo_name))
            fig.tight_layout()
            fig.subplots_adjust(top=.88)
 
            #plt.show()
            plt.savefig(data_filetree+"/{}_residuals.png".format(cosmo_name))
            plt.close()

        config_header.write("const int num_cosmologies = {};".format(len(cosmologies)))
        #///////////////////////////////////////////////////////////////////
        config_header.write("const int num_segments[{}] = {{".format(len(cosmologies)))
        config_header.write(" {}".format(int(num_segments[0])))
        for k in np.arange(len(cosmologies)-1):
            config_header.write(", {}".format(int(num_segments[k+1]))) 
        config_header.write("};\n")
        #///////////////////////////////////////////////////////////////////
        #///////////////////////////////////////////////////////////////////
        config_header.write("const int interp_degree[{}] = {{".format(len(cosmologies)))
        config_header.write(" {}".format(int(deg[0])))
        for k in np.arange(len(cosmologies)-1):
            config_header.write(", {}".format(int(deg[k+1]))) 
        config_header.write("};\n")
        #///////////////////////////////////////////////////////////////////

        #///////////////////////////////////////////////////////////////////
        config_header.write("const double boundaries_Z[{}][{}] = {{".format(len(cosmologies), int(num_segments.max())+1))
        for k in np.arange(len(cosmologies)): 
            config_header.write("{")
            config_header.write("{}\n".format(boundariesZ[k][0]))
            for i in np.arange(int(num_segments[k])):
                config_header.write(", {}\n".format(boundariesZ[k][i]))
            if k == (len(cosmologies)-1):
                config_header.write("}")
            else:
                config_header.write("},")
        config_header.write("};\n")
        #///////////////////////////////////////////////////////////////////

        #///////////////////////////////////////////////////////////////////
        config_header.write("const double boundaries_D[{}][{}] = {{".format(len(cosmologies), int(num_segments.max())+1))
        for k in np.arange(len(cosmologies)): 
            config_header.write("{")
            config_header.write("{}\n".format(boundariesD[k][0]))
            for i in np.arange(int(num_segments[k])):
                config_header.write(", {}\n".format(boundariesD[k][i]))
            if k == (len(cosmologies)-1):
                config_header.write("}")
            else:
                config_header.write("},")
        config_header.write("};\n")
        #///////////////////////////////////////////////////////////////////

        #///////////////////////////////////////////////////////////////////
        config_header.write("const double COEFF_VEC_DZ[{}][{}][{}] =  {{".format(
                    len(cosmologies), int(num_segments.max()), int(deg.max()) ) )
        for k in np.arange(len(cosmologies)): 
            config_header.write("{")
            for i in np.arange(int(num_segments[k])):
                config_header.write("{{ {}\n".format(coeffsZ[k][i][0]))
                for j in np.arange(int(deg[k])-1):
                    config_header.write(", {}\n".format(coeffsZ[k][i][j+1]))
                if i == num_segments[k]-1:
                    config_header.write("}")
                else:
                    config_header.write("},")
            if k == (len(cosmologies)-1):
                config_header.write("}")
            else:
                config_header.write("},")
        config_header.write("};\n")
        #///////////////////////////////////////////////////////////////////

        #///////////////////////////////////////////////////////////////////
        config_header.write("const double COEFF_VEC_ZD[{}][{}][{}] =  {{".format(
                    len(cosmologies), int(num_segments.max()), int(deg.max()) ) )
        for k in np.arange(len(cosmologies)): 
            config_header.write("{")
            for i in np.arange(int(num_segments[k])):
                config_header.write("{{ {}\n".format(coeffsD[k][i][0]))
                for j in np.arange(int(deg[k])-1):
                    config_header.write(", {}\n".format(coeffsD[k][i][j+1]))
                if i == num_segments[k]-1:
                    config_header.write("}")
                else:
                    config_header.write("},")
            if k == (len(cosmologies)-1):
                config_header.write("}")
            else:
                config_header.write("},")
        config_header.write("};\n")
        #///////////////////////////////////////////////////////////////////
            
        config_header.write("\n")
        config_header.write("#endif")
        config_header.close()





        #for i in np.arange(len(coeffsDZ)):
        #    #config_header.write("#define {}_DL_MIN_{} {}\n".format(cosmo_name, i, lds[i][0]))
        #    #config_header.write("#define {}_DL_MAX_{} {}\n".format(cosmo_name, i, lds[i][-1]))
        #    config_header.write("#define {}_DL_MIN_{} {}\n".format(k, i, lds[i][0]))
        #    config_header.write("#define {}_DL_MAX_{} {}\n".format(k, i, lds[i][-1]))



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
