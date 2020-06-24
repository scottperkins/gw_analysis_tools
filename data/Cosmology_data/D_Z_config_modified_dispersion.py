##########################################################
#This code writes a header file to store numerical fit 
#data to calculate redshift from modified luminosity distance 
#quicker and more efficiently. To do so, I use the 
#scipy.optimize.curve_fit function to fit coefficients 
#to the function defined below (power series in 1/2 power 
#increments) which is then written out to a header file
#which is in turn included in the util.cpp file.

#This file shouldn't need to be run once the header file is 
#created once, but its available to update cosmology information
#or to do cosmological studies

#This is run only for Planck15 cosmology, with a range of 
#dispersion powers
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
from functools import partial
from scipy import optimize
c = 299792458.
MPC_SEC = 3.085677581491367278913937957796471611e22/c



#Function being fit to the data -- power series in steps of 1/2 power
#def func(x, l,a, b ,c, d, e, f, g, h, i, j ,k):
#    return  l+(a*np.sqrt(x)  + b*x + c*x**(3./2.) + d*x*x + e*x**(5./2.) + f*x*x*x  +
#            g*x**(7./2.) + h*x*x*x*x + i*x**(9./2.) + j*x**(10./2.) + k*x**(11./2.) )
#Function being fit to the data -- power series in steps of 1/2 power
#def func(x, a,b,c,d,e,f,g,h,i,j,l,m,n,o,p,q,r,s,t,u,v):
def func(x, d,e,f,g,h,i,j,l,m,n,o,p,q,r,s,t,u):
#def func(x, e,f,g,h,i,j,l,m,n,o,p,q,r):
    #parameters = [a,b,c,d,e,f,g,h,i,j,l,m,n,o,p,q,r,s,t,u]
    #returnval = np.sum(np.asarray( [parameters[y]*x**(-5+y*.5) for y in range(len(parameters))] ))
    #returnval = a*x**(-5)
    #returnval = b*x**(-4.5)
    #returnval = c*x**(-4.0)
    returnval = d*x**(-3.5)
    returnval += e*x**(-3.)
    returnval += f*x**(-2.5)
    returnval += g*x**(-2.0)
    returnval += h*x**(-1.5)
    returnval += i*x**(-1.0)
    returnval += j*x**(-0.5)
    returnval += l
    returnval += m*x**(0.5)
    returnval += n*x**(1.0)
    returnval += o*x**(1.5)
    returnval += p*x**(2.0)
    returnval += q*x**(2.5)
    returnval += r*x**(3.0)
    returnval += s*x**(3.5)
    returnval += t*x**(4.0)
    returnval += u*x**(4.5)
    #returnval += v*x**(5)
    return returnval
            
if __name__=="__main__":
    testing_cosmo = cosmos.FlatLambdaCDM(H0=70,Om0=0.3)
    data_filetree = os.path.dirname(os.path.realpath(__file__))
    #include_filetree = data_filetree[:-4] + "include"
    include_filetree = data_filetree[:-19] + "include"
    #config_header = open(include_filetree+'/D_Z_Config.h','w')
    with open(include_filetree+'/gwat/D_Z_Config_modified_dispersion.h','w') as config_header:
        alphas = [0,0.5, 1,1.5, 2,2.5,3,3.5,4]
        dispersion_names = ["alpha_0", "alpha_0.5", "alpha_1","alpha_1.5","alpha_2","alpha_2.5","alpha_3","alpha_3.5","alpha_4"]
        #alphas = [0]
        #dispersion_names = [ "alpha_0"]
        deg = np.ones((len(alphas),))*17 #Number of parameters in the fitting function
        num_segments = np.ones((len(alphas),))*3 #This is how many chunks the interpolation range is broken up into 
        pts_seg = []
        for i in np.arange(len(alphas)):
            pts_seg.append([])
            for j in np.arange(num_segments[i]):
                pts_seg[-1] .append(1e3)#points in each segment
                #pts_seg[-1] .append(1e2)#points in each segment
        boundariesZ_base = []
        for i in np.arange(len(alphas)):
            boundariesZ_base .append(np.logspace(np.log10(1e-6),np.log10(20), num_segments[i]+1))#points in each segment
        #boundariesZ_base = np.logspace(np.log10(1e-6),np.log10(20), num_segments+1) #Written based on chunks of Z, then boundaries of DL will be dependent on cosmology

        config_header.write("/*! \\file \n \* Header file for the cosmological interpolation \n */ \n")
        config_header.write("// This header file contains the coefficients for the interpolation \n")
        config_header.write("// function for the conversion from redshift to modified dispersion distance and vice versa.  \n")
        config_header.write("// Adding cosmologies and adding more segments for accuracy should be solely done in the   \n")
        config_header.write("// python file D_Z_config_modified_dispersion.py file. Add the name to be used in the entire package,\n")
        config_header.write("// change the segments, and points per segement as needed, and rerun the python script to interpolate.\n")
        config_header.write("// The output is written to the header file, and the option is now available anywhere in the C++ code.\n")
        config_header.write("#ifndef D_Z_CONFIG_MODIFIED_DISPERSION_H\n")
        config_header.write("#define D_Z_CONFIG_MODIFIED_DISPERSION_H\n")
        config_header.write("\n")
        config_header.write("// I am aware this isn't very human-readable -- python has issues with line sizes, so the lines had to be truncated in some fashion\n")

        config_header.write("const char * MD_alphas_names[{}] =  {{".format(len(dispersion_names)))
        config_header.write("\"{}\"".format(dispersion_names[0]))
        for i in np.arange(len(dispersion_names)-1):
            config_header.write(", \"{}\"".format(dispersion_names[i+1]))
        config_header.write("};\n")
        config_header.write("const double MD_alphas[{}] =  {{".format(len(alphas)))
        config_header.write("{}".format(alphas[0]))
        for i in np.arange(len(alphas)-1):
            config_header.write(", {}".format(alphas[i+1]))
        config_header.write("};\n")
        
        boundariesZ =[]
        boundariesD =[]
        coeffsZ =[]
        coeffsD =[]
        for i in np.arange(len(alphas)):
            #boundariesZ.append([])
            #boundariesD.append([])
            #coeffsZ.append([])
            #coeffsD.append([])
            boundariesZ.append(boundariesZ_base[i])
            boundariesD.append(np.zeros((int(num_segments[i])+1,)))
            coeffsZ.append(np.zeros((int(num_segments[i]),int(deg[i]))))
            coeffsD.append(np.zeros((int(num_segments[i]),int(deg[i]))))

        for k in np.arange(len(alphas)):
            cosmo = cosmos.Planck15
            dispersion_name = dispersion_names[k]
            local_alpha = alphas[k]
            def H(z):
                return cosmo.H(z).to('Hz').value
            def helper(alpha,y):
                return (1+y)**(1-alpha)*quad(lambda x: (1+x)**(alpha-2) /H(x),0,y)[0]
            #def helperold(y):
            #    return Distance(y,unit=u.Mpc).compute_z(cosmology = cosmo)
            #def helper2(zfinal):
            #    return quad(lambda x: 1/((1+x)**2*H(x)),0,zfinal )[0]
            pool = mp.Pool(processes=mp.cpu_count())
            helper_wrapper = partial(helper,local_alpha)

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
                lds.append(np.asarray(pool.map(helper_wrapper, zs[i]))/MPC_SEC)
            ##########################################################
            ##### redshift to luminosity distance ###############
            ##########################################################
            #coeffsDZ = []
            #for v in np.arange(len(zs)):
            #    #popt, pcov = scipy.optimize.curve_fit(lambda x, a,b,c,d,e,f,g,h,i,j,l,m,n,o,p,q,r,s,t,w,aa: func(x,a,b,c,d,e,f,g,h,i,j,l,m,n,o,p,q,r,s,t,w,aa),lds[v], zs[v])
            #    popt, pcov = scipy.optimize.curve_fit(func,lds[v], zs[v])
            #    coeffsDZ.append(popt)

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
            for v in np.arange(int(num_segments[k])):
                #popt, pcov = scipy.optimize.curve_fit(lambda x, a,b,c,d,e,f,g,h,i,j,l,m,n,o,p,q,r,s,t,w,aa: func(x,a,b,c,d,e,f,g,h,i,j,l,m,n,o,p,q,r,s,t,w,aa),zs[v], lds[v])
                popt, pcov = scipy.optimize.curve_fit(func,zs[v], lds[v])
                #popt, pcov = scipy.optimize.curve_fit(func,zs[i], lds[i])
                coeffsZD.append(popt)

            boundariesZ[k][0] = zs[0][0]
            boundariesD[k][0] = lds[0][0]
            for i in np.arange(int(num_segments[k])):
                boundariesZ[k][i+1] = zs[i][-1]
                boundariesD[k][i+1] = lds[i][-1]
                for j in np.arange(int(deg[k])):
                    coeffsD[k][i][j] = coeffsZD[i][j]
                    #coeffsZ[k][i][j] = coeffsDZ[i][j]
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
            #for i in np.arange(len(coeffsDZ)):
            #    zs_test = np.logspace(np.log10(boundariesZ[k][i]),np.log10(boundariesZ[k][i+1]),100)
            #    lds_test = np.asarray(pool.map(helper_wrapper,zs_test))/MPC_SEC
            #    #print(zs_test,lds_test)
            #    interpolated = func(lds_test, *coeffsDZ[i]) 
            #    #print(interpolated)
            #    #print((interpolated-zs_test)/zs_test)
            #    #ax[0].plot(lds_test/1e3, abs(interpolated-zs_test)/zs_test,color="blue")
            #    ax[0].set_title("Modified Distance to Redshift")
            #    ax[0].set_xlabel("Modified Distance (GPC)")
            #    ax[0].set_ylabel(r"Fractional Difference $(Z_{true} - Z_{interp})/Z_{true}$")
            for i in np.arange(len(coeffsZD)):
                zs_test = np.logspace(np.log10(boundariesZ[k][i]),np.log10(boundariesZ[k][i+1]),100)
                lds_test = np.asarray(pool.map(helper_wrapper,zs_test))/MPC_SEC
                interpolated = np.asarray([func(zs_test[l], *coeffsZD[i]) for l in range(len(zs_test))] )
                #[print(coeffsZD[i+1][m],alphas[k]-20./4+m*.5) for m in range(len(coeffsZD[i+1]))]
                ax[1].loglog(zs_test, abs(interpolated-lds_test)/lds_test,color="blue")
                #print(abs(interpolated-lds_test)/lds_test)
                #ax[1].loglog(zs_test,interpolated,color="red")
                #ax[1].loglog(zs_test,lds_test,color="blue")
                #ax[1].set_yscale("log")
                ax[1].set_title("Redshift to Modified Distance")
                ax[1].set_xlabel("Redshift Z")
                ax[1].set_ylabel(r"Fractional Difference $(D_{true} - D_{interp})/D_{true}$")
            fig.suptitle(r"{} Residuals".format(dispersion_names[k]))
            fig.tight_layout()
            fig.subplots_adjust(top=.88)
 
            #plt.show()
            plt.savefig(data_filetree+"/{}_residuals.png".format(dispersion_names[k]))
            plt.close()

        config_header.write("const int num_MD_alphas = {};".format(len(alphas)))
        #///////////////////////////////////////////////////////////////////
        config_header.write("const int num_MD_segments[{}] = {{".format(len(alphas)))
        config_header.write(" {}".format(int(num_segments[0])))
        for k in np.arange(len(alphas)-1):
            config_header.write(", {}".format(int(num_segments[k+1]))) 
        config_header.write("};\n")
        #///////////////////////////////////////////////////////////////////
        #///////////////////////////////////////////////////////////////////
        config_header.write("const int interp_MD_degree[{}] = {{".format(len(alphas)))
        config_header.write(" {}".format(int(deg[0])))
        for k in np.arange(len(alphas)-1):
            config_header.write(", {}".format(int(deg[k+1]))) 
        config_header.write("};\n")
        #///////////////////////////////////////////////////////////////////

        #///////////////////////////////////////////////////////////////////
        config_header.write("const double MD_boundaries_Z[{}][{}] = {{".format(len(alphas), int(num_segments.max())+1))
        for k in np.arange(len(alphas)): 
            config_header.write("{")
            config_header.write("{}\n".format(boundariesZ[k][0]))
            for i in np.arange(int(num_segments[k])):
                config_header.write(", {}\n".format(boundariesZ[k][i+1]))
            if k == (len(alphas)-1):
                config_header.write("}")
            else:
                config_header.write("},")
        config_header.write("};\n")
        #///////////////////////////////////////////////////////////////////

        #///////////////////////////////////////////////////////////////////
        config_header.write("const double MD_boundaries_D[{}][{}] = {{".format(len(alphas), int(num_segments.max())+1))
        for k in np.arange(len(alphas)): 
            config_header.write("{")
            config_header.write("{}\n".format(boundariesD[k][0]))
            for i in np.arange(int(num_segments[k])):
                config_header.write(", {}\n".format(boundariesD[k][i+1]))
            if k == (len(alphas)-1):
                config_header.write("}")
            else:
                config_header.write("},")
        config_header.write("};\n")
        #///////////////////////////////////////////////////////////////////

        #///////////////////////////////////////////////////////////////////
        #config_header.write("const double MD_COEFF_VEC_DZ[{}][{}][{}] =  {{".format(
        #            len(alphas), int(num_segments.max()), int(deg.max()) ) )
        #for k in np.arange(len(alphas)): 
        #    config_header.write("{")
        #    for i in np.arange(int(num_segments[k])):
        #        config_header.write("{{ {}\n".format(coeffsZ[k][i][0]))
        #        for j in np.arange(int(deg[k])-1):
        #            config_header.write(", {}\n".format(coeffsZ[k][i][j+1]))
        #        if i == num_segments[k]-1:
        #            config_header.write("}")
        #        else:
        #            config_header.write("},")
        #    if k == (len(alphas)-1):
        #        config_header.write("}")
        #    else:
        #        config_header.write("},")
        #config_header.write("};\n")
        #///////////////////////////////////////////////////////////////////

        #///////////////////////////////////////////////////////////////////
        config_header.write("const double MD_COEFF_VEC_ZD[{}][{}][{}] =  {{".format(
                    len(alphas), int(num_segments.max()), int(deg.max()) ) )
        for k in np.arange(len(alphas)): 
            config_header.write("{")
            for i in np.arange(int(num_segments[k])):
                config_header.write("{{ {}\n".format(coeffsD[k][i][0]))
                for j in np.arange(int(deg[k])-1):
                    config_header.write(", {}\n".format(coeffsD[k][i][j+1]))
                if i == num_segments[k]-1:
                    config_header.write("}")
                else:
                    config_header.write("},")
            if k == (len(alphas)-1):
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
