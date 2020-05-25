import matplotlib.pyplot as plt
import matplotlib
import gwatpy.config as cf
from matplotlib.colors import ListedColormap, to_rgba_array
import scipy.stats as st
import numpy as np
#from gwatpy.util import GWAT_ROOT_DIRECTORY_PY

def set():
    #fpath = "/Users/sperkins/git-repos/gw_analysis_tools/gwatpy/fonts/static/EBGaramond-Regular.ttf"
    #prop = matplotlib.font_manager.FontProperties(fname=fpath)
    
    #matplotlib.style.use(GWAT_ROOT_DIRECTORY_PY.decode('utf-8')+"gwatpy/gwat_rc")

    #matplotlib.style.use(GWAT_ROOT_DIRECTORY_PY.decode('utf-8')+"gwatpy/gwat_rc")
    matplotlib.style.use(cf.GWATPY_ROOT_DIRECTORY+"/gwat_rc")

    #matplotlib.rcParams["font.sans-serif"] = ["Apple Symbols","Microsoft Sans Serif","Tahoma"]
    matplotlib.rcParams["font.family"]="serif"
    #matplotlib.rcParams['font.family'] = prop.get_name()
    #matplotlib.rcParams['font.family'] = prop.get_name()
    #matplotlib.rcParams['font'] = prop.get_name()

def show_minor_grid():
    matplotlib.rcParams["axes.grid.which"]="both"
def show_major_grid():
    matplotlib.rcParams["axes.grid.which"]="major"

def _plot_kernel_density(x_data, y_data, ax, grid_points,weights, log_contours,  cmap, alpha, fill,verbose,quantiles):
    xmin = np.amin(x_data)
    xmax = np.amax(x_data)
    ymin = np.amin(y_data)
    ymax = np.amax(y_data)
    grid = 1j*grid_points
    #ymin = ax[x,y].get_xlim()[0]
    #ymax = ax[x,y].get_xlim()[1]
    #xmax = ax[x,y].get_ylim()[1]
    #xmin = ax[x,y].get_ylim()[0]
    xx, yy = np.mgrid[xmin:xmax:grid,ymin:ymax:grid]
    positions = np.vstack([xx.ravel(),yy.ravel()])
    values = np.vstack([x_data,y_data])
    kernel = st.gaussian_kde(values,weights=weights)
    f = np.reshape(kernel(positions).T,xx.shape)
    #cfset = ax[x,y].contourf(yy,xx,f,cmap=cmap,alpha=alpha)
    levels  = np.quantile(f,quantiles)
    levels = np.append(levels,np.quantile(f,1)) 
    if verbose:
        print("Levels: ",levels)
    if fill:
        if verbose:
            print("Plotting filled contours")
        cfset = ax.contourf(xx,yy,f,levels=levels,cmap=cmap,alpha=alpha)
    else:
        if verbose:
            print("Plotting contours only")
        cfset = ax.contour(xx,yy,f,levels=levels,cmap=cmap)

def gwatpy_corner(data, data2 = None,fig=None, density=True, bins = None,bins2 = None,weights = None,weights2 = None,alpha = None,alpha2 = None,figsize=None,labels = None,color = "black",color2 = "blue",cmap=None,cmap2=None,kernel_density=False,grid_points=10,log_contours=False, verbose=False,fill = True,quantiles=[.5,.9]):
    #N = 256
    #rgb1 = to_rgba_array("white")
    #rgb = to_rgba_array(color)
    #rgb2 = to_rgba_array("black")
    #print(rgb)
    #vals = np.ones((N,4))
    #vals[:, 0] = np.linspace(90/256, 1, N)
    #vals[:, 1] = np.linspace(39/256, 1, N)
    #vals[:, 2] = np.linspace(41/256, 1, N)
    #newcmp = ListedColormap(rgb)
    #ccmap = newcmp

    alpha_diag = 1 
    alpha_offdiag = 1 
    if isinstance(alpha,list):
        alpha_diag = alpha[0] 
        alpha_offdiag = alpha[1] 
    elif isinstance(alpha,float):
        alpha_diag = alpha 
        alpha_offdiag = alpha 

    rows = len(data[0]) 
    cols = len(data[0]) 

    if fig is None :
        if verbose:
            print("Creating new figure")
        fig, ax = plt.subplots(nrows=rows,ncols=cols,sharex="col",figsize=figsize)
        if verbose:
            print("Linking axes")
        for y in np.arange(cols):
            if y != 0 and y != 1 :
                ax[0,1].get_shared_y_axes().join(ax[0,1],ax[0,y])
        for x in np.arange(1,rows):
            for y in np.arange(cols):
                if y !=0 and y != x :
                    ax[x,0].get_shared_y_axes().join(ax[x,0],ax[x,y])
        plt.tight_layout()
        fig.subplots_adjust(hspace=0.0,wspace=0.0)
    else:
        if verbose:
            print("Retrieving axes from figure")
        ax = fig.axes
        ax = np.reshape(ax,(rows,cols))

    for x in np.arange(rows):
        for y in np.arange(cols):
            if verbose:
                print("Plotting: ",x,y)

            if y> x and data2 is None:
                ax[x,y].axis("off")
                continue

            if x == y:
                ax[x,y].hist(data[:,x], density=density,weights=weights,bins=bins,color=color,alpha=alpha)
                if data2 is not None:
                    ax[x,y].hist(data2[:,x], density=density,weights=weights2,bins=bins2,color=color2,alpha=alpha2)
                #if labels is not None and data2 is None:
                #    ax[x,y].set_title(labels[x])

            elif data2 is None or y<x:
                if not kernel_density:
                    ax[x,y].hexbin(data[:,y],data[:,x], mincnt=1, C = weights,reduce_C_function = np.sum,alpha=alpha,cmap =cmap ,gridsize=grid_points)
                else:
                    _plot_kernel_density(data[:,y], data[:,x], ax[x,y], grid_points,weights, log_contours,  cmap, alpha, fill,verbose,quantiles)

            elif data2 is not None:
                if not kernel_density:
                    ax[x,y].hexbin(data2[:,y],data2[:,x], mincnt=1, C = weights2,reduce_C_function = np.sum,alpha=alpha2,cmap =cmap2 ,gridsize=grid_points)
                else:
                    _plot_kernel_density(data2[:,y], data2[:,x], ax[x,y], grid_points,weights2, log_contours,  cmap2, alpha2, fill,verbose,quantiles)

                #if labels is not None and x == 0:
                #    ax[x,y].set_ylabel(labels[x])
            if y != 0 or x==y:
                ax[x,y].set_yticklabels([])
    for x in np.arange(rows):
        for y in np.arange(cols):
            if y == x:
                if labels is not None and data2 is None:
                    ax[x,y].set_title(labels[x])
            if data2 is not None :
                if x == 0:
                    if labels is not None:
                        ax[x,y].set_title(labels[y])
                    ax[x,y].xaxis.tick_top()
    return fig
    
if __name__=="__main__":
    test = np.random.multivariate_normal(mean=np.asarray([-1,150,1540]),cov=np.asarray([[1.,.5,.2],[.5,1.,.2],[.1,.1,5]]),size=10000)
    test2= np.random.multivariate_normal(mean=np.asarray([0,145,1550]),cov=np.asarray([[1.,.2,.2],[.2,1.,.5],[.1,.1,5]]),size=1000)
    #weights = [ 1 if x[0]>0 else 0 for x in test2]
    weights = None
    alpha = 0.9
    labels = ["x",'y',"z"]
    gridsize = 20
    #quantiles = np.logspace(np.log10(.1),np.log10(.9),10)
    quantiles = np.linspace(.1,.99,100)
    print(quantiles)
    fig = gwatpy_corner(test,bins=50,alpha = alpha,figsize=[10,5] ,labels=labels,color="black",cmap = 'Greys',kernel_density=True,log_contours=False,verbose=True,grid_points=gridsize,fill=True,quantiles = quantiles)
    alpha = 0.7
    gwatpy_corner(test2,bins=50,alpha = alpha,figsize=[10,5] ,labels=labels,color="red",cmap = 'Reds',fig=fig,kernel_density=True,weights=weights,verbose=True,grid_points=gridsize,fill=True,quantiles=quantiles)
    plt.savefig("test.pdf")

    fig = gwatpy_corner(test,test2,bins=50,bins2=50,alpha = alpha,alpha2 = alpha, color2="red",cmap2="Reds",figsize=[10,5] ,labels=labels,color="black",cmap = 'Greys',kernel_density=True,log_contours=True,verbose=True,grid_points=gridsize,fill=True,quantiles=quantiles)
    plt.savefig("test2.pdf")
