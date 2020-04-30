import matplotlib.pyplot as plt
import matplotlib
import gwatpy.config as cf
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


