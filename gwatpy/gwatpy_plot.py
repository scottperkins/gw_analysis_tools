import matplotlib.pyplot as plt
import matplotlib

def set():
    #fpath = "/Users/sperkins/git-repos/gw_analysis_tools/gwatpy/fonts/static/EBGaramond-Regular.ttf"
    #prop = matplotlib.font_manager.FontProperties(fname=fpath)
    
    matplotlib.style.use("/Users/sperkins/git-repos/gw_analysis_tools/gwatpy/gwat_rc")
    #matplotlib.rcParams["font.sans-serif"] = ["Apple Symbols","Microsoft Sans Serif","Tahoma"]
    matplotlib.rcParams["font.family"]="serif"
    #matplotlib.rcParams['font.family'] = prop.get_name()
    #matplotlib.rcParams['font.family'] = prop.get_name()
    #matplotlib.rcParams['font'] = prop.get_name()

def show_minor_grid():
    matplotlib.rcParams["axes.grid.which"]="both"
def show_major_grid():
    matplotlib.rcParams["axes.grid.which"]="major"


