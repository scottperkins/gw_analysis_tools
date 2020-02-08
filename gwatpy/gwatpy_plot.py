import matplotlib.pyplot as plt
import matplotlib
matplotlib.style.use("/Users/sperkins/git-repos/gw_analysis_tools/gwatpy/gwat_rc")
matplotlib.rcParams["font.sans-serif"] = ["Apple Symbols","Microsoft Sans Serif","Tahoma"]

def show_minor_grid():
    matplotlib.rcParams["axes.grid.which"]="both"
def show_major_grid():
    matplotlib.rcParams["axes.grid.which"]="major"


