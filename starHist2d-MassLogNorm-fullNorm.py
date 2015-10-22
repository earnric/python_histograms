#
# 20 Oct 2015
# Rick Sarmento
#
# Purpose:
#  Reads star particle data and creates phase plots
#  Place histograms of x and y axis along axes
#  Uses pcolormesh norm=LogNorm(vmin=1,vmax=8)
#
# Method:
#  Main plot uses np.hist2d then takes log of result
#
# Revision history
#

# ##########################################################
# Generate colors for histogram bars based on height
# ##########################################################
def colorHistOnHeight(N, patches):
    fracs = np.log10(N.astype(float)) # Need to take log N since it is 
    # norm = mpl.colors.Normalize(fracs.min(), fracs.max())
    norm = mpl.colors.LogNorm(1.0, 14.0)
    # NOTE this color mapping is different from the one below.
    for thisfrac, thispatch in zip(fracs, patches):
        color = mpl.cm.jet(thisfrac)
        thispatch.set_facecolor(color)

    return


# ##########################################################
# Generate a density plot in log-log space, put histograms
# of the x and y values along axes
# ##########################################################
def genDensityPlot(x, y, mass, pf, z, minX, maxX, minY, maxY, filename, xaxislabel):
    nullfmt = NullFormatter() # Needed for empty labels... 

    # Plot location and size
    fig = plt.figure(figsize=(20, 20))
    ax2dhist = plt.axes(rect_2dhist)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # Fix any "log10(0)" points...
    x[x == np.inf] = 0.0
    y[y == np.inf] = 0.0
    y[y > 1.0] = 1.0  # Fix any minor numerical errors that could result in y>1

    # Compute comoving volume of sim
    volcm = 27.0 / (0.71)**3 # We're 'per h' so the volume is bigger than 27 Mpc^3

    # Bin data in log-space
    xrange = np.logspace(np.log10(minX), np.log10(maxX), xbins)
    yrange = np.logspace(np.log10(minY), np.log10(maxY), ybins)
    
    # Note axis order: y then x
    # H is the binned data... weighted by polluted mass of the sp's
    H, xedges, yedges = np.histogram2d(y, x, weights=mass * (1.0 - pf),  # We have log bins, so we take
                                       bins=(yrange, xrange))
    print("Raw H max, scaled by polluted mass %.2lf"%H.max())
    H = H / volcm # Normalize by comoving volume (now we're per Mpc)
    
    # size of each bin in x and y dimensions
    dx = np.diff(xrange)
    dy = np.diff(yrange) 
    area = dx[:,  None] * dy # compute the area of each bin using broadcasting
    H = H / area # Normalize by bin area
    H = np.log10(H)
    H = np.ma.array(H, mask=np.isnan(H))
    print("Normalized log H max %.2lf"%H.max())
    X, Y = np.meshgrid(xrange, yrange)  # Create a grid over the range of bins for the plot

    # Fix colors -- white for values of 1.0. 
    cmap = copy.copy(mpl.cm.jet)
    cmap.set_bad('w', 1.)  # w is color, for values of 1.0

    # Create a plot of the binned data, use a fixed lognormal scale
    #ma = ma/volcm # normalize by sim physical volume at z
    cax = (ax2dhist.pcolormesh(X, Y, H, cmap=cmap, norm=LogNorm(vmin=1,vmax=14)))

    # Setup the color bar
    cbar = fig.colorbar(cax, ticks=[1,2,4,10,14])
    cbar.ax.set_yticklabels(['1', '2', '4', '10', '14'], size=24)
    cbar.set_label('$log\, M_{sp, pol,\odot}\, / d\, ($ ' + xaxislabel
                   + " $) \, / d\, (log\, Z_{pri}/Z)\, /\, V)$ ", size=30)

    ax2dhist.tick_params(axis='x', labelsize=22)
    ax2dhist.tick_params(axis='y', labelsize=22)
    ax2dhist.set_xlabel(xaxislabel, size=30)
    ax2dhist.set_ylabel('$log\, Z_{pri}/Z$', size=30)

    ax2dhist.set_xlim([minX, maxX])
    ax2dhist.set_ylim([minY, maxY])
    ax2dhist.set_xscale('log')
    ax2dhist.set_yscale('log')
    ax2dhist.grid(color='0.75', linestyle=':', linewidth=2)

    # Generate the xy axes histograms
    xlims = ax2dhist.get_xlim()
    ylims = ax2dhist.get_ylim()

    # Note that even with log=True, the values in N are still masses (solar units)
    N, bins, patches = axHistx.hist(x, bins=xrange, log=True, weights=mass * (1.0 - pf))
    axHistx.set_xscale("log")
    colorHistOnHeight(N, patches) # Normalizes colors to our colorbar
    N, bins, patches = axHisty.hist(y, bins=yrange, log=True, weights=mass * (1.0 - pf),
                                     orientation='horizontal')
    axHisty.set_yscale('log')
    colorHistOnHeight(N, patches) # Normalizes colors to our colorbar

    # Setup format of the histograms
    axHistx.set_xlim(xlims)  # Match the x range on the horiz hist
    axHistx.set_ylim([ylims[-1], 1.e2])  # Constant range for all histograms
    axHistx.tick_params(labelsize=22)
    axHistx.yaxis.set_ticks([1e2, 1e4, 1e6, 1e8])
    axHistx.grid(color='0.75', linestyle=':', linewidth=2)

    axHisty.set_xlim([xlims[-1], 1.e2])  # We're rotated, so x axis is the value
    axHisty.set_ylim(ylims)  # Match the y range on the vert hist
    axHisty.tick_params(labelsize=22)
    axHisty.xaxis.set_ticks([1e2, 1e4, 1e6, 1e8])
    axHisty.grid(color='0.75', linestyle=':', linewidth=2)

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    if z[0] == '0': z = z[1:]
    axHistx.set_title('z=' + z, size=40)

    plt.savefig(filename + "-z_" + z + ".png", dpi=fig.dpi)
    #    plt.show()
    plt.close(fig)  # Release memory assoc'd with the plot
    return


# ##########################################################
# ##########################################################
##
## Main program
##
# ##########################################################
# ##########################################################
import matplotlib as mpl
import matplotlib.pyplot as plt
# import matplotlib.colors as colors # For the colored 1d histogram routine
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
import numpy as np
import copy as copy

files = [
    "18.00",
    "17.00",
    "16.00",
    "15.00",
    "14.00",
    "13.00",
    "12.00",
    "11.00",
    "10.00",
    "09.00",
    "08.50",
    "08.00",
    "07.50",
    "07.00",
    "06.50",
    "06.00",
    "05.50",
    "05.09"
]
# Plot parameters - global
left, width = 0.1, 0.63
bottom, height = 0.1, 0.63
bottom_h = left_h = left + width + 0.01

xbins = ybins = 100

rect_2dhist = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.15]
rect_histy = [left_h, bottom, 0.2, height]

prefix = "./"
# prefix="20Sep-BIG/"
for indx, z in enumerate(files):
    spZ = np.loadtxt(prefix + "spZ_" + z + ".txt", skiprows=1)
    spPZ = np.loadtxt(prefix + "spPZ_" + z + ".txt", skiprows=1)
    spPF = np.loadtxt(prefix + "spPPF_" + z + ".txt", skiprows=1)
    spMass = np.loadtxt(prefix + "spMass_" + z + ".txt", skiprows=1)

    print ("Generating phase diagram for z=%s" % z)
    minX = 1.e-8; maxX = 2.0
    minY = 1.e-4; maxY = 2.0
    genDensityPlot(spZ, spPZ / spZ, spMass, spPF, z, minX, maxX, minY, maxY,
                   "Z_PMassZ-MassHistLogNorm-Newmat", "$log\, Z_{\odot}$")
    minX = 1.e-5
    genDensityPlot(spZ / (1.0 - spPF), spPZ / spZ, spMass, spPF, z, minX, maxX, minY, maxY,
                   "Z_PMassZ_1-PGF-MassHistLogNorm-Newmat", "$log\, Z_{\odot}/f_{pol}$")
