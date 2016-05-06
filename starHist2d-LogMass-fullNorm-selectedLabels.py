#
# 20 Oct 2015
# Rick Sarmento
#
# Purpose:
#  Reads star particle data and creates prob density plots
#  Places histograms of x and y axis along axes
#  The data is normalized by bin size (area for hist2d,
#  width for hist) and comoving sim volume (in Mpc). The
#  data is also plotted on log scales.
#
# Revision history
#  29 Oct 2015 - Initial version
#  04 Nov 2015 - Fixed normalization

# File globals...
# Max color range value, log
global minCV,maxCV
# Max value for axis histograms
global histMax
# Plotting range limits, log
global minY, maxY
global minX, maxX

# ##########################################################
# Generate colors for histogram bars based on a normalized
# height. Normalize by bin width and comoving vol
# Method:
#  Take log of the histogram values (weighted counts)..
#  Create a LogNorm mapping between 1->max
#  Use the norm to map scalar values between 1 & max to rgb
# ##########################################################
def colorHistOnHeight(N, bins, patches, cmvol):
    cleanN = np.ma.masked_where(N == 0.0, N)
    widths = np.diff(np.log10(bins))
    fracs  = np.log10(cleanN/widths/cmvol)

    # normalize colors to the top of our scale
    norm   = mpl.colors.LogNorm(vmin=1, vmax=maxCV) 
    sm     = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.jet)
    sm.set_clim([minCV, maxCV]) # Force to use the whole range
    for thisfrac, thispatch in zip(fracs, patches):
        color = sm.to_rgba(thisfrac)
        thispatch.set_facecolor(color)
    return 

# ##########################################################
# Normalizes the histogram's bar height by the bin width
# and comoving volume of the simulation
# ##########################################################
def normBarHeight(bins, patches, cmvol, rotated=False):
    widths = np.diff(np.log10(bins))
    #print ("patches %i, bins %i"%(len(patches),len(widths)))
    for item,dbx in zip(patches,widths):
        #print ("Starting height: %.2f bin width: %.2e"%(item.get_height(),dbx))
        if not rotated:
            item.set_height(item.get_height()/dbx/cmvol)
        else:
            item.set_width(item.get_width()/dbx/cmvol)
        #print ("Ending width: %.2f"%item.get_width())
    return

# ##########################################################
# Format tick labels using their exponent only... 
# ##########################################################
def formatter(x,pos):
    return format(np.log10(x), '.0f')

# ##########################################################
# Generate a density plot in log-log space, put histograms
# of the x and y values along axes
# ##########################################################
def genDensityPlot(x, y, mass, pf, z, filename, xaxislabel, yaxislabel, normByPMass=True):
    import matplotlib.patches as patches
    from matplotlib.ticker import FuncFormatter

    zint = int(float(z))

    # If we had Z/f_pol = inf or nan it means with had f_pol = 0 or both = 0
    x[np.isinf(x)] = 10**minX
    x[np.isnan(x)] = 10**minX

    y[np.isinf(y)] = 10**minY
    y[np.isnan(y)] = 10**minY

    x[x<=10**minX] = 10**minX # catch the 0's... 
    y[y<=10**minY] = 10**minY
    
    # Fix any incoming bogus points...
    y[y > 1.0] = 1.0  # spPZ / spZ - Fix any numerical errs that could result in y>1

    # ################################################################
    # For weighting of PopIII Stars, ensure fpol = 1 since we've
    # already fixed x and y
    pf[pf > 1.0 - 1e-6] = 0.0
    # ################################################################
    
    ## print("genDensityPlot after X Inf #",len(x[np.isinf(x)]))
    ## print("genDensityPlot after X Nan #",len(x[np.isnan(x)]))
    ## print("genDensityPlot after X 10**minX #",len(x[x==10**minX]),minX)
    ## print("genDensityPlot after X 0 #",len(x[x==0]),minX)
    ## print("genDensityPlot after y Inf #",len(y[np.isinf(y)]))
    ## print("genDensityPlot after y Nan #",len(y[np.isnan(y)]))
    ## print("genDensityPlot after y 10**minY #",len(y[y == 10**minY]),minY)

    labelsize = 42
    nullfmt = NullFormatter()

    #FuncFormatter class instance defined from the function above
    custom_formatter = FuncFormatter(formatter)

    # Plot location and size
    fig = plt.figure(figsize=(20, 20))
    ax2dhist = plt.axes(rect_2dhist)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)

    # Compute comoving volume of sim
    cmvol = 27.0 / (0.71)**3 # We're 'per h' so the volume is bigger than 27 Mpc^3
    
    # Bin data in log-space
    # logsace expects min/max to be exponents... returns a log-spaced range
    xrange = np.logspace(minX-0.1,maxX,xbins)
    yrange = np.logspace(minY-0.1,maxY,ybins)

    # Note axis order: y then x
    # H is the binned data... weighted by polluted mass of the sp's
    # If we're looking at x = log Z, don't weight by mass * f_p... just mass!
    if normByPMass:
        H, xedges, yedges = np.histogram2d(y, x, weights=mass * (1.0 - pf), 
                                            bins=(yrange,xrange))
    else:
        H, xedges, yedges = np.histogram2d(y, x, weights=mass, 
                                            bins=(yrange,xrange))
    # Normalize the histogram data... 
    H = H / cmvol # Normalize by comoving volume (now we're per Mpc)
    
    # size of each bin in x and y dimensions, in log
    dx = np.diff(np.log10(xrange))
    dy = np.diff(np.log10(yrange))
    area = dx[:, None] * dy # compute the area of each bin using broadcasting

    H = H / area # Normalize by bin area
    H = np.log10(np.ma.masked_where(H == 0.0,H))
    H = np.ma.masked_invalid(H)

    X, Y = np.meshgrid(xrange, yrange)  # Create a grid over the range of bins for the plot

    # Fix colors -- white for values of 1.0. 
    cmap = copy.copy(mpl.cm.jet)
    cmap.set_bad('w', 1.)  # w is color, for values of 1.0

    # Create a plot of the binned data, use a fixed lognormal scale
    cax = (ax2dhist.pcolormesh(X, Y, H, cmap=cmap, norm=LogNorm(vmin=minCV,vmax=maxCV)))

    # Setup the color bar
    cbarticks = [1,2,4,6,8,maxCV]
    cbar = fig.colorbar(cax, spacing='uniform', ticks=[1,2,4,6,8,maxCV])
    cbar.ax.set_yticklabels(cbarticks, size=labelsize)

    clabel = r'${\rm log}\, \left(M/M_{\odot}\, / d\, \left(' + xaxislabel + r'\right) \, / d\, \left(' + yaxislabel + r'\right)\, /\, {\rm Mpc}^{3}\right)$'

    ## print("Color label: ",clabel)

    cbar.set_label(clabel, size=labelsize - 4, y=0.6)

    ax2dhist.tick_params(axis='x', labelsize=labelsize)
    ax2dhist.tick_params(axis='y', labelsize=labelsize)
    ax2dhist.set_xlabel('$'+xaxislabel+'$', size=labelsize)
    ax2dhist.set_ylabel('$'+yaxislabel+'$', size=labelsize)

    ax2dhist.set_xlim([10**(minX-0.1),10**maxX])
    ax2dhist.set_ylim([10**(minY-0.1),10**maxY])
    ax2dhist.set_xscale('log')
    ax2dhist.set_yscale('log')
    ax2dhist.grid(color='0.75', linestyle=':', linewidth=2)

    ax2dhist.xaxis.set_major_formatter(custom_formatter)
    ax2dhist.yaxis.set_major_formatter(custom_formatter)
    xpos = 1e-5;ypos = 10**-3.77
    ax2dhist.text(xpos,ypos,'z = %d'%zint, fontsize=28)
    # Generate the xy axes histograms
    xlims = ax2dhist.get_xlim()
    ylims = ax2dhist.get_ylim()

    ##########################################################
    # Create the histograms for the x and y axes
    # Note that even with log=True, the array N is NOT log of the
    # weighted counts. 
    # Normalize the histogram heights by the bin width & comoving volume
    if normByPMass:
        N, bins, patches = axHistx.hist(x, bins=xrange, log=True, weights=mass * (1.0 - pf))
        normBarHeight(bins, patches, cmvol)
    else:
        N, bins, patches = axHistx.hist(x, bins=xrange, log=True, weights=mass)
        normBarHeight(bins, patches, cmvol)

    axHistx.set_xscale("log")
    colorHistOnHeight(N, bins, patches, cmvol)

    ##########################################################
    # Create the histograms for the x and y axes
    if normByPMass:
        N, bins, patches = axHisty.hist(y, bins=yrange, log=True, weights=mass * (1.0 - pf),
                                        orientation='horizontal')
        normBarHeight(bins, patches, cmvol, rotated=True)

    else:
        N, bins, patches = axHisty.hist(y, bins=yrange, log=True, weights=mass,
                                        orientation='horizontal')        
        normBarHeight(bins, patches, cmvol, rotated=True)
        xpos = 1e-5;ypos = 10**-3.95
        ax2dhist.text(xpos,ypos,'Uncorrected', fontsize=28)

    axHisty.set_yscale('log')
    colorHistOnHeight(N, bins, patches, cmvol) # Normalizes colors to our colorbar

    # Setup format of the histograms
    axHistx.set_xlim(ax2dhist.get_xlim())  # Match the x range on the horiz hist
    axHistx.set_ylim([100.0,10.0**histMax])     # Constant range for all histograms
    axHistx.tick_params(labelsize=labelsize)
    if histMax == 10:
        axHistx.yaxis.set_ticks([1e3, 1e6, 1e9])
    else:
        axHistx.yaxis.set_ticks([1e3, 1e6, 1e9, 1e12])
    axHistx.grid(color='0.75', linestyle=':', linewidth=2)
    axHistx.yaxis.set_major_formatter(custom_formatter) # Use integers for ticks

    axHisty.set_xlim([100.0,10.0**histMax])     # We're rotated, so x axis is the value
    axHisty.set_ylim([10**(minY-0.1),10**maxY])  # Match the y range on the vert hist
    axHisty.tick_params(labelsize=labelsize)
    if histMax == 10:
        axHisty.xaxis.set_ticks([1e3, 1e6, 1e9])
    else:
        axHisty.xaxis.set_ticks([1e3, 1e6, 1e9, 1e12])
    axHisty.grid(color='0.75', linestyle=':', linewidth=2)
    axHisty.xaxis.set_major_formatter(custom_formatter) # Use integers for ticks

    # no labels
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    ## titlez = z
    ## if z[0] == '0': titlez = z[1:]
    ## axHistx.set_title('z=' + titlez, size=40)

    plt.savefig(filename + '-z_%s.pdf'%z, dpi=fig.dpi,bbox_inches='tight')
    # plt.show()
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
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
import numpy as np
import copy as copy
import gc

files = [
    ## "18.00",
    ##"17.00",
    "16.00",
    ## "15.00",
    ## "14.00",
    ## "13.00",
    "12.00",
    ## "11.00",
    "10.00",
    ## "09.00",
    ## "08.50",
    "08.00",
    ## "07.50",
    ## "07.00",
    ## "06.50",
    "06.00",
    ## "05.50",
    "05.00"
]

# Plot area sizes...
left, width = 0.1, 0.61
bottom, height = 0.1, (0.61 * 4.0/5.1)
bottom_h = left + height + 0.0
left_h = left + width + 0.0

rect_2dhist = [left, bottom, width, height]
rect_histx = [left, bottom_h, width, 0.11]
rect_histy = [left_h, bottom, 0.11, height]

# Number of bins
xbins = ybins = 100

# Process files and generate plots
prefix = "./"
minCV = 1; maxCV = 10
for indx, z in enumerate(files):
    spZ = np.loadtxt(prefix + "spZ_" + z + ".txt", skiprows=1)
    spPZ = np.loadtxt(prefix + "spPZ_" + z + ".txt", skiprows=1)
    spPF = np.loadtxt(prefix + "spPPF_" + z + ".txt", skiprows=1)
    spMass = np.loadtxt(prefix + "spMass_" + z + ".txt", skiprows=1)

    # RS - 24 Feb, don't let polluted fraction be less than 1:1e6
    print("\nFixing really small polluted fractions...")
    print("P > 1-1e-6: ", len(spPF[(spPF > 1.0-1e-10) & (spPF < 1.0)]))
    spPF[spPF > 1.0-1e-6] = 1.0

    print ("\nGenerating phase diagram for z=%s\n" % z)
    print("Z < 5.1", len(spZ[spZ<1e-5]))
    print("Z == 0", len(spZ[spZ==0]))
    f_pol = (1.0 - spPF)  # The polluted fraction
    print("f_pol < 1e-6: ", len(f_pol[f_pol < 1e-6]))
    print("f_pol == 0: ", len(f_pol[f_pol == 0.0]))

    # Put 0 metallicity stars just below critical... 
    spZ_corr = spZ / f_pol
    PZ_Z = spPZ / spZ
    ## print("spZ / f_pol Inf #",len(spZ_corr[np.isinf(spZ_corr)]))
    ## print("spZ / f_pol Nan #",len(spZ_corr[np.isnan(spZ_corr)]))
    ## print("spZ / f_pol < 1e-5 #",len(spZ_corr[spZ_corr < 1e-5]))
    ## print("---")
    ## print("spPZ / spZ Inf #",len(PZ_Z[np.isinf(PZ_Z)]))
    ## print("spPZ / spZ Nan #",len(PZ_Z[np.isnan(PZ_Z)]))
    ## print("spPZ / spZ < 1e-5 #",len(PZ_Z[PZ_Z < 1e-5]))

    spZ[spZ < 1e-5] = 0.0 # we'll fix this in the plot routine ... map to minX
    spZ_corr[spZ_corr < 1e-5] = 0.0
    
    # Set plot limits, log space
    minY = -4.0; maxY = 0.1
    minX = -5.1; maxX = 0.1
    histMax = 10
    print("***** Z_star vs PZ/Z")
    genDensityPlot(spZ_corr, # x-axis
                   PZ_Z,  # y-axis
                   spMass, spPF, z,
                   'Z_fpol-vs-Z_P-MassHistLogFullNorm',
                   r'{\rm log}\, \left(Z_{\star}/Z_{\odot}\right)',
                   r'{\rm log}\, \left(Z_{{\rm P},\star} / Z_{\star} \right)')

    if z == '05.00':
        print('***** <Z> vs PZ/Z')
        genDensityPlot(spZ, # x-axis
                    PZ_Z, # y-axis
                    spMass, spPF, z,
                    'Z-vs-Z_P-MassHistLogFullNorm',
                    r'{\rm log}\, \left(\langle Z_{\star}\rangle/Z_{\odot}\right)',
                    r'{\rm log}\, \left(\langle Z_{{\rm P},\star}\rangle /\langle Z_{\star}\rangle \right)',
                    normByPMass=False)
    gc.collect()
    ## if z == '05.00':
    ##     print('Making comparison plots ...')
    ##     minX = -8.1
    ##     histMax = 10
    ##     print('***** Z_star vs PZ/Z')
    ##     genDensityPlot(spZ_corr, # x-axis
    ##                 PZ_Z,  # y-axis
    ##                 spMass, spPF, z,
    ##                 'Z_fpol-vs-Z_P-MassHistLogFullNorm-fullX',
    ##                 r'{\rm log}\, \left(Z_{\star}/Z_{\odot}\right)',
    ##                 r'{\rm log}\, \left(\langle Z_{{\rm P},\star}\rangle /\langle Z_{\star}\rangle \right)')
        
    ##     print('***** <Z> vs PZ/Z')
    ##     genDensityPlot(spZ, # x-axis
    ##                PZ_Z, # y-axis
    ##                spMass, spPF, z,
    ##                'Z-vs-Z_P-MassHistLogFullNorm-fullX',
    ##                r'{\rm log}\, \left(\langle Z_{\star}\rangle/Z_{\odot}\right)',
    ##                r'{\rm log}\, \left(\langle Z_{{\rm P},\star}\rangle /\langle Z_{\star}\rangle \right)', normByPMass=False)
