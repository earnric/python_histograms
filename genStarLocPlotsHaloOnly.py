#
# Compute histogram of metallicity (weighted by sp mass)
# 31 Aug 2015
# Rick Sarmento
#
# Purpose:
#  Read in RAMSES data and compute a Z histogram weighted
#  by the particles' pristine fraction
#
# Method:
#
# Revision history
#

import os
import re # Regular expression matching
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import sys, os, glob, pickle, gc
import math
#import mmap

mpl.rcParams['figure.figsize'] = (12,9)
mpl.rcParams['font.size'] = 26

## ## ## Setup colormap ## ## ##
# define the colormap
cmap = plt.cm.jet
# extract all colors from the .jet map
cmaplist = [cmap(i) for i in range(cmap.N)]
# force the first color entry to be grey
cmaplist[0] = (.5,.5,.5,1.0)
# create the new map
cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)

# define the bins and normalize
bounds = np.linspace(-7.5,0,16)
ticks  =[-7,-6,-5,-4,-3,-2,-1,0]
norm   = mpl.colors.BoundaryNorm(bounds, cmap.N)
## ## ##

# #######################################################
# Initialize variables
# #######################################################

dotNorm = 10.0  # For dot-size scaling
comovbox = 5.0 
z = 5.0
i = 2562828
sbox = comovbox / (1.0 + z) * 0.71 # Create a box that's sbox kpc physical

x,y,zz = np.loadtxt("z05.00_SpCoord_2562828.txt")
locs = np.loadtxt("spLocHaloOnly_05.00_2562828.txt",skiprows=1) - np.array([x,y,zz]) # Normalize
mass = np.loadtxt("spMassHaloOnly_05.00_2562828.txt",skiprows=1)
Z    = np.loadtxt("spZHaloOnly_05.00_2562828.txt",skiprows=1) # SOLAR units in the file... 
#ppf  = np.loadtxt("spPPFHaloOnly_05.00_2562828.txt",skiprows=1)
#pzf  = np.loadtxt("spPZHaloOnly_05.00_2562828.txt",skiprows=1)

print(len(Z[Z<1.0e-5]) # = 1.0e-10
#pzf[pzf<1.0e-5] = 1.0e-10

rng1 = (Z < 1.e-5)
rng2 = ((Z >= 1.e-5) & (Z < 1.e-3))
rng3 = ((Z >= 1.e-3) & (Z < 1.e-1))
rng4 = (Z >= 1.e-1)
print(Z[rng4])
print(locs[rng4])
z1=np.log10(Z[rng1])
z2=np.log10(Z[rng2])
z3=np.log10(Z[rng3])
z4=np.log10(Z[rng4])

# Create the starlocation Z plot
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
plt.setp([a.get_xticklabels() for a in [ax1,ax2]], visible=False)
plt.setp([a.get_yticklabels() for a in [ax2,ax4]], visible=False)
plt.setp([a.set_xlim([-sbox/2.0,sbox/2.0]) for a in [ax1,ax2,ax3,ax4]])
plt.setp([a.set_ylim([-sbox/2.0,sbox/2.0]) for a in [ax2,ax4,ax3,ax4]])

xcoord = locs[:,0]; ycoord = locs[:,1]
ax1.scatter(xcoord[rng1], ycoord[rng1], s=mass[rng1]/dotNorm, c=z1, cmap=cmap,vmin=-7.5, vmax=0)
ax2.scatter(xcoord[rng2], ycoord[rng2], s=mass[rng2]/dotNorm, c=z2, cmap=cmap,vmin=-7.5, vmax=0)
ax3.scatter(xcoord[rng3], ycoord[rng3], s=mass[rng3]/dotNorm, c=z3, cmap=cmap,vmin=-7.5, vmax=0)
ax4.scatter(xcoord[rng4], ycoord[rng4], s=mass[rng4]/dotNorm, c=z4, cmap=cmap,vmin=-7.5, vmax=0)
# create a second axes for the colorbar
ax5 = fig.add_axes([0.85, 0.1, 0.025, 0.85])
cb = mpl.colorbar.ColorbarBase(ax5, cmap=cmap, norm=norm, spacing='proportional',
                                ticks=ticks, boundaries=bounds, format='%1i')
        
## fig.suptitle('Star Particles z=%.1lf, %.1lf kpc comoving\nCenter: [%.2lf,%.2lf,%.2lf]'%
##              (z,comovbox,xo,yo,zo),size=24)
xpos = ax1.get_xlim()[0] - 0.07 * ax1.get_xlim()[0]
ypos = ax1.get_ylim()[0] - 0.08 * ax1.get_ylim()[0]
bbox = {'facecolor':'white', 'alpha':0.75, 'pad':3}
ax1.text(xpos,ypos,'$Z_{\odot} <\, 10^{-5}$',bbox=bbox, fontsize=20)
ax2.text(xpos,ypos,'$10^{-5} \leq\, Z_{\odot}\, <\, 10^{-3}$',bbox=bbox,fontsize=20)
ax3.text(xpos,ypos,'$10^{-3} \leq\, Z_{\odot}\, <\, 10^{-1}$',bbox=bbox,fontsize=20)
ax4.text(xpos,ypos,'$10^{-1} \leq\, Z_{\odot}$',bbox=bbox,fontsize=20)

startx, endx = ax1.get_xlim(); dx = (endx-startx) * 0.1
starty, endy = ax1.get_ylim(); dy = (endy-starty) * 0.1
formatter = FormatStrFormatter('%.2f')
ax1.yaxis.set_ticks([starty+dy, 0,endy-dy]); ax1.yaxis.set_major_formatter(formatter)
ax3.yaxis.set_ticks([starty+dy, 0,endy-dy]); ax3.yaxis.set_major_formatter(formatter)
ax3.xaxis.set_ticks([startx+dx, 0,endx-dx]); ax3.xaxis.set_major_formatter(formatter)
ax4.xaxis.set_ticks([startx+dx, 0,endx-dx]); ax4.xaxis.set_major_formatter(formatter)

ax3.set_xlabel('x kpc')
ax4.set_xlabel('x kpc')
ax1.set_ylabel('y kpc')
ax3.set_ylabel('y kpc')

# Control number of ticks
#ax3.locator_params(nbins=3)
#ax4.locator_params(nbins=3)
#ax1.locator_params(nbins=3)

ax5.set_ylabel('$log\; Z_{\odot}$', size=34)
plt.subplots_adjust(left=0.15, bottom=0.1, right=0.84, top=.95, wspace=.05, hspace=.05)
plt.savefig("SP_Z_locs4Panel_z=%.1lf-%i.pdf"%(z,indx))
    #        plt.show()
print( "Files saved: SP_Z_locs4Panel_z=%.1lf-%i.pdf"%(z,indx))
print( "\n")
plt.close()
np.savetxt("4PanelPlotMassTotals_z=%.1lf-%i.txt"%(z,indx),massFilt,comments='')
np.savetxt("4PanelPlotPopIIIMassTotals_z=%.1lf-%i.txt"%(z,indx),massFilt * ppfFilt,comments='')
np.savetxt("4PanelPlotPZMassTotals_z=%.1lf-%i.txt"%(z,indx), (1.0-ppfFilt * massFilt) * (ppfFilt / Zfilt) * massFilt,comments='')

