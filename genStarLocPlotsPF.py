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
# Create list of files to process
## locFiles  = [f for f in os.listdir('.') if re.match(r'spLoc_[0-9]*.*.txt', f)]
## massFiles = [f for f in os.listdir('.') if re.match(r'spMass_[0-9]*.*.txt', f)]
## ZFiles    = [f for f in os.listdir('.') if re.match(r'spZ_[0-9]*.*.txt', f)]
## ppfFiles  = [f for f in os.listdir('.') if re.match(r'spPPF_[0-9]*.*.txt', f)]
## ppzFiles  = [f for f in os.listdir('.') if re.match(r'spPZ_[0-9]*.*.txt', f)]
locFiles  = ['spLoc_16.00.txt','spLoc_08.00.txt' ]
massFiles = ['spMass_16.00.txt','spMass_08.00.txt']
ZFiles    = ['spZ_16.00.txt','spZ_08.00.txt']
ppfFiles  = ['spPPF_16.00.txt','spPPF_08.00.txt']
ppzFiles  = ['spPZ_16.00.txt','spPZ_08.00.txt']

# Sort the files
# Not needed if manually create list.
#locFiles.sort();massFiles.sort();ZFiles.sort();ppfFiles.sort();ppzFiles.sort()

#zs  = np.loadtxt("zKeysForSPfiles.txt",skiprows=1)
zs  = np.loadtxt("PaperzKeysForSPfiles.txt",skiprows=1)

if len(locFiles) != len(zs):
    print("Diff # of files to process than I have redshifts")
    os._exit(1)

dotNorm = 10.0  # For dot-size scaling
comovbox = 5.0 
i = 0
for locF,massF,ZF,ppfF,pzfF in zip(locFiles,massFiles,ZFiles,ppfFiles,ppzFiles):
    z    = zs[i][0]
    bs   = zs[i][1] 
    sbox = comovbox / (1.0 + z) * 0.71 # Create a box that's sbox kpc physical
    print ("z=%.3lf"%z)
    spPosExp = r'z%05.2lf_SpCoord_[0-9]*.txt' %z # These are the x,y,z locations of interest from the genStarPlots6.py
    poisFiles  = [f for f in os.listdir('.') if re.match(spPosExp, f)]
    poisFiles.sort()

    for indx, poi in enumerate(poisFiles):
        print( "Star locs: ", locF)
        print( "Coord file: ", poi)
        locs = np.loadtxt(locF,skiprows=1) 
        mass = np.loadtxt(massF,skiprows=1)
        #Z    = np.loadtxt(ZF,skiprows=1)
        ppf  = np.loadtxt(ppfF,skiprows=1)
        #pzf  = np.loadtxt(pzfF,skiprows=1)


        x,y,zz = np.loadtxt(poi)  # star particle location of interest
        xo,yo,zo = x,y,zz
        print( "Star offset [%.2lf %.2lf %.2lf]"%(x,y,zz))
        print( "Boxsize at this z (physical) %.4lf"%bs)
        print( "Our box at this z (physical) %.4lf"%sbox)
        # Ensure we center such that we can depict the sbox kpc box
        if (abs(x) > bs/2.0 - sbox/2.0):
            x = np.sign(x) * (bs/2.0 - sbox/2.0)
        if (abs(y) > bs/2.0 - sbox/2.0):
            y = np.sign(y) * (bs/2.0 - sbox/2.0)
        if (abs(zz) > bs/2.0 - sbox/2.0):
            zz = np.sign(zz) * (bs/2.0 - sbox/2.0)
        print( "Star offset adjusted for boxsize [%.2lf %.2lf %.2lf]"%(x,y,zz))
            
        xmin = x - sbox/2.0
        xmax = x + sbox/2.0
        ymin = y - sbox/2.0
        ymax = y + sbox/2.0
        zmin = zz - sbox/2.0 
        zmax = zz + sbox/2.0
        print( "Box range x(%.2lf,%.2lf), y(%.2lf,%.2lf), z(%.2lf,%.2lf)"%(xmin,xmax,ymin,ymax,zmin,zmax))
        xcond = ((locs[:,0] >= xmin) & (locs[:,0] <= xmax))
        ycond = ((locs[:,1] >= ymin) & (locs[:,1] <= ymax))
        zcond = ((locs[:,2] >= zmin) & (locs[:,2] <= zmax))

        locsFilt = locs[xcond & ycond & zcond] # filter: just get stars in plot region
        if len(locsFilt) == 0: continue
        massFilt = mass[xcond & ycond & zcond] # filter: just get stars in plot region
        #Zfilt    = Z[xcond & ycond & zcond] # filter: just get stars in plot region
        ppfFilt  = ppf[xcond & ycond & zcond] # filter: just get stars in plot region
        #pzfFilt  = pzf[xcond & ycond & zcond] # filter: just get stars in plot region

        # Fix zeros
        ppfFilt[ppfFilt == 0.0] = 1e-10
        #pzfFilt[pzfFilt == 0.0] = 1e-10

        # Normalize locations
        locsFilt = locsFilt - np.array([x,y,zz])

        rng1 = (ppfFilt < 1.e-5)
        rng2 = ((ppfFilt >= 1.e-5) & (ppfFilt < 1.e-3))
        rng3 = ((ppfFilt >= 1.e-3) & (ppfFilt < 1.e-1))
        rng4 = (ppfFilt >= 1.e-1)
        p1=np.log10(ppfFilt[rng1])
        p2=np.log10(ppfFilt[rng2])
        p3=np.log10(ppfFilt[rng3])
        p4=np.log10(ppfFilt[rng4])

        # Create the starlocation Z plot
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
        # Fine-tune figure; hide x ticks for top plots and y ticks for right plots
        plt.setp([a.get_xticklabels() for a in [ax1,ax2]], visible=False)
        plt.setp([a.get_yticklabels() for a in [ax2,ax4]], visible=False)
        plt.setp([a.set_xlim([-sbox/2.0,sbox/2.0]) for a in [ax1,ax2,ax3,ax4]])
        plt.setp([a.set_ylim([-sbox/2.0,sbox/2.0]) for a in [ax2,ax4,ax3,ax4]])

        xcoord = locsFilt[:,0]; ycoord = locsFilt[:,1]
        ax1.scatter(xcoord[rng1], ycoord[rng1], s=massFilt[rng1]/dotNorm, c=p1, cmap=cmap,vmin=-7.5, vmax=0)
        ax2.scatter(xcoord[rng2], ycoord[rng2], s=massFilt[rng2]/dotNorm, c=p2, cmap=cmap,vmin=-7.5, vmax=0)
        ax3.scatter(xcoord[rng3], ycoord[rng3], s=massFilt[rng3]/dotNorm, c=p3, cmap=cmap,vmin=-7.5, vmax=0)
        ax4.scatter(xcoord[rng4], ycoord[rng4], s=massFilt[rng4]/dotNorm, c=p4, cmap=cmap,vmin=-7.5, vmax=0)
        # create a second axes for the colorbar
        ax5 = fig.add_axes([0.85, 0.1, 0.025, 0.85])
        cb = mpl.colorbar.ColorbarBase(ax5, cmap=cmap, norm=norm, spacing='proportional',
                            ticks=ticks, boundaries=bounds, format='%1i')
        
        ## fig.suptitle('Star Particles z=%.1lf, %.1lf kpc comoving\nCenter: [%.2lf,%.2lf,%.2lf]'%
        ##              (z,comovbox,xo,yo,zo),size=24)
        xpos = ax1.get_xlim()[0] - 0.07 * ax1.get_xlim()[0]
        ypos = ax1.get_ylim()[0] - 0.08 * ax1.get_ylim()[0]
        bbox = {'facecolor':'white', 'alpha':0.75, 'pad':3}
        ax1.text(xpos,ypos,'$PF <\, 10^{-5}$',
                 bbox=bbox, fontsize=20)
        ax2.text(xpos,ypos,'$10^{-5} \leq\, PF\, <\, 10^{-3}$',
                 bbox=bbox,fontsize=20)
        ax3.text(xpos,ypos,'$10^{-3} \leq\, PF\, <\, 10^{-1}$',
                 bbox=bbox,fontsize=20)
        ax4.text(xpos,ypos,'$10^{-1} \leq\, PF$',
                 bbox=bbox,fontsize=20)

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
    
        ax5.set_ylabel('$log\; PF$', size=34)
        plt.subplots_adjust(left=0.15, bottom=0.1, right=0.84, top=0.95, wspace=.05, hspace=.05)
        plt.savefig("SP_PF_locs4Panel_z=%.1lf-%i.pdf"%(z,indx))
#        plt.show()
        print( "Files saved: SP_PF_locs4Panel_z=%.1lf-%i.pdf"%(z,indx))
        print( "\n")
        plt.close()

    i=i+1
