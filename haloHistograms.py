#
# 11 Nov 2015
# Rick Sarmento
#
# Purpose:
#  Reads star particle data and creates histograms.
#  The data files are for halo's
#
# Revision history
#  11 Nov 2015 - Initial version

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
import numpy as np

prefix = "./"
z = "16.0"; i = "1540"
totSPM = np.loadtxt(prefix + "totSM_z" + z + "-" + i +".txt", skiprows=1)
ppfSPM = np.loadtxt(prefix + "pristSM_z" + z + "-" + i +".txt", skiprows=1)
ppzSPM = np.loadtxt(prefix + "primordSM_z" + z + "-" + i +".txt", skiprows=1)

# Remember, the data is already binned in the files... range is
##{-10., -9.8, -9.6, -9.4, -9.2, -9., -8.8, -8.6, -8.4, -8.2, -8., \
## -7.8, -7.6, -7.4, -7.2, -7., -6.8, -6.6, -6.4, -6.2, -6., -5.8, -5.6, \
## -5.4, -5.2, -5., -4.8, -4.6, -4.4, -4.2, -4., -3.8, -3.6, -3.4, -3.2, \
## -3., -2.8, -2.6, -2.4, -2.2, -2., -1.8, -1.6, -1.4, -1.2, -1., -0.8, \
## -0.6, -0.4, -0.2, 0.}
# Slightly offset bins so we can see the data.
xrange1 = np.logspace(-10,0,51) # 1e-10, 1e-9.8, 1e-9.6, ... 
xrange2 = np.logspace(-9.96,0.04,51) # 1e-9.96, 1e-9.76, 1e-9.56, ... 
xrange3 = np.logspace(-9.92,0.08,51) # 1e-9.92, 1e-9.72, 1e-9.52, ...

# We're gonna cheat a bit and just move the histogram value to the lower edge of the bin... 
# Our bins are small so this shouldn't even show... 
histRect = [0.1, 0.1, .85, 0.85]
axHist   = plt.axes(histRect)

axHist.set_xscale('log')
axHist.set_yscale('log')
axHist.set_ylim([100,1e6])

axHist.set_xlabel("log $Z_{\odot}$", size=24)
axHist.set_ylabel('log $M_{\odot}$', size=24)

axHist.plot(xrange1[:50],totSPM,linestyle='-',  linewidth = 2.0, marker = 'd')
axHist.plot(xrange2[:50]+0,ppfSPM,linestyle='--', linewidth = 2.0, marker = 's' )
axHist.plot(xrange3[:50],ppzSPM,linestyle='-.', linewidth = 2.0, marker = '^' )
axHist.legend(['Total','Pop III','Primoridal Z'],loc='upper left')

plt.show()
