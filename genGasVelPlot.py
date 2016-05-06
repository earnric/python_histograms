import gc
from matplotlib.pylab import *
import matplotlib.pyplot as plt 
import pynbody
import pynbody.plot.sph as sph
import mmap

pynbody.ramses.multiprocess_num = 12
pynbody.config['number_of_threads'] = 24

rcParams['figure.figsize'] = (12,9)
rcParams['font.size'] = 36


##s = pynbody.load('output_00016')  # z=16
s = pynbody.load('output_00121')  # z=8
print ("Loading {}".format('file...'))
s['pos']
s['pos'] -= 0.5
s.physical_units();

z = 1/s.properties['a']-1
print ("Redshift = {:.0f}".format(z))
boxsizestring = "%.2f" % s.properties['boxsize'].in_units('kpc')
boxsizestring += " kpc"
print (boxsizestring)

gc.collect()

s.g['zsolar']  = s.g['metal'] * 50.0         # Solar units
s.g['pzsolar'] = s.g['pzf'] * 50.0           # Solar units

s.g['zsolar'][s.g['zsolar']< 1e-5] = 1e-6
s.g['pzsolar'][s.g['pzsolar']< 1e-5] = 0.0   # Since we are looking at ratio


sbox = 40.0 / (1.0 + z) * 0.71 # 40 kpc comoving box
smallbox = str(sbox) + " kpc"
print(smallbox)

# z=16, i=0 : -20.26, -122.12, 40.22
# z=16, i=770 : 43.96, 26.56, 121.94
# z=16, i=1540 : 44.54, 19.03, 118.77

## rx,ry,rz = 44.54, 19.03, 118.77 # 00016
## tic = 0.5
## i=1540 # 00016, 1540

# z=8, i=0 : -121.77, -104.01, -203.10
# z=8, i=265793 : -0.88, 215.56, 143.27
# z=8, i=531586 : 77.55, 54.08, 218.96
# z=8, i=797379 : 137.45, -200.82, -220.27

rx,ry,rz = -121.77, -104.01, -203.10
tic = 0.8
i=0  # 00121, 0

print(rx,ry,rz)

impData = s[pynbody.filt.Cuboid(str((rx-sbox/2.0)) + " kpc", str((ry-sbox/2.0)) + " kpc",str((rz-sbox/4.0)) + " kpc",
                                str((rx+sbox/2.0)) + " kpc", str((ry+sbox/2.0)) + " kpc",str((rz+sbox/4.0)) + " kpc")]

rect = [0.15,0.15,0.85,0.9]
coords= [-rx,-ry,-rz] # Translation requires negative of the coord
with pynbody.transformation.translate(impData,coords):
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(True)
    ax.yaxis.set_visible(True)
    fileOut="img_log_ax_vel-z=%.1lf-%i.pdf"% (z,i)
    titleStr = "v - z = %.1lf" % z
    print (titleStr)
    sph.image(impData.g,qty="vel",width=smallbox,cmap="RdYlBu_r", denoise=True ,av_z=False, # units="km s^-1",
                      log=True, approximate_fast=False, subplot=ax, qtytitle=r"${\rm log}\, v\, {\rm km/s}$"
                      ); #vmin=0.006, vmax=1.0,
    plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight') 
    plt.close(fig)
    

print("Fraction of PM ",
      np.sum(impData.g['rho'] * impData.g['x']**3 * impData.g['pzsolar']/impData.g['zsolar'])/np.sum(impData.g['rho'] * impData.g['x']**3))
