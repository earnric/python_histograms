import gc
from matplotlib.pylab import *
import matplotlib.pyplot as plt 
import pynbody
import pynbody.plot.sph as sph
import sys

pynbody.ramses.multiprocess_num = 8
pynbody.config['number_of_threads'] = 24

rcParams['figure.figsize'] = (12,9)
rcParams['font.size'] = 36


## s = pynbody.load('output_00016')  # z=16
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
s.g['pzsolar'][s.g['pzsolar']< 1e-5] = 1e-6

# ##############################
#  BOX SIZE
# ##############################
sbox = 120.0 / (1.0 + z) * 0.71 # 40 kpc comoving box
smallbox = str(sbox) + " kpc"
print(smallbox)

# z=16, i=0 : -20.26, -122.12, 40.22
# z=16, i=770 : 43.96, 26.56, 121.94
# z=16, i=1540 : 44.54, 19.03, 118.77

## rx,ry,rz = 44.54, 19.03, 118.77 # 00016
## tic = 0.5
## i=1540 # 00016, 1540

## The one that matches the scatter plot example... 
## rx,ry,rz =  43.96, 26.56, 121.94 # 00016
## tic = 0.5
## i=770 # 00016, 1540

# #################################################
# #################################################

# z=8, i=0 : -121.77, -104.01, -203.10
# z=8, i=265793 : -0.88, 215.56, 143.27
# z=8, i=531586 : 77.55, 54.08, 218.96
# z=8, i=797379 : 137.45, -200.82, -220.27

## rx,ry,rz = -121.77, -104.01, -203.10
## tic = 0.8
## i=0  # 00121, 0

## The one at z=8 that matches the scatter plot
rx,ry,rz = -0.88, 215.56, 143.27
tic = 2.5
i=265793

print(rx,ry,rz)

impData = s[pynbody.filt.Cuboid(str((rx-sbox/2.0)) + " kpc", str((ry-sbox/2.0)) + " kpc",str((rz-sbox/4.0)) + " kpc",
                                str((rx+sbox/2.0)) + " kpc", str((ry+sbox/2.0)) + " kpc",str((rz+sbox/4.0)) + " kpc")]

minPGF = impData.g['pgf'].min()
maxPGF = impData.g['pgf'].max()
if minPGF == maxPGF:
    print("MIN = MAX PGF ... aborting")
    sys.exit(0)

rect = [0.15,0.15,0.85,0.9]
coords= [-rx,-ry,-rz] # Translation requires negative of the coord
with pynbody.transformation.translate(impData,coords):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.yaxis.set_visible(False)
    #ax.set_xlabel(fontsize=40)
    fileOut = "img_log_ax_Z-z=%.1lf-%i.pdf"% (z,i)
    titleStr = r"$Z_{\odot}$ - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
    print (titleStr)
    sph.image(impData.g,qty="zsolar",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False, subplot=ax,
                      log=True, vmax=1.0, vmin=1e-7, qtytitle=r"${\rm log}\, \langle Z\rangle/Z_{\odot}$",
                      approximate_fast=False
                      ); #vmin=0.006, vmax=1.0,
    ax.set_xticks([-tic, 0, tic])
    plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight')
    plt.close(fig)
    del(ax)
    gc.collect()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.yaxis.set_visible(False)
    fileOut = "img_log_ax_PZ-z=%.1lf-%i.pdf"% (z,i)
    titleStr = r"$Z_{P, \odot}$ - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
    print (titleStr)
    sph.image(impData.g,qty="pzsolar",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False,subplot=ax,
                      log=True, vmax=1.0, vmin=1e-7, qtytitle=r"${\rm log}\,\langle Z_{P}\rangle /Z_{\odot}$",
                      approximate_fast=False );
    ax.set_xticks([-tic, 0, tic])
    plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight')
    plt.close(fig)
    del(ax)

    ## fig = plt.figure()
    ## ax = fig.add_subplot(111)
    ## ax.xaxis.set_visible(True)
    ## ax.yaxis.set_visible(True)
    ## fileOut = "img_log_ax_v_r-z=%.1lf-%i.pdf"% (z,i)
    ## titleStr = "Radial Vel - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
    ## print (titleStr)
    ## sph.image(impData.g,qty="vr",width=smallbox,cmap="RdYlBu_r", denoise=True ,av_z=False, units="km s^-1",
    ##                   log=True, approximate_fast=False, #vmin=10**-4.5, vmax=10**2.0,
    ##                   subplot=ax,qtytitle=r"${\rm log}\, v_{\rm r}\, {\rm km/s}$"); 
    ## ax.set_yticks([-tic, 0,tic])
    ## plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight')
    ## plt.close(fig)
    ## del(ax)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    fileOut = "img_log_ax_Density-z=%.1lf-%i.pdf"% (z,i)
    titleStr = "Density - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
    print (titleStr)
    sph.image(impData.g,qty="rho",width=smallbox,cmap="terrain", denoise=True ,av_z=False, units="m_p cm^-3",
                      log=True, approximate_fast=False, vmin=10**-4.5, vmax=10**2.0,
                      subplot=ax,qtytitle=r"${\rm log}\,m_{\rm p}/{\rm cm}^{3}$"); 
    ax.set_yticks([-tic, 0,tic])
    plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight')
    plt.close(fig)
    del(ax)
    gc.collect()

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    fileOut = "img_log_ax_cs-z=%.1lf-%i.pdf"% (z,i)
    titleStr = "cs - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
    print (titleStr)
    sph.image(impData.g,qty="cs",width=smallbox,cmap="RdYlBu_r", denoise=True ,av_z=False, units="km s^-1",
                      log=True, approximate_fast=False,subplot=ax,  vmin=10**0.5, vmax=1e3,
                      qtytitle=r"${\rm log}\, c_{\rm s}\, {\rm km/s}$"); 
    plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight') 
    plt.close(fig)
    del(ax)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    fileOut = "img_log_ax_vt-z=%.1lf-%i.pdf"% (z,i)
    titleStr = "$v_{t}$ @ z = %.1f" % z + "\nThin slice @ %s" % str(coords)
    print(titleStr)
    sph.velocity_image(impData.g, qty="tv", width=smallbox, cmap = "RdYlBu_r", units="km s^-1",
                   mode='quiver', quiverkey=False,qtytitle=r"${\rm log}\, v_{t}\, {\rm km/s}$",
                   density = 1.0, vector_resolution=40, vmin=10**0.5, vmax=1e3, subplot=ax,
                   show_cbar=True, vector_color='black') # gist_ncar
    plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight') 
    plt.close(fig)
    del(ax)
    gc.collect()
                    
    try: 
        fig = plt.figure()
        ax = fig.add_subplot(111)
        fileOut = "img_log_ax_PGF-z=%.1lf-%i.pdf"% (z,i)
        titleStr = "PGF - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
        print (titleStr)
        sph.image(impData.g,qty="pgf",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False, subplot=ax,
                log=True, vmax = 1.0, vmin=1e-7, qtytitle=r"${\rm log}\, P}$",approximate_fast=False);
        ax.set_xticks([-tic, 0, tic])
        ax.set_yticks([-tic, 0, tic])
        plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight')
        plt.close(fig)
        del(ax)
    except:
        print("Unable to make PGF plot")
        pass
    
s.g['pzsolar'][s.g['pzsolar']< 1e-5] = 0.0   # Since we are looking at ratio

print("Fraction of PM ", np.sum(impData.g['rho'] * impData.g['x']**3 * impData.g['pzsolar']/impData.g['zsolar'])/np.sum(impData.g['rho'] * impData.g['x']**3))
## print("Largest V_t: ",impData['tv'].max)
