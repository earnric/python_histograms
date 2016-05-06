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

# z 16, 8
files =  [ 'output_00016','output_00121' ]
coords = [ (44.54, 19.03, 118.77), (-121.77, -104.01, -203.10) ]
tics =   [ 0.5, 0.8 ]
theis =  [ 1540, 0 ]

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

## rx,ry,rz = -121.77, -104.01, -203.10
## tic = 0.8
## i=0  # 00121, 0

for file, (rx,ry,rz), tic, i in zip(files,coords,tics,theis):
    print ("Loading {}".format(file))
    s = pynbody.load(file)  # z=16
    s['pos']
    s['pos'] -= 0.5
    s.physical_units();
    
    z = 1/s.properties['a']-1
    print ("Redshift = {:.0f}".format(z))
    boxsizestring = "%.2f" % s.properties['boxsize'].in_units('kpc')
    boxsizestring += " kpc"
    print ("Box size: {}".format(boxsizestring))
    
    gc.collect()
    
    s.g['zsolar']  = s.g['metal'] * 50.0         # Solar units
    s.g['pzsolar'] = s.g['pzf'] * 50.0           # Solar units
    s.g['zsolarOrig']=s.g['zsolar']
    
    s.g['zsolar'][s.g['zsolar']< 1e-5] = 1e-6
    s.g['pzsolar'][s.g['pzsolar']< 1e-5] = 1e-6
    # compute enhanced Z
    s.g['zEnhanced'] = s.g['zsolarOrig']/(1.0-s.g['pgf'])
    s.g['zEnhanced'][s.g['zEnhanced'] == np.nan] = 1e-6
    s.g['zEnhanced'][s.g['zEnhanced'] == np.inf] = 1e-6
    s.g['zEnhanced'][s.g['zEnhanced']< 1e-5] = 1e-6

    print("Number of nan: {:d}".format(len(s.g['zEnhanced'][s.g['zEnhanced'] == np.nan])))
    print("Number of inf: {:d}".format(len(s.g['zEnhanced'][s.g['zEnhanced'] == np.inf])))
    
    sbox = 40.0 / (1.0 + z) * 0.71 # 40 kpc comoving box
    smallbox = str(sbox) + " kpc"
    print(smallbox)
    print(rx,ry,rz)
    
    impData = s[pynbody.filt.Cuboid(str((rx-sbox/2.0)) + " kpc", str((ry-sbox/2.0)) + " kpc",str((rz-sbox/4.0)) + " kpc",
                    str((rx+sbox/2.0)) + " kpc", str((ry+sbox/2.0)) + " kpc",str((rz+sbox/4.0)) + " kpc")]
    
    rect = [0.15,0.15,0.85,0.9]
    coords= [-rx,-ry,-rz] # Translation requires negative of the coord
    with pynbody.transformation.translate(impData,coords):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.yaxis.set_visible(False)
        #ax.set_xlabel(fontsize=40)
        fileOut="img_log_ax_Z-z=%.1lf-%i.pdf"% (z,i)
        titleStr = "$Z_{\odot}$ - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
        print (titleStr)
        sph.image(impData.g,qty="zsolar",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False, subplot=ax,
            log=True, vmax=1.0, vmin=1e-6, qtytitle=r"${\rm log}\, \langle Z\rangle/Z_{\odot}$", approximate_fast=False
            ); #vmin=0.006, vmax=1.0,
        ax.set_xticks([-tic, 0, tic])
        ## ax.set_xticklabels(['-0.7','0','0.7'])
        plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight')
        plt.close(fig)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        #ax.yaxis.set_visible(False)
        #ax.set_xlabel(fontsize=40)
        fileOut="img_log_ax_Z-f_pol-z=%.1lf-%i.pdf"% (z,i)
        titleStr = r"$Z_{\odot}/f_{pol}$ - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
        print (titleStr)
        sph.image(impData.g,qty="zEnhanced",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False, subplot=ax,
            log=True, vmax=1.0, vmin=1e-6, qtytitle=r"${\rm log}\, Z/Z_{\odot}$", approximate_fast=False
            ); #vmin=0.006, vmax=1.0,
        ax.set_xticks([-tic, 0, tic])
        ax.set_yticks([-tic, 0, tic])
        plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight')
        plt.close(fig)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        fileOut="img_log_ax_PGF-z=%.1lf-%i.pdf"% (z,i)
        titleStr = "PGF - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
        print (titleStr)
        sph.image(impData.g,qty="pgf",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False,subplot=ax,
            log=True, vmax = 1.0, vmin=1e-6, qtytitle=r"${\rm log}\, P}$",approximate_fast=False
            ); #vmin=0.006, vmax=1.0,
        ax.set_xticks([-tic, 0, tic])
        ax.set_yticks([-tic, 0, tic])
        plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight')
        plt.close(fig)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.yaxis.set_visible(False)
        fileOut="img_log_ax_PZ-z=%.1lf-%i.pdf"% (z,i)
        titleStr = r"$Z_{P, \odot}$ - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
        print (titleStr)
        sph.image(impData.g,qty="pzsolar",width=smallbox,cmap="nipy_spectral", denoise=True ,av_z=False,subplot=ax,
            log=True, vmax=1.0, vmin=1e-6, qtytitle=r"${\rm log}\,\langle Z_{P}\rangle /Z_{\odot}$",
            approximate_fast=False);

        ax.set_xticks([-tic, 0, tic])
        ## ax.set_xticklabels(['-0.7','0','0.7'])
        plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight')
        plt.close(fig)
    
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.xaxis.set_visible(False)
        fileOut="img_log_ax_Density-z=%.1lf-%i.pdf"% (z,i)
        titleStr = "Density - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
        print (titleStr)
        sph.image(impData.g,qty="rho",width=smallbox,cmap="terrain", denoise=True ,av_z=False, units="m_p cm^-3",
            log=True, approximate_fast=False,subplot=ax,qtytitle=r"${\rm log}\,m_{\rm p}/{\rm cm}^{3}$"
            ); #vmin=0.006, vmax=1.0,
        ax.set_yticks([-tic, 0,tic])
        ## ax.set_yticklabels(['-0.7','0','0.7'])
        plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight')
        plt.close(fig)
    
        del(ax)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        fileOut="img_log_ax_Temp-z=%.1lf-%i.pdf"% (z,i)
        titleStr = "Temp - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
        print (titleStr)
        sph.image(impData.g,qty="temp",width=smallbox,cmap="RdYlBu_r", denoise=True ,av_z=False,
            log=True, approximate_fast=False,subplot=ax,qtytitle=r"${\rm log}\,Temp\, {\rm K}$"
            ); #vmin=0.006, vmax=1.0,
        plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight') 
        plt.close(fig)
        
        del(ax)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        fileOut="img_log_ax_cs-z=%.1lf-%i.pdf"% (z,i)
        titleStr = "cs - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
        print (titleStr)
        sph.image(impData.g,qty="cs",width=smallbox,cmap="RdYlBu_r", denoise=True ,av_z=False, units="km s^-1",
            log=True, approximate_fast=False,subplot=ax,qtytitle=r"${\rm log}\, c_{\rm s}\, {\rm km/s}$"
            ); #vmin=0.006, vmax=1.0,
        plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight') 
        plt.close(fig)
        
        ## del(ax)
        ## fig = plt.figure()
        ## ax = fig.add_subplot(111)
        ## ax.xaxis.set_visible(False)
        ## ax.yaxis.set_visible(False)
        ## fileOut="img_log_ax_v-z=%.1lf-%i.pdf"% (z,i)
        ## titleStr = "v - z = %.1lf" % z# + "\n[%.2lf %.2lf %.2lf]"%(rx,ry,rz)
        ## print (titleStr)
        ## sph.image(impData.g,qty="v",width=smallbox,cmap="RdYlBu_r", denoise=True ,av_z=False, units="km s^-1",
        ##                   log=True, approximate_fast=False,subplot=ax,qtytitle=r"${\rm log}\, v\, {\rm km/s}$"
        ##                   ); #vmin=0.006, vmax=1.0,
        ## plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight') 
        ## plt.close(fig)

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        fileOut="img_log_ax_vt-z=%.1lf-%i.pdf"% (z,i)
        titleStr = "$v_{t}$ @ z = %.1f" % z + "\nThin slice @ %s" % str(coords)
        print(titleStr)
        sph.velocity_image(impData.g, qty="tv", width=smallbox, cmap = "gist_ncar", units="km s^-1",
                       mode='quiver', quiverkey=False,qtytitle=r"${\rm log}\, v_{t}\, {\rm km/s}$",
                       density = 1.0, vector_resolution=40, vmax=1e3, subplot=ax,
                       show_cbar=True, vector_color='black')
        plt.savefig(fileOut,dpi=fig.dpi,bbox_inches='tight') 
        plt.close(fig)

    gc.collect()
    del(s)
    gc.collect()
