#
# Compute sfr for a set of RAMSES output files, store
# in a np array and output to a text file.
# 12 Nov 2015
# Rick Sarmento
#
# Purpose:
#  Reads in RAMSES output files and outputs masses for sp.
#
# Method:
#


# ##########################################################
# Output the data needed to generate the Mathematica histograms
# ##########################################################
def outputHistogramData(z,i,halo):
    bins    = np.logspace(-10, 0, 51) # Log bins for histogram data 
    psm     = np.zeros(len(bins)-1) # Pristine Stellar mass in that bin
    tsm     = np.zeros(len(bins)-1) # total mass
    primsm  = np.zeros(len(bins)-1) # Primordial stellar mass
#    primsm2  = np.zeros(len(bins)-1) # Primordial stellar mass, corrected 
#    polsm    = np.zeros(len(bins)-1) # Polluted stellar mass

    # By using zsolar, we are binning in units of solar Z.
    print("Number of bin edges %i"%len(bins))
    for indx2,j in enumerate(bins):
        if indx2 < len(bins)-1 :
            cond = (halo.s['zsolar'] >= j) & (halo.s['zsolar'] < bins[indx2+1]) # This selects for sp's in the bin
            psm[indx2] = np.sum(halo.s['ppf'][cond] * halo.s['mass'][cond]) # Pop III
            #polsm[indx2] = (1.0-np.sum(halo.s['ppf'][cond]) * halo.s['mass'][cond]) # Polluted
            tsm[indx2] = np.sum(halo.s['mass'][cond]) # Total star particle mass in this bin

            # For sp's that are in our bin (above):
            #   Compute the polluted fraction (1-ppf)
            #   Compute the fraction of pristine metals: pzf/Z
            #   Compute the mass of polluted stars that are polluted only by pristine metals
            # Note that 'metal' maybe 1e-10, but that means pzf = 0.0 and we get the correct
            # mass
            primsm[indx2] = np.sum((1.0-halo.s['ppf'][cond]) * 
                                   (halo.s['pzf'][cond] / halo.s['metal'][cond]) * halo.s['mass'][cond])

            print("Bin left edge %.2e"%j)
            print("Bin total mass %.0lf"%tsm[indx2])
            print("Bin pop III mass %.0lf"%psm[indx2])
            print("Bin primord mass %.0lf"%primsm[indx2])
            print("pzf > metal:")

            if psm[indx2] > tsm[indx2] :
                print("Somehow POP III particle mass > total star mass", psm[indx2], tsm[indx2])
            if primsm[indx2] > tsm[indx2] :
                print("Somehow primordial mass > total star mass",primsm[indx2],tsm[indx2])
                print("Number of pzf > Z = %i"%len(halo.s['pzf'][halo.s['metal'] < halo.s['pzf']]))
                print("Points pzf, Z")
                print(halo.s['pzf'][halo.s['metal'] < halo.s['pzf']])
                print(halo.s['metal'][halo.s['metal'] < halo.s['pzf']])
            print("\n")

            # allPZ = cond & (halo.s['pzf'] > halo.s['metal'])
            # halo.s['pzf'][allPZ] = halo.s['metal'][allPZ]
            # primsm2[indx2] = np.sum((1.0-halo.s['ppf'][cond]) * (halo.s['pzf'][cond] / halo.s['metal'][cond]) * halo.s['mass'][cond])

    pristMass   = np.sum(psm)
    totMass     = np.sum(tsm)
    primordMass = np.sum(primsm)
        
    # Save the total masses
#    np.savetxt("polSM_z%.1lf-%i.txt"%(z,i),polsm,comments='') # comments='' gets rid of "#" in output
#    np.savetxt("primord2SM_z%.1lf-%i.txt"%(z,i),primsm2,comments='') # comments='' gets rid of "#" in output
    np.savetxt("primordSM_z%.1lf-%i.txt"%(z,i),primsm,comments='') # comments='' gets rid of "#" in output
    np.savetxt("pristSM_z%.1lf-%i.txt"%(z,i),psm,comments='') # comments='' gets rid of "#" in output
    np.savetxt("totSM_z%.1lf-%i.txt"%(z,i),tsm,comments='') # comments='' gets rid of "#" in output
    print "Saving data..."
    
    # Update this first one with the center coords so these are all relative to the center pt.
    np.savetxt("spLocHaloOnly_%05.2lf_%i.txt"%(z,i),halo.s['pos'],header="x\ty\tz",comments='') 
    np.savetxt("spMassHaloOnly_%05.2lf_%i.txt"%(z,i),halo.s['mass'],header="mass",comments='') 
    np.savetxt("spZHaloOnly_%05.2lf_%i.txt"%(z,i),halo.s['zsolar'],header="Z_solar",comments='') 
    np.savetxt("spPPFHaloOnly_%05.2lf_%i.txt"%(z,i),halo.s['ppf'],header="ppf",comments='') 
    np.savetxt("spPZHaloOnly_%05.2lf_%i.txt"%(z,i),halo.s['pzsolar'],header="pz_solar",comments='') 
    np.savetxt("spBTHaloOnly_%05.2lf_%i.txt"%(z,i),halo.s['tform'].in_units('yr'),header="birth_time_yr",comments='') 
    print "Star Particle Data Files saved for %05.1lf, %i"%(z,i)
    return

# ##########################################################
# ##########################################################
##
## Main program
##
# ##########################################################
# ##########################################################
#%matplotlib inline
import os
import re # Regular expression matching
import matplotlib as mpl
mpl.use('Agg') # Needed if you don't have an X-server
import matplotlib.pyplot as plt
import numpy as np
import pynbody
import pynbody.plot as pp
import pynbody.plot.sph as sph
import pynbody.filt as filt
import pynbody.units as units
import pynbody.analysis.profile as profile
import sys, os, glob, pickle, gc
import math
#import mmap

mpl.rcParams['figure.figsize'] = (12,12)
mpl.rcParams['font.size'] = 22
pynbody.ramses.multiprocess_num = 32
pynbody.config['number_of_threads'] = 128

# #######################################################
# Initialize variables
# #######################################################
# Create list of files to process
files = [f for f in os.listdir('.') if re.match(r'output_000[0-9][0-9]', f)]
files =[
    ## "output_00007",
    ## "output_00008",
    ## "output_00011",
    "output_00016",
    ## "output_00020",
    "output_00026",
    ## "output_00033",
    "output_00037",
    ## "output_00043",
    "output_00050",
    ## "output_00058",
    "output_00066",
    ## "output_00068",
    "output_00073",
    "output_00084",
    ## "output_00097",
    ## "output_00107",
    "output_00121",
    ## "output_00136",
    "output_00152",
    ## "output_00169",
    "output_00191",
    ## "output_00214",
    "output_00241"
    ]
files.sort()
print( "Files to process %d"%len(files))

np.set_printoptions(precision=3,suppress=True)

key = np.empty([len(files),2]) # a map of z to boxsize, 

# Loop and load data
# If we wanted to just loop over filename - "for file in files:"
for indx, file in enumerate(files):
    gc.collect()
    print ("Processing file %s" % file)
    s = pynbody.load("./"+file)
    s['pos'] -= 0.5
    s.physical_units();
    
    z = 1/s.properties['a']-1
    print ("Redshift = %.2lf" % z)
    boxsizestring = "%.2lf kpc" % s.properties['boxsize'].in_units('kpc')
    print ("Boxsize @ this redshift %s"%boxsizestring)
    print ("aexp %.3lf" % s.properties['a'])
    bs = float(s.properties['boxsize'])

    z_crit = 1.0e-5 # solar units
    print ("Normalize Z data...")

    s.s['ppf'][s.s['ppf'] > (1.0-1e-6)]  = 1.0  # Don't keep tiny polluted fractions... 
    s.s['zsolar']  = s.s['metal'] * 50.0        # Solar units
    s.s['pzsolar'] = s.s['pzf'] * 50.0          # Solar units

    condition = (s.s['pzf'] > s.s['metal'])
    if ( len(s.s['pzf'][condition]) ):
        print("WHAT!!! There are star particles with pzf > metal (next pzf, then metal)")
        print(s.s['pzf'][condition])
        print(s.s['metal'][condition])
    # Fix birth times
    print (s.derivable_keys()[1:5])
    pynbody.analysis.ramses_util.get_tform(s,"/home1/02744/earnric/bin/part2birth")
    
    sbox = 5.0 / (1.0 + z) * 0.71 # 5 kpc comoving box
    print ("5 kpc comoving Small box @ z=%.2lf is %.2lf"%(z,sbox))

    step = len(s.s)/3 # Note integer division
    if step == 0: step = 1
    for i in range(0,len(s.s),step):
        x,y,zz = - s.s['pos'][i] # get a star, this is an offset (-1 * pos)
        print ("Star offset [%.2lf %.2lf %.2lf]"%(x,y,zz))
        print ("Boxsize @ this z ", bs)
        # Ensure we center such that we can depict a 40 kpc box
        if (abs(x) > bs/2.0 - sbox/2.0):
            x = np.sign(x) * (bs/2.0 - sbox/2.0) # closest x coord to sp
        if (abs(y) > bs/2.0 - sbox/2.0):
            y = np.sign(y) * (bs/2.0 - sbox/2.0) # closest y coord to sp
#        if (abs(zz) > bs/2.0 - sbox/2.0): # Commented out so we get the exact z-coord of the sp
#            zz = np.sign(zz) * (bs/2.0 - sbox/2.0) # we could go off the end... ok
                    
        # Create the size of our box for sph image plots
        # May not be centered on sp - we have adjusted x,y such that we
        # don't go 'off the edge' for a 'smallbox' sized image
        coords= [x,y,zz] # NOTE WE ARE TRANSLATING... coords is -1 * location of sp of interest
        if sbox > 1.0: smallbox= str(sbox) + ' kpc'
        else: smallbox= str(sbox*1000.0) + ' pc'
        print ("Plotsize: ",smallbox)

        # ##########################################################
        # Save the index and coordinates of the center of our image...
        np.savetxt("z%05.2lf_SpCoord_%i.txt"%(z,i), -1.0*np.array(coords),fmt='%.2lf',comments='') 
        # ##########################################################

        rx,ry,rz = -1 * np.array([x,y,zz]) # For filters, we need the actual coords, not translation coords.
        print ("Cuboid center [%.2lf %.2lf %.2lf]"%(rx,ry,rz))
        # Note that the thickness in z is only 1/2 that of the other dimensions
        halo = s[pynbody.filt.Cuboid(str((rx-sbox/2.0)) + " kpc", str((ry-sbox/2.0)) + " kpc",str((rz-sbox/4.0)) + " kpc",
                                        str((rx+sbox/2.0)) + " kpc", str((ry+sbox/2.0)) + " kpc",str((rz+sbox/4.0)) + " kpc")]
        
        # Compute stellar masses and save for the halo histogram plot
        outputHistogramData(z,i,halo)

    del s
    gc.collect()
    gc.collect()

print ("Done")
