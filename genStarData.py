#
# Compute sfr for a set of RAMSES output files, store
# in a np array and output to a text file.
# 09 Sep 2015
# Rick Sarmento
#
# Purpose:
#  Reads in RAMSES output files and computes the various
#  pynbody as well as star-location plots.
#
# Method:
#  Sample 3 sps from the s.s collection and plot Z, PGF,
#  Temp, Rho and locations for each... 
#
# Revision history
#

# ##########################################################
# Output the star particle data for all stars in the snap
# ##########################################################
def outputStarData(z,s):
    print ("Saving data...")
    np.savetxt("spLoc_%05.2lf.txt"%z,s.s['pos'],header="x\ty\tz",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spMass_%05.2lf.txt"%z,s.s['mass'],header="mass",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spZ_%05.2lf.txt"%z,s.s['zsolar'],header="Z_solar",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spPPF_%05.2lf.txt"%z,s.s['ppf'],header="ppf",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spPZ_%05.2lf.txt"%z,s.s['pzsolar'],header="pz_solar",comments='') # comments='' gets rid of "#" in output
    np.savetxt("spBT_%05.2lf.txt"%z,s.s['tform'].in_units('yr'),header="birth_time_yr",comments='') # comments='' gets rid of "#" in output
    print ("Star Particle Data Files saved for %05.1lf"%z)
    return




# ##########################################################
# Generate and save stellar mass totals for z
# ##########################################################
def outputMassTots(file,s,data,indx):
    # Make a local copy of the stars, set pzf to 0.0
    # for small values.
    temp = s.s
    temp['pzf'][temp['pzf'] < 1e-10] = 0.0
    #s.s['pzf'][s.s['pzf'] <= 1e-10] = 0.0

    # Get total stellar mass
    totalStarMass = np.sum(s.s['mass'].in_units('Msol'))
    print ("Total stellar mass is %.2e Msol @ z~%.2lf" % (totalStarMass,z))
    
    # Total POP3 star mass
    totalPop3StarMass = np.sum(s.s['mass'].in_units('Msol') * s.s['ppf'])
    print ("Total POPIII stellar mass is %.2e Msol @ z~%.2lf" % (totalPop3StarMass,z))

    # Total of only totally pristine star particles ***** Entire SP must be < z_crit!! ****** 
    # Note that we have to use mass-fraction "metal" here... 
    totalPop3SubcritZStarMass = np.sum(s.s['mass'][s.s['metal'] < z_crit])
    print ("Total sp mass for spZ < Z_crit is %.2e Msol @ z~%.2lf" % (totalPop3SubcritZStarMass,z))
    
    # Total polluted stars
    totalPolStarMass = np.sum(s.s['mass'].in_units('Msol') * (1.0-s.s['ppf']))
    print ("Total polluted stellar mass is %.2e Msol @ z~%.2lf" % (totalPolStarMass,z))

    # Compute the mass of stars' primordial metals
    # 1 - ppf is the polluted fraction of stars, by mass
    # pzf/Z is then the fraction of primordial metals
    # (1-ppf) * pzf / Z is hence the fraction of stars polluted only by primordial metals
    totalPriStarMass = np.sum(s.s['mass'].in_units('Msol') * (1.0 - s.s['ppf']) * temp['pzf']/ s.s['metal'])
    print ("Total primordial stellar mass is %.2e Msol @ z~%.2lf" % (totalPriStarMass,z))

    totalGasMass  = np.sum(s.g['mass'].in_units('Msol'))
    print ("Total gas mass is %.2e Msol @ z~%.2lf" % (totalGasMass,z))

    totalPristGasMass = np.sum(s.g['mass'].in_units('Msol') * s.g['pgf'])
    print ("Total pristing gas mass is %.2e Msol @ z~%.2lf" % (totalPristGasMass,z))

    # Compute total star mass for stars NOT made of primordial Z
    totNonPrimordStarMass = np.sum(s.s['mass'].in_units('Msol') * (1.0 - s.s['ppf']) * (1.0 - temp['pzf']/s.s['metal']))
    print ("Total non-primordial stellar mass is %.2e Msol @ z~%.2lf" % (totNonPrimordStarMass, z))

    data[indx] = [z, min(s.s['tform'].in_units('yr')), max(s.s['tform'].in_units('yr')),
                  totalStarMass,totalPop3StarMass,totalPolStarMass,totalPriStarMass,
                  totalGasMass,totalPristGasMass,totalPop3SubcritZStarMass,totNonPrimordStarMass]
    print ("Writing out totals for z=%.2lf"%z)
    np.savetxt("%s_new-data3Mpc.txt"%file,data[0:indx],comments='',header=headerStr) # comments='' gets rid of "#" in output
    del temp
    return data

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
#import matplotlib as mpl
#mpl.use('Agg') # Needed if you don't have an X-server
#import matplotlib.pyplot as plt
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

#mpl.rcParams['figure.figsize'] = (12,8)
#mpl.rcParams['font.size'] = 18
pynbody.ramses.multiprocess_num = 24
pynbody.config['number_of_threads'] = 96

# #######################################################
# Initialize variables
# #######################################################
# Create list of files to process
files = [f for f in os.listdir('.') if re.match(r'output_000[0-9][0-9]', f)]
files =[
    "output_00004",
    "output_00005",  # z = 20 
    "output_00006",
    "output_00008", 
    "output_00010", 
    "output_00012", 
    "output_00016", 
    "output_00021",  # z = 14
    "output_00029",  # z = 13
    "output_00034", 
    "output_00041", 
    "output_00048",  # z = 11.5 
    "output_00053", 
    "output_00062", 
    "output_00070", 
    "output_00080", 
    "output_00090",  # z = 9 
    "output_00101", 
    "output_00113"   # z = 8
    ]

files.sort()
print ("Files to process %d"%len(files))

np.set_printoptions(precision=3,suppress=True)

key = np.empty([len(files),2]) # a map of z to boxsize, 

data      = np.empty([len(files),11])
headerStr = 'z\ttStart\ttEnd\ttotStarMass\ttotPop3StarMass\ttotPollStarMass\ttotPrimordStarMass\ttotGasMass\ttotPristGasMass\ttotSubcritStarMass\ttotNonPrimordStarMass'

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

    z_crit = 2.0e-7 # Mass fraction
    print ("Normalize Z data...")
    s.g['zsolar'] = s.g['metal'] * 50.0         # Solar units
    s.g['pzsolar'] = s.g['pzf'] * 50.0          # Solar units

    print("handle small values in Z, PZ")
    s.s['metal'][s.s['metal']<1e-10]    = 1e-10 # Avoid 0's for now...
#    s.s['ppf'][s.s['ppf']>(1.0-1e-6)]  = 1.0   # Set limits in post-processing, not here!
    s.s['pzf'][s.s['pzf']<1e-10]        = 1e-10 # Avoid 0's for now...
    s.s['zsolar'] = s.s['metal'] * 50.0         # Solar units
    s.s['pzsolar'] = s.s['pzf'] * 50.0          # Solar units

    # Fix birth times
    print (s.derivable_keys()[1:5])
    pynbody.analysis.ramses_util.get_tform(s,"/home/earnric/bin/part2birth");

    # ##########################################################
    # Output SP data
    outputStarData(z,s)
    # Save a list of redshift and the boxsize for the files we're processing
    key[indx] = ["%.2lf"%z,"%.2lf"%s.properties['boxsize'].in_units('kpc')]

    # ##########################################################
    # Output mass totals for Mathematica overall histogram
    data = outputMassTots(file,s,data,indx)

    del s
    gc.collect()
    gc.collect()

# Sort the keys... Since we'll process the files in Z order
# from low to high, we need to reverse the entries
keys2=key[np.argsort(key[:, 0])]
print ("Keys unsorted: \n",key)
print ("Keys sorted: \n"  ,keys2)
# ##########################################################
np.savetxt("zKeysForSPfiles2.txt",keys2,header="z,boxsize_kpc",comments='')
# ##########################################################
# Write out the array to a file
print ("Mass totals\n",data)
np.savetxt("StarParicle-data3Mpc.txt",data,comments='',header=headerStr) # comments='' gets rid of "#" in output
del data
print ("Done")
