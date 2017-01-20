import sys, gc
sys.path.append('/Users/earnric/Google Drive/ASU/Codes/PythonCode/modules')

import matplotlib as mpl
import matplotlib.pyplot as plt
from itertools import chain
import numpy as np
import math as ma

from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
import astropy 

import loadSP as lsp
import cosmo
import halos

plt.rcParams['figure.figsize'] = (13,11)
plt.rcParams['font.size'] = 32
plt.rcParams.update({'font.size': 22, 'font.family': 'STIXGeneral', 'mathtext.fontset': 'stix'})

#
# Define function for string formatting of scientific notation
#
def exp_tex(float_number):
    """
    Returns a string representation of the scientific
    notation of the given number formatted for use with
    LaTeX or Mathtext.
    """
    neg = False
    if float_number == 0.0:
        return r"$0.00$"
    elif float_number < 0.0:
        neg = True

    exponent = np.floor(np.log10(abs(float_number)))
    mantissa = float_number/10**exponent
    if neg:
        mantissa = -mantissa
    mantissa_format = "{:.2f}".format(mantissa)
    return r"${0}\times10^{{{1}}}$".format(mantissa_format, str(int(exponent)))


#
# Compute a virial radius using star particles and the critical density
# of baryons at z
#
def getRv(z, locs, mass, spAges, Z, PZ, PPF, haloPos, distance=15, radiusScaleFac=1.25):
    """ Returns the virial radius of the star forming region based on the mass in stars
    We use a scale factor of 1.25 since we do not account for gas, only star particles.
    Add parameter distances are expected to be consistent.
    """
    # Kinda kludgy: use a fixed distance to grab sps, then use their mass to compute a r_v
    haloStars,halomass,haloages,haloZ,haloPZ,haloPPF = halos.getHaloStars(locs,mass,spAges,Z,PZ,PPF,
                                                                          haloPos, distance)
    totMass = np.sum(halomass) * u.solMass
    radius = ((totMass * 3.0/(4.0 * ma.pi) * 1.0/(myCosmo.critical_density(z) * myCosmo.Ob(z) * 100))**(1.0/3.0)).to(u.kpc) * (1.0 + z)
    radius = radius.value * radiusScaleFac
    return radius


#
# Plot the halos and star particles
#
def pltHalos(z, locs, halosPos, haloMassTots, filename="starsAndHalos_", labelHalos=3, alpha=0.05):
    fig, ax = plt.subplots(1, 1, figsize=(12,12))

    labels = ['halo{0}'.format(i+1) for i in range(len(halosPos))]
    
    # Legend
    meanMassStr = exp_tex(haloMassTots.mean())
    ax.scatter(6,6.2,s=120,facecolors='none', edgecolors='r')
    ax.annotate("Mean mass @ z={} is {}".format(z,meanMassStr) + r" $M_{\odot}$",
                    xy = (5.65, 6), textcoords = 'data', ha = 'right', va = 'bottom',fontsize=18)

    # Plot Data
    ax.scatter(locs[:,0], locs[:,1],s=15, c='b', edgecolors='none', alpha=alpha) 
    ax.scatter(halosPos[:,0], halosPos[:,1], s=120 * haloMassTots/haloMassTots.mean(),
                facecolors='none', edgecolors='r')

    # Plot settings
    ax.grid(b=True, which='major', color='k', linestyle='--',alpha=0.5)
    ax.set_xlabel('Mpc/h')
    ax.set_ylabel('Mpc/h')
    ax.set_xlim([-6.5,6.5])
    ax.set_ylim([-6.5,6.5])

    if labelHalos > 0:
        sign = 1
        for label, x, y in zip(labels[:labelHalos], halosPos[:,0][:labelHalos], halosPos[:,1][:labelHalos]):
            # for label, x, y in zip(chain(labels[:10],labels[-3:]), 
            #                        chain(halosPos[:,0][:10],halosPos[:,0][-3:]),
            #                        chain(halosPos[:,1][:10],halosPos[:,1][-3:])):
            ax.annotate(
                label, 
                xy = (x, y), xytext = (-20 * sign, 20),
                textcoords = 'offset points', ha = 'right', va = 'bottom',fontsize=18,
                bbox = dict(boxstyle = 'round,pad=0.05', fc = 'yellow', alpha = 0.75),
                arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
            sign *= -1
    plt.savefig(prefix+filename+"{}.pdf".format(z),dpi=fig.dpi,bbox_inches='tight')
    del ax
    gc.collect()
    return

#
# Plot a sample of galaxies
# Note that the routine expects units of Mpc!
#
def pltGalaxies(z, locs, mass, spAges, Z, PZ, PPF, halosPos, theRange=[0,1,2], distance=15, filename="SP-galaxy_z_"):
    print("arcsec/kpc = {:.3f} @ z={}".format(arcPerKpcP,z))
    for halo in theRange:
        fig, ax = plt.subplots(1, 1, figsize=(6,6))
        
        print("********************** Halo {} **********************".format(halo))
        radius = getRv(z, locs * 1000.0, mass, spAges, Z, PZ, PPF, halosPos[halo] * 1000.0, distance)
        haloStars,halomass,haloages,haloZ,haloPZ,haloPPF = halos.getHaloStars(locs * 1000.0,mass,spAges,Z,PZ,PPF,
                                                                            halosPos[halo] * 1000.0, radius)
        totMass = np.sum(halomass)
        
        print("Stellar mass in r_v halo {} is {:.2e}, radius {:.2f}".format(halo, np.sum(halomass * u.solMass),radius))

        # Plot Data
        ax.scatter(haloStars[:,0],haloStars[:,1],s=20 * halomass/halomass.mean(),facecolors='b', edgecolors='None') 
        ax.annotate('halo{0}'.format(halo+1), xy = (0, 0), xytext = (-80, 80),
            textcoords = 'offset points', ha = 'right', va = 'bottom',fontsize=14,
            bbox = dict(boxstyle = 'round,pad=0.05', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))

        # Legend
        massStr = exp_tex(totMass)
        ax.annotate('Stellar mass {}'.format(massStr) + r" $M_{\odot}$", xy = (0, 0), xytext = (.5, 0.01),
                    textcoords = 'axes fraction', ha = 'left', va = 'bottom',fontsize=14)

        # Plot settings
        ax.grid(b=True, which='major', color='b', linestyle='--', alpha=0.5)
        ax.set_xlabel('x [kpc]')
        ax.set_ylabel('y [kpc]')
        ax.set_xlim([-radius,radius])
        ax.set_ylim([-radius,radius])
        hax2 = ax.twinx()
        hax2.set_ylim([0,radius/(1+z)/h0 * 2.0 * arcPerKpcP])
        hax2.set_ylabel(r'$arcsec$')
        hax2.yaxis.labelpad = -3
        plt.savefig(prefix+filename + "{}_{}.pdf".format(z,halo),dpi=fig.dpi,bbox_inches='tight')
        del ax
        gc.collect()
    return


#
# Compute the center of mass of a group of masses
#
def centMass(relLocs,masses):
    totMass = np.sum(masses)
    x = np.sum(masses * relLocs[:,0])/totMass
    y = np.sum(masses * relLocs[:,1])/totMass
    z = np.sum(masses * relLocs[:,2])/totMass
    return (x,y,z)


# Setup cosmology for distance computations.
# Cosmology from my paper
myCosmo = FlatLambdaCDM(H0=71.0, Om0=0.267, Ob0=0.0449, name='myCosmo')
comovingBox = 12.0 # Mpc/h
Odm = 0.222 # Om - Ob = Odm
Ob  = 0.0449
h0  = 0.71
rho0dm = myCosmo.critical_density0 * Odm 
rho0b  = myCosmo.critical_density0 * Ob   
dmPartMass = 5.589e5 * u.Msun # From 12 Mpc sim... 1024 init conditions... 

#
# In case we need it ... 
#
def totMass(numParts):
    return dmPartMass * numParts

################################################
################################################
## MAIN
################################################
################################################

#
# Map redshift to output file directory
#
files = {17.0: 16, 16.0: 21, 15.0: 28, 14.0: 35, 13.0: 44, 12.0: 54, 11.0: 65, 10.5: 72, 10.0: 78, 9.5: 86, 9.0: 93, 8.5: 102}

if len(sys.argv) < 2:
    print("Please supply redshift")
    sys.exit()
z   = float(sys.argv[1])

if z in files:
    num = files[z]
else:
    print("No files for z={}".format(z))
    sys.exit()
    
prefix = '/Users/earnric/Research/RAMSES-Data/12Mpc-n.48-eps20/'
infoSize = lsp.getBoxSize(num,prefix) ### REMEMBER to Load from info_xxxx ###
print("Using info size, astro convert = {:.3f}".format((infoSize * u.cm).to(u.kpc)) )
boxsizekpc = (infoSize * u.cm).to(u.kpc).value

rho_critDMz = myCosmo.critical_density(z) * myCosmo.Odm(z) # rho_crit @ z for DM only
print("rho_crit,dm,{:.1f} = {:.3e}".format(z,rho_critDMz))
rho_critBz = myCosmo.critical_density(z) * myCosmo.Ob(z) # rho_crit @ z for baryons only
print('rho_crit,b,{:.1f} = {:.3e}'.format(z,rho_critBz))
gc.collect()


##############################################
# LOAD the star particle data
##############################################
locsRaw, mass, bt, Z, PZ, PPF = lsp.loadSPgz(z)
##############################################
# Need to correct the Z and PZ for the pristine
# fraction of the SP... 
##############################################
fpol = 1.0-PPF
fpol[fpol < 1e-5] = 0.0 # Don't form polluted stars when fpol is very small

Z_corr = Z/fpol
Z_corr[np.isinf(Z_corr)] = 0.9e-5 # Make the Z of the particles < Z_crit (1e-5)
Z_corr[Z_corr < 1e-5]    = 0.9e-5 # Place them all just below Z_crit
gc.collect()
locsRaw.size

# Convert to comoving Mpc ... 
locs = (locsRaw/boxsizekpc) * comovingBox

# Compute age of the star particles...
# spAges is in Myr ... 

baseAge = cosmo.ageAtz(71,z) # ageAtZ return Myr, parameters are H and z
print("My code - Age at z={} is {:.2f} Myr".format(z,baseAge))
spAges = baseAge - bt/1e6    # bt from SP file is in yr, convert to Myr... 
print("Star particle age range {:.2f} Myr - {:.2f} Myr".format(min(spAges),max(spAges)))
print("Star particle Z range   {:.3e} Solar - {:.4f} Solar".format(min(Z),max(Z)))

arcPerKpcC = myCosmo.arcsec_per_kpc_comoving(z=z)
print('AstroPy - Comoving @ {}, {:.4f}'.format(z,arcPerKpcC))
arcPerKpcP = (myCosmo.arcsec_per_kpc_proper(z=z)).value
print('AstroPy -  Proper @  {}, {:.4f}'.format(z,arcPerKpcP))
#print('AstroPy -  Proper @  {}, {:.4f}'.format(z,1/arcPerKpcP))

kpcPerArcSec = cosmo.arcAtz(71,z=z,om=0.267) # size of an arcsec at the specified reshift in kpc
print('My code - Comoving @ {}, {:.4f} arcsec / kpc'.format(z,1/kpcPerArcSec/(1+z)))
print('My code -  Proper @  {}, {:.4f} arcsec / kpc'.format(z,1/kpcPerArcSec))
#print('My code -  Proper @  {}, {:.4f} kpc / arcsec'.format(z,kpcPerArcSec))


print("Age of the universe: {:.3f}".format(myCosmo.age(z)))

##
## Load the halo locations
## hop halo pos file key:    #   npart,mass,cont.frac,xc,yc,zc,uc,vc,wc
## Scale locations to co-moving coordates

##
## THESE HALOS ARE BASED ON STAR PARTICLE LOCATIONS ##
##
prefix = '/Users/earnric/Research/RAMSES-Data/12Mpc-n.48-eps20/12Mpc-n.48-eps20-output/'
halosRawPos = lsp.loadHaloGrps(num,prefix=prefix+'hopStarData/') # Returns positions normalized to (-0.5, 0.5)
halossizes = lsp.loadHaloSizes(num,prefix=prefix+'hopStarData/') # Returns positions normalized to (-0.5, 0.5)

#halosPos = halosRawPos * boxsizekpc
halosPos = halosRawPos * comovingBox  # Comoving box

print("Num of halos @ {} = {}".format(z,len(halosRawPos)))


#
# Compute the masses and virial radii of our halos
#
haloMassTots = []   # Halo masses
halorv       = []   # Halo r_v using baryons
for h in halosPos * 1000.0:
    radius = getRv(z, locs * 1000.0, mass, spAges, Z, PZ, PPF, h, 15)

    haloStars,halomass,haloages,haloZ,haloPZ,haloPPF = halos.getHaloStars(locs * 1000.0, mass, spAges, Z, PZ, PPF,
                                                                          h, radius)
    haloMassTots.append(np.sum(halomass)) # Append to array of halo masses
    halorv.append(np.sum(radius))         # Append to array of halo r_v
    if np.sum(halomass) == 0.0:
        print('ZERO mass halo',h)

#
# Clean up the halos -- get rid of empties
#
haloMassTots = np.array(haloMassTots) # Masses of our halos
halorv       = np.array(halorv)       # Virial radius of our halos
emptyHaloCond    = (haloMassTots == 0.0)
nonEmptyHaloCond = (haloMassTots > 0.0)
print("Number of empty halos {}".format(len(haloMassTots[emptyHaloCond])))
# haloMassTots & halorv are lists with same len as halosPos.. Must modify all at once
#haloMassTots = haloMassTots[nonEmptyHaloCond]
#halorv       = halorv[nonEmptyHaloCond]
#halosPos     = halosPos[nonEmptyHaloCond]
# etc.

# Plot the HOP halos...    
pltHalos(z, locs, halosPos, haloMassTots, filename="starsAndHalos2_")

# Plot some galaxies...
# The routine expects units of Mpc!
numHalos = len(halosPos)
#pltGalaxies(z, locs, mass, spAges, Z, PZ, PPF, halosPos,
#                chain(range(3),range(numHalos-3,numHalos)), distance=15, filename="SP-galaxy2_z_")
pltGalaxies(z, locs, mass, spAges, Z, PZ, PPF, halosPos, range(3), distance=15, filename="SP-galaxy2_z_")

####################################################
# Find all the 'unfound' halos...
# Systematically remove all of the stars in a
# complete list of stars that are within v_r of a halo.
####################################################
#
# Find all the star particles that are not
# within a halo's r_v
#
stars   = locs * 1000.0 # Copy stars into a list we can manipulate. Work in kpc...
masses  = mass
spAges2 = spAges
Zs      = Z
PZs     = PZ
PPFs    = PPF
haloXMassTots = []   # Halo masses
haloXrv       = []   # Halo r_v using baryons
print("Total # sps:",len(stars))
for aHalo in halosPos * 1000:
    # Get a r_v and a mass in stars
    radius = getRv(z, stars, masses, spAges2, Zs, PZs, PPFs, aHalo, 15)
    haloStars,halomass,haloages,haloZ,haloPZ,haloPPF = halos.getHaloStars(stars, masses, spAges2,
                                                                              Zs, PZs, PPFs, aHalo, radius)
    totalHaloMass = np.sum(halomass) * u.solMass

    centrd    = stars - aHalo 
    dists     = np.linalg.norm(centrd,axis=1)
    # Keep stars (and assoc'd info) for stars NOT in this halo
    stars     = stars[dists > radius]
    masses    = masses[dists > radius]
    spAges2   = spAges2[dists > radius]
    Zs        = Zs[dists > radius] 
    PZs       = PZs[dists > radius] 
    PPFs      = PPFs[dists > radius]
    
# At the end of the loop, any stars that remain are not within the r_v of a halo... 
print("Remaining orphaned star particles:",len(stars))
orphans = stars


####################################################
# At this point _stars_ are the set of star
# particles that are not within r_v kpc of the center of a halo
####################################################
# We can now group these stars into new halos. 
#  
# We'll assume they are within a 15 kpc radius, comoving
#
## Use the orphaned stars from above!
rad = 15 # kpc
i = 0
minMass       = 1e3 # Solar masses
newHalos      = []
newHaloMasses = []
newHaloRv     = []
while len(stars):
    aStar   = stars[0]      # Grab a star... 
    centrd  = stars - aStar # Compute distances between aStar and all others... 
    dists   = np.linalg.norm(centrd,axis=1) # Linear distances between stars and aStar
    nearby  = (dists <= rad) # Nearby criteria
    
    haloNStars = centrd[nearby] # Get nearby stars locations normalized to 'aStar'...
    haloStars  = stars[nearby]  # Get nearby stars locations
    haloMasses = masses[nearby] # Get masses
    haloSAges  = spAges2[nearby]# Get ages
    haloZs     = Zs[nearby]     # Get Z
    haloPZs    = PZs[nearby]    # Get masses
    haloPPFs   = PPFs[nearby]   # Get masses
    
    haloMass = np.sum(haloMasses)
    if  haloMass >= minMass:
        # print('Found halo with mass {:.3e} and {} sps'.format(haloMass * u.solMass, len(haloStars)))
        localCMass = centMass(haloStars,haloMasses)
        newHalos.append(localCMass)    # Halo's location
        newHaloMasses.append(haloMass) # Halo's mass
        radius = getRv(z, stars,masses,spAges2,Zs,PZs,PPFs, aStar, 15)
        # print("Halo rv = {}".format(radius))
        newHaloRv.append(radius)       # Halo's r_v

    stars   = stars[dists > radius]  # Remove the stars from the list
    masses  = masses[dists > radius] # Get masses
    spAges2 = spAges2[dists > radius] # Get Z
    Zs      = Zs[dists > radius] # Get Z
    PZs     = PZs[dists > radius] # Get masses
    PPFs    = PPFs[dists > radius] # Get masses
#     print("Halo {:2d} has {:3d} sps and mass {:.3e}".format(i,len(haloStars),np.sum(haloMasses) * u.solMass))
    i += 1
    gc.collect()
print("New halos count: {}".format(len(newHalos)))
print("Stars left: {}".format(len(stars)))
newHalos = np.array(newHalos)
newHaloMasses = np.array(newHaloMasses)

#
# Plot the new halos.... 
#
pltHalos(z, orphans/1000.0, newHalos/1000.0, newHaloMasses, filename="extraHalos2_",labelHalos=0,alpha=0.5)

print("Saving extra halos... {}".format(len(newHalos)))

## Noramlize locations to (0,1.0) and save to be compat with HOP format
#HOP header-> "   #   npart       mass  cont.frac         xc         yc         zc         uc         vc         wc"
header = "   x       y       z"
newHalos = newHalos/1000.0/comovingBox # Convert to range -.5, .5
#newHalos += 0.5 # Convert to range 0, 1
# newHalos now in inter
np.savetxt(prefix+"extraHalos2_{:.1f}.txt".format(z),newHalos,fmt="% .3e",header=header, delimiter=' ', comments='')

#
# Combine and save
#
header="   xc         yc         zc         mass          rv"
finalHaloPoss   = np.concatenate((halosRawPos, newHalos))
finalHaloMasses = np.concatenate((haloMassTots, newHaloMasses))
lenMass = len(finalHaloMasses)
finalHaloRv     = np.concatenate((halorv, newHaloRv))
lenRv   = len(finalHaloRv)
final = np.hstack((finalHaloPoss,finalHaloMasses.reshape(lenMass,1),finalHaloRv.reshape(lenRv,1)))

# Save final halo information... 
np.savetxt(prefix+"finalHalos2_{:.1f}.txt".format(z),final,fmt="% .3e",header=header, delimiter=' ', comments='')

# Print out all halos and stars... 
pltHalos(z, locs, (np.array(finalHaloPoss)) * comovingBox, finalHaloMasses, filename="AllHalos_",labelHalos=0,alpha=0.05)

print("Done...")
