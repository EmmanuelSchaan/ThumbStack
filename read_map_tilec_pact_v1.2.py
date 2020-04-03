import numpy as np, time
import matplotlib.pyplot as plt

import flat_map
reload(flat_map)
from flat_map import *

import cmb
reload(cmb)
from cmb import *

#from enlib import enmap, utils, powspec
from pixell import enmap, utils, powspec, enplot, reproject

import healpy as hp

# copy rotfuncs.py somewhere on your python path,
# so you can import it
import rotfuncs


# for cori
#plt.switch_backend('Agg')


#########################################################################
#########################################################################

mapId = sys.argv[1]


if mapId=='cmbksz_d56':
   pathFig = "./figures/cmb_map/tilec_pact_cmbksz_d56_v1.2.0/"
   pathOut = "./output/cmb_map/tilec_pact_cmbksz_d56_v1.2.0/"
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_deep56/tilec_single_tile_deep56_cmb_map_v1.2.0_joint.fits"
   pathTilecMask = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_deep56/tilec_mask.fits" # no hit count: use the mask
if mapId=='cmbksz_boss':
   pathFig = "./figures/cmb_map/tilec_pact_cmbksz_boss_v1.2.0/"
   pathOut = "./output/cmb_map/tilec_pact_cmbksz_boss_v1.2.0/"
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_boss/tilec_single_tile_boss_cmb_map_v1.2.0_joint.fits"
   pathTilecMask = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_boss/tilec_mask.fits"





#########################################################################
#########################################################################
# Create folders for output and figures if needed
if not os.path.exists(pathFig):
    os.makedirs(pathFig)
if not os.path.exists(pathOut):
    os.makedirs(pathOut)


#########################################################################
# CMB power spectra

cmb1_6 = StageIVCMB(beam=1.6, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)
cmb7_3 = StageIVCMB(beam=7.3, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)

#########################################################################
# load maps

# mask sits in the older folder
pathPSMask = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/source_mask_s16_simonecat_sn5_cross_20171105.fits"

print "Read maps"
baseMap = enmap.read_map(pathMap)
tilecMask = enmap.read_map(pathTilecMask)
psMask = enmap.read_map(pathPSMask)

print("Map properties:")
print("nTQU, nY, nX = "+str(baseMap.shape))
print("WCS attributes: "+str(baseMap.wcs))


#########################################################################
#########################################################################
# convert map to healpix

print("Convert from CAR to healpix")
hMap = enmap.to_healpix(baseMap)
hTilecMask = enmap.to_healpix(tilecMask)
hPSMask = enmap.to_healpix(psMask)

# re-threshold the PS mask
hPSMask = (hPSMask>0.95).astype(np.float)

# Save healpix maps to files
hp.write_map(pathOut+"healpix_f150_prelim_map_mono.fits", hMap, overwrite=True)
hp.write_map(pathOut+"healpix_f150_prelim_div_mono.fits", hTilecMask, overwrite=True)
hp.write_map(pathOut+"healpix_pointsourcemask.fits", hPSMask, overwrite=True)

nSide = hp.get_nside(hMap)
print("Nside = "+str(nSide))


#########################################################################
# For debugging: downgrade the map for speed
'''
hMap = hp.ud_grade(hMap, 1024)
hTilecMask = hp.ud_grade(hTilecMask, 1024)
hPSMap = hp.ud_grade(hPSMap, 1024)
nSide = hp.get_nside(hMap)
'''

#########################################################################
# plot maps
'''
# plot PS mask
fig=plt.figure(0)
hp.mollview(hPSMask, fig=0, title="Point source mask", max=np.max(hPSMask), coord=None, cbar=False, unit='')
fig.savefig(pathFig+"psmask.pdf")
fig.clf()

# Plot tilec mask
fig=plt.figure(0)
hp.mollview(hTilecMask, fig=0, title="Tilec mask", max=0.5*np.max(hTilecMask), coord=None, cbar=True, unit='')
fig.savefig(pathFig+"tilecmask.pdf")
fig.clf()

# plot unmasked temperature map
mean = 0.   #np.mean(hMap) 
sigma = 110.   #np.std(hMap)
#
fig=plt.figure(0)
hp.mollview(hMap, fig=0, min=mean-3.*sigma, max=mean+3.*sigma, title="T", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"T_unmasked.pdf")
fig.clf()
'''


#########################################################################
# Generate a footprint mask from the hit map

# copy the masked temperature map
footMask = hTilecMask.copy()

# generate a mask from all the null pixels,
footMask = (footMask>0.99).astype(np.float)

fig=plt.figure(0)
hp.mollview(footMask, fig=0, title="Footprint", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"footmask.pdf")
fig.clf()

# save the footprint to file
hp.write_map(pathOut+"hp_footprint_mask.fits", footMask, overwrite=True)

# read the footprint mask
footMask = hp.read_map(pathOut+"hp_footprint_mask.fits")


#########################################################################
# Mask the Milky Way with a Planck mask
'''
# Read the Planck Galactic mask
pathPlanckMask = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/HFI_Mask_GalPlane-apo0_2048_R2.00.fits"
# field, fsky:
#0, 0.20
#1, 0.40
#2, 0.60
#3, 0.70
#4, 0.80
#5, 0.90
#6, 0.97
#7, 0.99
planckMask = hp.read_map(pathPlanckMask, field=3)

# plot the Planck mask in Galactic coordinates
fig=plt.figure(0)
hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"planck_mask_G.pdf")
fig.clf()

# rotate it to equatorial coordinates
rot = hp.Rotator(coord=['G','C'])  # Transforms galactic to equatorial coordinates
print "rotate the map"
planckMask = rot.rotate_map_alms(planckMask)
#planckMask = rot.rotate_map_alms(planckMask, use_pixel_weights=False)
#planckMask = rot.rotate_map_pixel(planckMask)
print "rotation done"

# upgrade to Nside=4096
planckMask = hp.ud_grade(planckMask, nSide)

fig=plt.figure(0)
hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"planck_mask_C.pdf")
fig.clf()

# rethreshold the mask, to make it zeros and ones
planckMask = (planckMask>0.999).astype(np.float)

fig=plt.figure(0)
hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"planck_mask_C_thresh.pdf")
fig.clf()

# save the Planck mask, rotated
hp.write_map(pathOut+"Planck_mask_4096_equatorial_coord.fits", planckMask, overwrite=True)

# read the Planck mask, rotated
planckMask = hp.read_map(pathOut+"Planck_mask_4096_equatorial_coord.fits", planckMask, overwrite=True)
'''

# read the Planck mask
pathPlanckMask = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/Planck_mask_4096_equatorial_coord.fits"
planckMask = hp.read_map(pathPlanckMask)


#########################################################################
# Combine the masks into one, and extend it around the edges

fullMask = footMask * planckMask

print "fsky =", np.sum(fullMask) / len(fullMask)


# Gaussian smooth the mask
sigma = 4. * np.pi/180. # sigma in rad. We will truncate at around 2 sigmas
fullMask = hp.smoothing(fullMask, sigma=sigma)

# threshold the mask
fullMask = (fullMask>0.95).astype(np.float)

# combine it with the point source mask
fullMask *= hPSMask

# save the full mask
hp.write_map(pathOut+"mask_foot_planck_ps.fits", fullMask, overwrite=True)

# read the full mask
fullMask = hp.read_map(pathOut+"mask_foot_planck_ps.fits")

fSky = np.sum(fullMask) / len(fullMask)
print "unmasked fsky =", fSky


# convert to CAR
carFullMask = reproject.enmap_from_healpix(fullMask, tilecMask.shape, tilecMask.wcs, rot=None)
# re-threshold
carFullMask = (carFullMask>0.95).astype(np.float)
# save
enmap.write_map(pathOut+"mask_foot_planck_ps_car.fits", carFullMask)
# check that fsky is the same
checkFSky = np.sum(carFullMask) / (carFullMask.shape[1] * carFullMask.shape[2])
print "checking unmasked fSky =", checkFSky


