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

# path for figures and output
pathFig = "./figures/cmb_map/tilec_pact_cmbksz/"
pathOut = "./output/cmb_map/tilec_pact_cmbksz/"

if not os.path.exists(pathFig):
    os.makedirs(pathFig)
if not os.path.exists(pathOut):
    os.makedirs(pathOut)


#########################################################################
print "Read maps"

# Deep56
pathMapD56 = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_deep56/" + "tilec_single_tile_deep56_cmb_map_v1.0.0_rc_joint.fits"
pathTilecMaskD56 = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_deep56/" + "tilec_mask.fits"   # no hit count: use the mask
d56 = enmap.read_map(pathMapD56)
tilecMaskD56 = enmap.read_map(pathTilecMask56)

# BOSS North
pathMapBN = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_boss/" + "tilec_single_tile_boss_cmb_map_v1.0.0_rc_joint.fits"
pathTilecMaskBN = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_boss/" + "tilec_mask.fits"   # no hit count: use the mask
bn = enmap.read_map(pathMapBN)
tilecMaskBN = enmap.read_map(pathTilecMaskBN)

print("Map properties:")
print("nTQU, nY, nX = "+str(baseMap.shape))
print("WCS attributes: "+str(baseMap.wcs))


#########################################################################
# ACT point source mask

# mask sits in the older folder
pathPSMask = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/source_mask_s16_simonecat_sn5_cross_20171105.fits"
psMask = enmap.read_map(pathPSMask)

#########################################################################
# combine Deep56 and BOSS North maps, on the total ACT footprint

# reproject the D56 and BN maps to the overall ACT footprint
d56 = enmap.project(d56, psMask.shape, psMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
tilecMaskD56 = enmap.project(tilecMaskD56, psMask.shape, psMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
#
bn = enmap.project(bn, psMask.shape, psMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
tilecMaskBN = enmap.project(tilecMaskBN, psMask.shape, psMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)

# combine the masks
tilecMask = tilecMaskD56 + tilecMaskBN
# combine the maps
tilecMap = d56 * tilecMaskD56 + bn * tilecMaskBN


# Plot the combined maps to check


#########################################################################
# Mask the Milky Way with a Planck mask

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
hpPlanckMaskGal = hp.read_map(pathPlanckMask, field=3)

# plot the Planck mask in Galactic coordinates
#fig=plt.figure(0)
#hp.mollview(hpPlanckMaskGal, fig=0, title="Planck Galactic mask (gal coord)", coord=None, cbar=True, unit='')
#fig.savefig(pathFig+"planck_mask_G.pdf")
#fig.clf()

# rotate the Galactic mask to equatorial coord,
# and convert to CAR
planckMask = reproject.enmap_from_healpix_interp(hpPlanckMaskGal, tilecMask.shape, tilecMask.wcs, rot="gal,equ")

print("Planck mask values:")
print("min="+str(planckMask.min())+" max="+str(planckMask.max()))

# Re-threshold the masks
# Galactic mask
planckMask = (planckMask>0.95).astype(np.float)
# footprint
tilecMask = (tilecMask>0.95).astype(np.float)
# point source mask
psMask = (psMask>0.95).astype(np.float)

# combined mask: footprint, point sources, galactic mask
carFullMask = tilecMask * psMask * planckMask


#########################################################################
# save map and mask

# save the masks
enmap.write_map(pathOut+"mask_foot_car.fits", tilecMask)
enmap.write_map(pathOut+"mask_ps_car.fits", psMask)
enmap.write_map(pathOut+"mask_planck_car.fits", planckMask)
enmap.write_map(pathOut+"mask_foot_planck_ps_car.fits", carFullMask)
enmap.write_map(pathOut+"tilec_map.fits", tilecMap)


