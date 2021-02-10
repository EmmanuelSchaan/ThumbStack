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


mapId = sys.argv[1]

pathOut = "./output/cmb_map/tilec_pact_" + mapId + "_v1.2.0/"
pathFig = "./figures/cmb_map/tilec_pact_" + mapId + "_v1.2.0/"

if mapId=="cmbksz":
   pathMapD56 = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_deep56_cmb_map_v1.2.0_joint.fits"
   pathMapBN = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_boss_cmb_map_v1.2.0_joint.fits"
elif mapId=="cmbksznoy":
   pathMapD56 = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_deep56_cmb_deprojects_comptony_map_v1.2.0_joint.fits"
   pathMapBN = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_boss_cmb_deprojects_comptony_map_v1.2.0_joint.fits"
elif mapId=="cmbksznocib":
   pathMapD56 = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_deep56_cmb_deprojects_cib_map_v1.2.0_joint.fits"
   pathMapBN = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_boss_cmb_deprojects_cib_map_v1.2.0_joint.fits"
elif mapId=="y": 
   pathMapD56 = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_deep56_comptony_map_v1.2.0_joint.fits"
   pathMapBN = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_boss_comptony_map_v1.2.0_joint.fits"
elif mapId=="ynocib":
   pathMapD56 = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_deep56_comptony_deprojects_cib_map_v1.2.0_joint.fits"
   pathMapBN = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_boss_comptony_deprojects_cib_map_v1.2.0_joint.fits"
elif mapId=="ynocmb":
   pathMapD56 = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_deep56_comptony_deprojects_cmb_map_v1.2.0_joint.fits"
   pathMapBN = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_boss_comptony_deprojects_cmb_map_v1.2.0_joint.fits"




#########################################################################
# Create folders for output and figures if needed

if not os.path.exists(pathFig):
    os.makedirs(pathFig)
if not os.path.exists(pathOut):
    os.makedirs(pathOut)


#########################################################################
print "Read maps and masks"

# Deep56
d56 = enmap.read_map(pathMapD56)
pathTilecMaskD56 = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/deep56_tilec_mask.fits"   # no hit count: use the mask
tilecMaskD56 = enmap.read_map(pathTilecMaskD56)

# BOSS North
bn = enmap.read_map(pathMapBN)
pathTilecMaskBN = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/boss_tilec_mask.fits"   # no hit count: use the mask
tilecMaskBN = enmap.read_map(pathTilecMaskBN)

print("Map properties:")
print("nTQU, nY, nX = "+str(d56.shape))
print("WCS attributes: "+str(d56.wcs))


#########################################################################
# ACT point source mask

pathPSMask = "/global/cscratch1/sd/eschaan/project_ucsc/data/act_ps_mask/source_mask_s16_simonecat_sn5_cross_20171105.fits"
psMask = enmap.read_map(pathPSMask)


#########################################################################
# combine Deep56 and BOSS North maps, on the total ACT footprint
# reproject the D56 and BN maps to the overall ACT footprint
# then crop the map to reduce size

# combine the masks
tilecMaskD56 = enmap.project(tilecMaskD56, psMask.shape, psMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
tilecMaskBN = enmap.project(tilecMaskBN, psMask.shape, psMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
tilecMask = tilecMaskD56 + tilecMaskBN
# crop the geometry, and reproject the partial masks on it
tilecMask = enmap.autocrop(tilecMask)
tilecMaskD56 = enmap.project(tilecMaskD56, tilecMask.shape, tilecMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
tilecMaskBN = enmap.project(tilecMaskBN, tilecMask.shape, tilecMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)

# combine the maps: use the cropped geometry
d56 = enmap.project(d56, tilecMask.shape, tilecMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
bn = enmap.project(bn, tilecMask.shape, tilecMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
tilecMap = d56 * tilecMaskD56 + bn * tilecMaskBN

# also reproject the PS mask
psMask = enmap.project(psMask, tilecMask.shape, tilecMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)


# Plot the combined maps to check
plot=enplot.plot(tilecMask, grid=True)
enplot.write(pathFig+"foot_mask", plot)

plot=enplot.plot(tilecMap, grid=True)
enplot.write(pathFig+"map", plot)

# Save the map now, to clear memory
print "saving map to " + pathOut + "tilec_maps.fits"
enmap.write_map(pathOut+"tilec_map.fits", tilecMap)


#########################################################################
# Mask the Milky Way with a Planck mask

# Read the Planck Galactic mask
pathPlanckMask = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_galactic_mask/HFI_Mask_GalPlane-apo0_2048_R2.00.fits"
# field, fsky:
#0, 0.20
#1, 0.40
#2, 0.60
#3, 0.70
#4, 0.80
#5, 0.90
#6, 0.97
#7, 0.99
#planckMask = hp.read_map(pathPlanckMask, field=3)
hpPlanckMaskGal = hp.read_map(pathPlanckMask, field=2)

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

plot=enplot.plot(carFullMask, grid=True)
enplot.write(pathFig+"foot_planck_ps_mask", plot)

#########################################################################
# save map and mask

# save the masks
enmap.write_map(pathOut+"mask_foot_car.fits", tilecMask)
enmap.write_map(pathOut+"mask_ps_car.fits", psMask)
enmap.write_map(pathOut+"mask_planck_car.fits", planckMask)
enmap.write_map(pathOut+"mask_full_foot_gal_ps.fits", carFullMask)


