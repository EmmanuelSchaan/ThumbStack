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


from shutil import copyfile



###################################################################################3
# Difference between 150GHz reconv to the 90GHz beam and the 90GHz map
# for null test (tSZ and kSZ estimators)

pathMap1 = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvto90.fits"
pathMap2 = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f090_daynight_map.fits"
pathMask = "./output/cmb_map/pact20200228_r2/" + "mask_full_foot_gal60_ps.fits"
pathDirOut = "./output/cmb_map/planck_act_coadd_2020_02_28_r2/"
if not os.path.exists(pathDirOut):
    os.makedirs(pathDirOut)
pathOut = pathDirOut + "act_planck_s08_s18_cmb_f150reconvto90_minus_f090_daynight_map.fits"

# save the difference map
diffMap = enmap.read_map(pathMap1)
diffMap -= enmap.read_map(pathMap2)
enmap.write_map(pathOut, diffMap)
# copy the mask
copyfile(pathMask, pathDirOut + "mask_full_foot_gal_ps.fits")


###################################################################################3
# Difference between TileC y and TileC y no CIB
# for CIB null test (tSZ estimator)

pathMap1 = "./output/cmb_map/tilec_pact_y_v1.2.0/" + "tilec_reconv2.4_map.fits"
pathMap2 = "./output/cmb_map/tilec_pact_ynocib_v1.2.0/" + "tilec_map.fits"
pathMask = "./output/cmb_map/tilec_pact_y_v1.2.0/" + "mask_full_foot_gal_ps.fits"
pathDirOut = "./output/cmb_map/tilec_pact_yminusynocib_v1.2.0/"
if not os.path.exists(pathDirOut):
    os.makedirs(pathDirOut)
pathOut = pathDirOut + "diff_map.fits"

# save the difference map
diffMap = enmap.read_map(pathMap1)
diffMap -= enmap.read_map(pathMap2)
enmap.write_map(pathOut, diffMap)
# copy the mask
copyfile(pathMask, pathDirOut + "mask_full_foot_gal_ps.fits")



###################################################################################3
# Difference between PACT 150 and TileC cmbksz 
# for kSZ null test

pathMap1 = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvtotilec.fits"
pathMap2 = "./output/cmb_map/tilec_pact_cmbksz_v1.2.0/" + "tilec_map.fits"
pathMask = "./output/cmb_map/tilec_pact_cmbksz_v1.2.0/" + "mask_full_foot_gal_ps.fits"
pathDirOut = "./output/cmb_map/pactf150daynight20200228maskgal60r2_minus_tilec_pact_cmbksz/"
if not os.path.exists(pathDirOut):
    os.makedirs(pathDirOut)
pathOut = pathDirOut + "diff_map.fits"

# read the maps
map1 = enmap.read_map(pathMap1)[0]
map2 = enmap.read_map(pathMap2)
# reproject to the smaller TileC geometry
map1 = enmap.project(map1, map2.shape, map2.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
# take the difference
diffMap = map1 - map2
enmap.write_map(pathOut, diffMap)
# copy the mask
copyfile(pathMask, pathDirOut + "mask_full_foot_gal_ps.fits")


###################################################################################3
# conversion from y to muK at 150 GHz
Tcmb = 2.726   # K
h = 6.63e-34   # SI
kB = 1.38e-23  # SI
def f(nu):
   """frequency dependence for tSZ temperature
   """
   x = h*nu/(kB*Tcmb)
   return x*(np.exp(x)+1.)/(np.exp(x)-1.) -4.
yTomuK = f(150.e9) * Tcmb * 1.e6  # [muK * sr]


###################################################################################3
# Difference between PACT 150 and TileC y 
# after converting the y map to muK at 150GHz
# for tSZ pipeline test (both maps will have some dust; the dust level may be different in both maps though)

pathMap1 = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvtotilec.fits"
pathMap2 = "./output/cmb_map/tilec_pact_y_v1.2.0/" + "tilec_map.fits"
pathMask = "./output/cmb_map/tilec_pact_y_v1.2.0/" + "mask_full_foot_gal_ps.fits"
pathDirOut = "./output/cmb_map/pactf150daynight20200228maskgal60r2_minus_tilec_pact_ymuk/"
if not os.path.exists(pathDirOut):
    os.makedirs(pathDirOut)
pathOut = pathDirOut + "diff_map.fits"

# read the maps
map1 = enmap.read_map(pathMap1)[0]
map2 = enmap.read_map(pathMap2) * yTomuK
# reproject to the smaller TileC geometry
map1 = enmap.project(map1, map2.shape, map2.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
# take the difference
diffMap = map1 - map2
enmap.write_map(pathOut, diffMap)
# copy the mask
copyfile(pathMask, pathDirOut + "mask_full_foot_gal_ps.fits")


###################################################################################3
# Difference between PACT 150 and TileC CMB no CIB
# to check for dust contamination
# for kSZ foreground test

pathMap1 = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvtotilecdeproj.fits"
pathMap2 = "./output/cmb_map/tilec_pact_cmbksznocib_v1.2.0/" + "tilec_map.fits"
pathMask = "./output/cmb_map/tilec_pact_cmbksznocib_v1.2.0/" + "mask_full_foot_gal_ps.fits"
pathDirOut = "./output/cmb_map/pactf150daynight20200228maskgal60r2_minus_tilec_pact_cmbksznocib/"
if not os.path.exists(pathDirOut):
    os.makedirs(pathDirOut)
pathOut = pathDirOut + "diff_map.fits"

# read the maps
map1 = enmap.read_map(pathMap1)[0]
map2 = enmap.read_map(pathMap2)
# reproject to the smaller TileC geometry
map1 = enmap.project(map1, map2.shape, map2.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
# take the difference
diffMap = map1 - map2
enmap.write_map(pathOut, diffMap)
# copy the mask
copyfile(pathMask, pathDirOut + "mask_full_foot_gal_ps.fits")


###################################################################################3
# Difference between 150GHz daynight and 150GHz night-only,
# for null tests

pathMap1 = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits"
pathMap2 = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_night_map.fits"
pathMask = "./output/cmb_map/pact20200228_r2/" + "mask_full_foot_gal60_ps.fits"
pathDirOut = "./output/cmb_map/planck_act_coadd_2020_02_28_r2/"
if not os.path.exists(pathDirOut):
    os.makedirs(pathDirOut)
pathOut = pathDirOut + "act_planck_s08_s18_cmb_f150_daynight_minus_night_map.fits"

# save the difference map
diffMap = enmap.read_map(pathMap1)
diffMap -= enmap.read_map(pathMap2)
enmap.write_map(pathOut, diffMap)
# copy the mask
copyfile(pathMask, pathDirOut + "mask_full_foot_gal_ps.fits")


