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




# Difference between 150GHz reconv to the 90GHz beam and the 90GHz map
# for null test (tSZ and kSZ estimators)
pathMap1 = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvto90.fits"
pathMap2 = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_map.fits"
pathMask = "./output/cmb_map/pact20200228/" + "mask_full_foot_gal60_ps.fits"
pathDirOut = "./output/cmb_map/planck_act_coadd_2020_02_28/"
if not os.path.exists(pathDirOut):
    os.makedirs(pathDirOut)
pathOut = pathDirOut + "act_planck_s08_s18_cmb_f150reconvto90_minus_f090_daynight_map.fits"

# save the difference map
diffMap = enmap.read_map(pathMap1)
diffMap -= enmap.read_map(pathMap2)
enmap.write_map(pathOut, diffMap)
# copy the mask
copyfile(pathMask, pathDirOut + "mask_full_foot_gal_ps.fits")



# Difference between TileC y and TileC y no CIB
# for CIB null test (tSZ estimator)
pathMap1 = "./output/cmb_map/tilec_pact_y_v1.2.0/" + "tilec_reconv14_map.fits"
pathMap2 = "./output/cmb_map/tilec_pact_ynocib_v1.2.0/" + "tilec_reconv14_map.fits"
pathMask = "./output/cmb_map/tilec_pact_y_v1.2.0/" + "mask_full_foot_gal_ps.fits"
pathDirOut = "./output/cmb_map/tilec_pact_yminusynocib_v1.2.0/"
if not os.path.exists(pathDirOut):
    os.makedirs(pathDirOut)
pathOut = pathDirOut + "tilec_reconv14_map.fits"

# save the difference map
diffMap = enmap.read_map(pathMap1)
diffMap -= enmap.read_map(pathMap2)
enmap.write_map(pathOut, diffMap)
# copy the mask
copyfile(pathMask, pathDirOut + "mask_full_foot_gal_ps.fits")


