#!/usr/bin/env python
# coding: utf-8

# In[21]:

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


##########################################################################
# # Input: choose one of the TileC maps

# beam of the 150 GHz map: 1.41'
# beam of the 90 GHz map: 2.32'
initialBeamFwhm = 1.41  # [arcmin]
finalBeamFwhm = 2.32 # [arcmin]

pathMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits"
pathMask = "./output/cmb_map/pact20200228/" + "mask_full_foot_gal60_ps.fits"
pathHit = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits"

pathOutMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvto90.fits"


##########################################################################
# # Do everything

# read maps
iMap = enmap.read_map(pathMap)
mask = enmap.read_map(pathMask)


# # Reconvolve

# beam sigmas in radians
# convert from fwhm to sigma
si = initialBeamFwhm * np.pi/(180.*60.) / np.sqrt(8.*np.log(2.))
sf = finalBeamFwhm * np.pi/(180.*60.) / np.sqrt(8.*np.log(2.))

# do the reconvolution in Fourier space
mapF = enmap.fft(iMap)
l2 = np.sum(iMap.lmap()**2,0)
mapF *= np.exp(-0.5*l2*(sf**2 - si**2))
oMap = enmap.ifft(mapF).real

# save it to file
enmap.write_map(pathOutMap, oMap)


