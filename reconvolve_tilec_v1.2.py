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
# reconvolve the no-deproj tilec maps
# to the beam of the deproj tilec maps

for mapId in ['cmbksz', 'y']:

   pathOut = "./output/cmb_map/tilec_pact_" + mapId + "_v1.2.0/"
   pathFig = "./figures/cmb_map/tilec_pact_" + mapId + "_v1.2.0/"

   # create output/figure folders if needed
   if not os.path.exists(pathOut):
       os.makedirs(pathOut)
   if not os.path.exists(pathFig):
       os.makedirs(pathFig)

   pathMap = pathOut+"tilec_map.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_v1.2.0/"+"mask_full_foot_gal_ps.fits"
   fwhmOld = 1.6 
   fwhmNew = 2.4

   # beam sigmas in radians
   # convert from fwhm to sigma
   sOld = fwhmOld * np.pi/(180.*60.) / np.sqrt(8.*np.log(2.))
   sNew = fwhmNew * np.pi/(180.*60.) / np.sqrt(8.*np.log(2.))

   # read maps
   tilec = enmap.read_map(pathMap)

   # do the reconvolution in Fourier space
   tilecReconvF = enmap.fft(tilec)
   l2 = np.sum(tilec.lmap()**2,0)
   tilecReconvF *= np.exp(-0.5*l2*(sNew**2 - sOld**2))
   tilecReconv = enmap.ifft(tilecReconvF).real

   # save it to file
   enmap.write_map(pathOut+"tilec_reconv"+str(fwhmNew)+"_map.fits", tilecReconv)






