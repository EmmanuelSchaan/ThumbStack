#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

import os

# copy rotfuncs.py somewhere on your python path,
# so you can import it
import rotfuncs


## to use matplotlib and latex with jupyter
#from matplotlib import rc
#rc('text', usetex=True)
#rc('font', family='serif')
#get_ipython().run_line_magic('matplotlib', 'inline')
#import matplotlib.pyplot as plt


# # Specify map and paths

# The 90 and 150GHz maps have the same mask, so no need to do it twice.

# In[2]:


iMap = '90'


# In[3]:


if iMap=='90':
    pathHit = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f090_daynight_ivar.fits"
if iMap=='150':
    pathHit = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits"


# Choose the Galactic mask to use: value of fsky

# In[4]:


iGalMaskFsky = sys.argv[1]
print "Galactic mask: fsky = "+iGalMaskFsky+" %"
#iGalMaskFsky = '60'


# In[5]:


if iGalMaskFsky=='20':
    iField = 0
if iGalMaskFsky=='40':
    iField = 1
if iGalMaskFsky=='60':
    iField = 2
if iGalMaskFsky=='70':
    iField = 3
if iGalMaskFsky=='80':
    iField = 4
if iGalMaskFsky=='90':
    iField = 5
if iGalMaskFsky=='97':
    iField = 6
if iGalMaskFsky=='99':
    iField = 7


# In[6]:


# Path to output and figures
pathOut = "./output/cmb_map/pact20200228_r2/"
pathFig = "./figures/cmb_map/pact20200228_r2/"


# In[7]:


# create output/figure folders if needed
if not os.path.exists(pathOut):
    os.makedirs(pathOut)
if not os.path.exists(pathFig):
    os.makedirs(pathFig)


# # Generate footprint mask from hit count

# In[ ]:


# read the hit count map
footMask = enmap.read_map(pathHit)[0]
# threshold it
footMask = (footMask<>0.).astype(np.float)


# In[ ]:


# # Check sky area
# fSky = np.sum(footMask) / np.prod(footMask.shape)
# print "sky area =", fSky * footMask.area() * (180./np.pi)**2, "deg2"

# # look at the footprint
# plot=enplot.plot(footMask, downgrade=16)
# enplot.show(plot)


# In[ ]:


enmap.write_map(pathOut+"mask_foot.fits", footMask)


# # Planck Galactic mask

# In[ ]:


# Read the Planck Galactic mask
pathPlanckMask = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_galactic_mask/HFI_Mask_GalPlane-apo0_2048_R2.00.fits"
# field, fsky:
#0, 0.20
#1, 0.40
#2, 0.60
#3, 0.70
#4, 0.80
#5, 0.90
#6, 0.97
#7, 0.99
# planckMask = hp.read_map(pathPlanckMask, field=3)
# planckMask = hp.read_map(pathPlanckMask, field=2)
planckMask = hp.read_map(pathPlanckMask, field=iField)


# In[ ]:


# # plot the Planck mask in Galactic coordinates
# fig=plt.figure(0)
# hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
# # fig.savefig(pathFig+"planck_mask_G.pdf")
# # fig.clf()
# plt.show()


# In[ ]:


# rotate it to equatorial coordinates
rot = hp.Rotator(coord=['G','C'])  # Transforms galactic to equatorial coordinates
planckMask = rot.rotate_map_pixel(planckMask)


# In[ ]:


# fig=plt.figure(0)
# hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
# # fig.savefig(pathFig+"planck_mask_C.pdf")
# # fig.clf()
# plt.show()


# In[ ]:


# rethreshold the mask, to make it zeros and ones
planckMask = (planckMask>0.999).astype(np.float)


# In[ ]:


# fig=plt.figure(0)
# hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
# # fig.savefig(pathFig+"planck_mask_C_thresh.pdf")
# # fig.clf()
# plt.show()


# In[ ]:


# convert to enmap
planckMask = reproject.enmap_from_healpix(planckMask, footMask.shape, footMask.wcs, rot=None)


# In[ ]:


# plot=enplot.plot(planckMask, downgrade=16)
# enplot.show(plot)


# # Read PS mask

# In[ ]:


pathPSMask = "/global/cscratch1/sd/eschaan/project_ucsc/data/act_ps_mask/source_mask_s16_simonecat_sn5_cross_20171105.fits"
psMask = enmap.read_map(pathPSMask)


# In[ ]:


# plot=enplot.plot(psMask, downgrade=16)
# enplot.show(plot)


# # Full mask: footprint + Galactic + point sources

# In[ ]:


fullMask = footMask * planckMask * psMask
fullMask = (fullMask>0.95).astype(np.float)


# In[ ]:


# plot=enplot.plot(fullMask, downgrade=16)
# enplot.show(plot)


# In[ ]:

print "Saving to:"
print pathOut+"mask_full_foot_gal"+iGalMaskFsky+"_ps.fits"
enmap.write_map(pathOut+"mask_full_foot_gal"+iGalMaskFsky+"_ps.fits", fullMask)


# # Check that the masked maps look good

# In[ ]:


# path90 = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_map.fits"
# path150= "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits"


# In[ ]:


# pact150 = enmap.read_map(path150)[0]
# pact90 = enmap.read_map(path90)[0]


# In[ ]:


# plot=enplot.plot(pact150, downgrade=16)
# enplot.show(plot)


# In[ ]:


# plot=enplot.plot(pact150 * (1. + 10.*fullMask)/11., downgrade=16)
# enplot.show(plot)


# In[ ]:


# plot=enplot.plot(pact90 * (1. + 10.*fullMask)/11., downgrade=16)
# enplot.show(plot)

