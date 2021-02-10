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



# # Specify map and paths

# In[2]:


mapId = sys.argv[1]

# In[3]:


if mapId=='cmbksz_d56':
   pathFig = "./figures/cmb_map/tilec_pact_cmbksz_d56_v1.2.0/"
   pathOut = "./output/cmb_map/tilec_pact_cmbksz_d56_v1.2.0/"
   pathMap = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_deep56_cmb_map_v1.2.0_joint.fits"
   pathTilecMask = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/deep56_tilec_mask.fits"
if mapId=='cmbksz_boss':
   pathFig = "./figures/cmb_map/tilec_pact_cmbksz_boss_v1.2.0/"
   pathOut = "./output/cmb_map/tilec_pact_cmbksz_boss_v1.2.0/"
   pathMap = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/tilec_single_tile_boss_cmb_map_v1.2.0_joint.fits"
   pathTilecMask = "/global/cscratch1/sd/eschaan/project_ucsc/data/tilec_v1.2.0/boss_tilec_mask.fits"


# In[4]:


# create output/figure folders if needed
if not os.path.exists(pathOut):
    os.makedirs(pathOut)
if not os.path.exists(pathFig):
    os.makedirs(pathFig)


# # Generate footprint mask from hit count

# In[5]:


# read the hit count map
footMask = enmap.read_map(pathTilecMask)
# threshold it
footMask = (footMask<>0.).astype(np.float)


# In[6]:


# # Check sky area
# fSky = np.sum(footMask) / np.prod(footMask.shape)
# print "sky area =", fSky * footMask.area() * (180./np.pi)**2, "deg2"

# # look at the footprint
# plot=enplot.plot(footMask, downgrade=16)
# enplot.show(plot)


# In[7]:


enmap.write_map(pathOut+"mask_foot.fits", footMask)


# # Planck Galactic mask

# In[8]:


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
planckMask = hp.read_map(pathPlanckMask, field=2)


# In[9]:


# # plot the Planck mask in Galactic coordinates
# fig=plt.figure(0)
# hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
# # fig.savefig(pathFig+"planck_mask_G.pdf")
# # fig.clf()
# plt.show()


# In[10]:


# rotate it to equatorial coordinates
rot = hp.Rotator(coord=['G','C'])  # Transforms galactic to equatorial coordinates
planckMask = rot.rotate_map_pixel(planckMask)


# In[11]:


# fig=plt.figure(0)
# hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
# # fig.savefig(pathFig+"planck_mask_C.pdf")
# # fig.clf()
# plt.show()


# In[12]:


# rethreshold the mask, to make it zeros and ones
planckMask = (planckMask>0.999).astype(np.float)


# In[13]:


# fig=plt.figure(0)
# hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
# # fig.savefig(pathFig+"planck_mask_C_thresh.pdf")
# # fig.clf()
# plt.show()


# In[20]:


# convert to enmap
planckMask = reproject.enmap_from_healpix(planckMask, footMask.shape, footMask.wcs, rot=None)


# In[22]:


# plot=enplot.plot(planckMask, downgrade=16, min=0., max=1.)
# enplot.show(plot)


# # Read PS mask

# In[15]:


pathPSMask = "/global/cscratch1/sd/eschaan/project_ucsc/data/act_ps_mask/source_mask_s16_simonecat_sn5_cross_20171105.fits"
psMask = enmap.read_map(pathPSMask)


# In[16]:


psMask = (psMask>0.95).astype(np.float)


# In[17]:


# project to TileC footprint
psMask = enmap.project(psMask, footMask.shape, footMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)


# In[18]:


# plot=enplot.plot(psMask, downgrade=1, min=0., max=1.)
# enplot.show(plot)


# # Full mask: footprint + Galactic + point sources

# In[23]:


fullMask = footMask * planckMask * psMask
fullMask = (fullMask>0.95).astype(np.float)


# In[25]:


# plot=enplot.plot(fullMask, downgrade=1)
# enplot.show(plot)


# In[26]:

print "Saving full mask to:"
print pathOut+"mask_full_foot_gal_ps.fits"
enmap.write_map(pathOut+"mask_full_foot_gal_ps.fits", fullMask)


# In[ ]:





# In[ ]:




