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

mapId = sys.argv[1]

pathOut = "./output/cmb_map/tilec_pact_" + mapId + "_v1.2.0/"
pathFig = "./figures/cmb_map/tilec_pact_" + mapId + "_v1.2.0/"


# CMB
if mapId=="cmbksz":
   pathMap = pathOut+"tilec_map.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_v1.2.0/"+"mask_foot_planck_ps_car.fits"
elif mapId=="cmbksz_d56":
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_deep56/tilec_single_tile_deep56_cmb_map_v1.2.0_joint.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_d56_v1.2.0/" + "mask_foot_planck_ps_car.fits"
elif mapId=="cmbksz_boss":
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_boss/tilec_single_tile_boss_cmb_map_v1.2.0_joint.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_boss_v1.2.0/" + "mask_foot_planck_ps_car.fits"

# CMB no y
elif mapId=="cmbksznoy":
   pathMap = pathOut+"tilec_map.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_v1.2.0/"+"mask_foot_planck_ps_car.fits"
elif mapId=="cmbksznoy_d56":
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_deep56/tilec_single_tile_deep56_cmb_deprojects_comptony_map_v1.2.0_joint.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_d56_v1.2.0/" + "mask_foot_planck_ps_car.fits"
elif mapId=="cmbksznoy_boss":
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_boss/tilec_single_tile_boss_cmb_deprojects_comptony_map_v1.2.0_joint.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_boss_v1.2.0/" + "mask_foot_planck_ps_car.fits"

# CMB no CIB
elif mapId=="cmbksznocib":
   pathMap = pathOut+"tilec_map.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_v1.2.0/"+"mask_foot_planck_ps_car.fits"
elif mapId=="cmbksznocib_d56":
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_deep56/tilec_single_tile_deep56_cmb_deprojects_cib_map_v1.2.0_joint.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_d56_v1.2.0/" + "mask_foot_planck_ps_car.fits"
elif mapId=="cmbksznoy_boss":
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_boss/tilec_single_tile_boss_cmb_deprojects_cib_map_v1.2.0_joint.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_boss_v1.2.0/" + "mask_foot_planck_ps_car.fits"

# y
elif mapId=="y":
   pathMap = pathOut+"tilec_map.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_v1.2.0/"+"mask_foot_planck_ps_car.fits"
elif mapId=="y_d56":
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_deep56/tilec_single_tile_deep56_comptony_map_v1.2.0_joint.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_d56_v1.2.0/" + "mask_foot_planck_ps_car.fits"
elif mapId=="y_boss":
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_boss/tilec_single_tile_boss_comptony_map_v1.2.0_joint.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_boss_v1.2.0/" + "mask_foot_planck_ps_car.fits"

# y no CIB
elif mapId=="ynocib":
   pathMap = pathOut+"tilec_map.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_v1.2.0/"+"mask_foot_planck_ps_car.fits"
elif mapId=="ynocib_d56":
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_deep56/tilec_single_tile_deep56_comptony_deprojects_cib_map_v1.2.0_joint.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_d56_v1.2.0/" + "mask_foot_planck_ps_car.fits"
elif mapId=="ynocib_boss":
   pathMap = "/global/cscratch1/sd/msyriac/data/depot/tilec/v1.2.0_20200324/map_v1.2.0_joint_boss/tilec_single_tile_boss_comptony_deprojects_cib_map_v1.2.0_joint.fits"
   pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_boss_v1.2.0/" + "mask_foot_planck_ps_car.fits"




##########################################################################
# # Do everything

# In[27]:


# create output/figure folders if needed
if not os.path.exists(pathOut):
    os.makedirs(pathOut)
if not os.path.exists(pathFig):
    os.makedirs(pathFig)


# In[28]:


# read maps
tilec = enmap.read_map(pathMap)
carFullMask = enmap.read_map(pathCarFullMask)


# # Reconvolve it to a 1.4' beam

# In[30]:


# beam sigmas in radians
s14 = 1.4 * np.pi/(180.*60.)
s16 = 1.4 * np.pi/(180.*60.)

# do the reconvolution in Fourier space
tilecReconvF = enmap.fft(tilec)
l2 = np.sum(tilec.lmap()**2,0)
tilecReconvF *= np.exp(-0.5*l2*(s14**2 - s16**2))
tilecReconv = enmap.ifft(tilecReconvF).real


# In[31]:


del tilec



# In[33]:


# save it to file
enmap.write_map(pathOut+"tilec_reconv14_map.fits", tilecReconv)


# # Convert maps to healpy (to take their power spectra)
'''
# In[34]:


# convert maps to healpix
hTilecReconv = enmap.to_healpix(tilecReconv)
hCarFullMask = enmap.to_healpix(carFullMask)
print hp.get_nside(hTilecReconv)


# In[35]:


del tilecReconv
del carFullMask



# # Compare power spectrum to single frequency maps

# In[37]:


# Binned power spectrum wrapper for healpy

def powerSpectrum(hMap, mask=None, theory=[], fsCl=None, nBins=51, lRange=None, plot=False, path="./figures/tests/test_power.pdf", save=True):
   """Compute the power spectrum of a healpix map.
   """
   
   nSide = hp.get_nside(hMap)
   if mask is not None:
      hMap *= mask

   # define ell bins
   lMax = 3 * nSide - 1
   if lRange is None:
      lEdges = np.logspace(np.log10(1.), np.log10(lMax), nBins, 10.)
   else:
      lEdges = np.logspace(np.log10(lRange[0]), np.log10(lRange[-1]), nBins, 10.)

   # Unbinned power spectrum
   power = hp.anafast(hMap)
   power = np.nan_to_num(power)

   # corresponding ell values
   ell = np.arange(len(power))

   # Bin the power spectrum
   Cl, lEdges, binIndices = stats.binned_statistic(ell, power, statistic='mean', bins=lEdges)
   
   # correct for fsky from the mask
   if mask is not None:
      fsky = np.sum(mask) / len(mask)
      Cl /= fsky

   # bin centers
   lCen, lEdges, binIndices = stats.binned_statistic(ell, ell, statistic='mean', bins=lEdges)
   # when bin is empty, replace lCen by a naive expectation
   lCenNaive = 0.5*(lEdges[:-1]+lEdges[1:])
   lCen[np.where(np.isnan(lCen))] = lCenNaive[np.where(np.isnan(lCen))]
   # number of modes
   Nmodes, lEdges, binIndices = stats.binned_statistic(ell, 2*ell+1, statistic='sum', bins=lEdges)
   Nmodes = np.nan_to_num(Nmodes)

   # 1sigma uncertainty on Cl
   if fsCl is None:
      sCl = Cl*np.sqrt(2)
   else:
      sCl = np.array(map(fsCl, lCen))
   sCl /= np.sqrt(Nmodes)
   sCl = np.nan_to_num(sCl)

   if plot:
      factor = lCen**2  # 1.
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.errorbar(lCen, factor*Cl, yerr=factor* sCl, c='b', fmt='.')
      ax.errorbar(lCen, -factor*Cl, yerr=factor* sCl, c='r', fmt='.')
      #
      for f in theory:
         L = np.logspace(np.log10(1.), np.log10(np.max(ell)), 201, 10.)
         ClExpected = np.array(map(f, L))
         ax.plot(L, L**2 * ClExpected, 'k')
      #
      ax.set_xscale('log', nonposx='clip')
      ax.set_yscale('log', nonposy='clip')
      #ax.set_xlim(1.e1, 4.e4)
      #ax.set_ylim(1.e-5, 2.e5)
      ax.set_xlabel(r'$\ell$')
      #ax.set_ylabel(r'$\ell^2 C_\ell$')
      ax.set_ylabel(r'$C_\ell$')
      #
      if save==True:
         print "saving plot to "+path
         fig.savefig(path, bbox_inches='tight')
#          fig.clf()
      else:
         plt.show()

   return lCen, Cl, sCl


# In[38]:


# CMB power spectra for reference
cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)
flensedTT = lambda l: cmb1_4.flensedTT(l) * cmb1_4.fbeam(l)**2
fdetectorNoise = lambda l: cmb1_4.fdetectorNoise(l) * cmb1_4.fbeam(l)**2
ftotal = lambda l: cmb1_4.ftotal(l) * cmb1_4.fbeam(l)**2


# In[39]:


# compare with the Sigurd map power spectrum, on the same mask
pathPact150 = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_night_map.fits"
pact150 = enmap.read_map(pathPact150)[0]
hPact150 = enmap.to_healpix(pact150)

del pact150

lCen, ClPact150, sClPact150 = powerSpectrum(hPact150, mask=hCarFullMask, theory=[flensedTT, ftotal], fsCl=None, nBins=101, lRange=None, plot=False, path=pathFig+"test", save=False)


# In[40]:


# TileC power spectrum
lCen, ClTilecReconv, sClTilecReconv = powerSpectrum(hTilecReconv, mask=hCarFullMask, theory=[flensedTT, ftotal], fsCl=None, nBins=101, lRange=None, plot=False, path=pathFig+"test", save=False)

del hTilecReconv

# # Save it to file
# data = np.zeros((len(lCen), 3))
# data[:,0] = lCen
# data[:,1] = np.nan_to_num(Cl)
# data[:,2] = np.nan_to_num(sCl)
# np.savetxt(pathIn+"f150_power_T_masked.txt", data)


# In[42]:


fig=plt.figure(0, figsize=(14, 10))
ax=fig.add_subplot(111)
#
f = lCen**2 / (2.*np.pi)
#
ax.loglog(lCen, f * ClTilecReconv, label=r'TileC 1.4')
ax.loglog(lCen, f * ClPact150, label=r'PACT150')
#
ax.loglog(lCen, f * np.array(map(flensedTT, lCen)), 'k', label=r'lensed CMB')
#
ax.legend(loc=3)
ax.set_ylim((5.e1, 1.e4))
ax.set_xlim((10., 7.e3))
ax.set_xlabel(r'$\ell$')
ax.set_xlabel(r'$\ell^2 C_\ell$')
#
fig.savefig(pathFig + "compare_power_tilec_vs_pact150.pdf")
fig.clf()
#plt.show()


# In[ ]:
'''



