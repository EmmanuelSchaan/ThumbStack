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

# copy rotfuncs.py somewhere on your python path,
# so you can import it
import rotfuncs


# to use matplotlib and latex with jupyter
from matplotlib import rc
rc('text', usetex=True)
rc('font', family='serif')
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt


# In[2]:


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


# In[3]:


# CMB power spectra for reference
cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)
flensedTT = lambda l: cmb1_4.flensedTT(l) * cmb1_4.fbeam(l)**2
fdetectorNoise = lambda l: cmb1_4.fdetectorNoise(l) * cmb1_4.fbeam(l)**2
ftotal = lambda l: cmb1_4.ftotal(l) * cmb1_4.fbeam(l)**2


# In[4]:


# Path to output and figures
pathOut = "./output/cmb_map/compare_maps"
pathFig = "./figures/cmb_map/compare_maps"


# In[5]:


# create output/figure folders if needed
if not os.path.exists(pathOut):
    os.makedirs(pathOut)
if not os.path.exists(pathFig):
    os.makedirs(pathFig)


# # Read maps and compute power spectra

# ## TileC CMB maps, reconvolved to 1.4'

# In[6]:


# TileC CMB: D56 + BN
pathMap = "./output/cmb_map/tilec_pact_cmbksz/" + "tilec_reconv14_map.fits" 
pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz/"+"mask_foot_planck_ps_car.fits"

# read maps
tMap = enmap.read_map(pathMap)
tMask = enmap.read_map(pathCarFullMask)

# convert maps to healpix
hTMap = enmap.to_healpix(tMap)
hTMask = enmap.to_healpix(tMask)
del tMap, tMask

# plot masked map
hp.mollview(hTMap * hTMask, title='TileC BN + D56')

# compute power spectra
lCen, ClTilec, sClTilec = powerSpectrum(hTMap, mask=hTMask, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
del hTMap, hTMask


# In[7]:


# TileC CMB: d56
pathMap = "./output/cmb_map/tilec_pact_cmbksz_d56/" + "tilec_reconv14_map.fits" 
pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_d56/" + "mask_foot_planck_ps_car.fits"

# read maps
tMap = enmap.read_map(pathMap)
tMask = enmap.read_map(pathCarFullMask)

# convert maps to healpix
hTMap = enmap.to_healpix(tMap)
hTMask = enmap.to_healpix(tMask)
del tMap, tMask

# plot masked map
hp.mollview(hTMap * hTMask, title='TileC D56')

# compute power spectra
lCen, ClTilecD56, sClTilecD56 = powerSpectrum(hTMap, mask=hTMask, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
del hTMap, hTMask


# In[8]:


# CMB: BOSS N
pathMap = "./output/cmb_map/tilec_pact_cmbksz_boss/" + "tilec_reconv14_map.fits" 
pathCarFullMask = "./output/cmb_map/tilec_pact_cmbksz_boss/" + "mask_foot_planck_ps_car.fits"

# read maps
tMap = enmap.read_map(pathMap)
tMask = enmap.read_map(pathCarFullMask)

# convert maps to healpix
hTMap = enmap.to_healpix(tMap)
hTMask = enmap.to_healpix(tMask)
del tMap, tMask

# plot masked map
hp.mollview(hTMap * hTMask, title='TileC BN')

# compute power spectra
lCen, ClTilecBN, sClTilecBN = powerSpectrum(hTMap, mask=hTMask, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
del hTMap, hTMask


# ## PACT maps

# In[9]:


# Masks to be used on all PACT maps, with different footprints
# All maps include point sources

# Full AdvACT
pathMask = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits"
mask = enmap.read_map(pathMask)[0]
hMask = enmap.to_healpix(mask)
del mask
hp.mollview(hMask, title='AdvACT')

# TileC: BN + D56
pathMaskTilec = "./output/cmb_map/tilec_pact_cmbksz/"+"mask_foot_planck_ps_car.fits"
maskTilec = enmap.read_map(pathMaskTilec)
hMaskTilec = enmap.to_healpix(maskTilec)
del maskTilec
hp.mollview(hMaskTilec, title='TileC BN+D56')

# TileC: BN
pathMaskTilecBN = "./output/cmb_map/tilec_pact_cmbksz_boss/" + "mask_foot_planck_ps_car.fits"
maskTilecBN = enmap.read_map(pathMaskTilecBN)
hMaskTilecBN = enmap.to_healpix(maskTilecBN)
del maskTilecBN
hp.mollview(hMaskTilecBN, title='TileC BN')

# TileC: D56
pathMaskTilecD56 = "./output/cmb_map/tilec_pact_cmbksz_d56/" + "mask_foot_planck_ps_car.fits"
maskTilecD56 = enmap.read_map(pathMaskTilecD56)
hMaskTilecD56 = enmap.to_healpix(maskTilecD56)
del maskTilecD56
hp.mollview(hMaskTilecD56, title='TileC D56')


# In[10]:


# PACT 150GHz
pathMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_night_map.fits"

tMap = enmap.read_map(pathMap)
hMap = enmap.to_healpix(tMap)
del tMap

# plot masked map
hp.mollview(hMap * hMask, title='PACT 150')

# compute power spectrum on various patches
lCen, ClPACT150, sClPACT150 = powerSpectrum(hMap, mask=hMask, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
lCen, ClPACT150BND56, sClPACT150BND56 = powerSpectrum(hMap, mask=hMaskTilec, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
lCen, ClPACT150BN, sClPACT150BN = powerSpectrum(hMap, mask=hMaskTilecBN, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
lCen, ClPACT150D56, sClPACT150D56 = powerSpectrum(hMap, mask=hMaskTilecD56, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
del hTMap, hMask


# In[ ]:


np.shape(hMap * hMask)


# In[ ]:





# In[ ]:


# PACT 90GHz
pathMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_night_map.fits"

tMap = enmap.read_map(pathMap)
hMap = enmap.to_healpix(tMap)
del tMap

# plot masked map
hp.mollview(hMap * hMask, title='PACT 90')

# compute power spectrum on various patches
# compute power spectrum on various patches
lCen, ClPACT90, sClPACT90 = powerSpectrum(hMap, mask=hMask, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
lCen, ClPACT90BND56, sClPACT90BND56 = powerSpectrum(hMap, mask=hMaskTilec, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
lCen, ClPACT90BN, sClPACT90BN = powerSpectrum(hMap, mask=hMaskTilecBN, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
lCen, ClPACT90D56, sClPACT90D56 = powerSpectrum(hMap, mask=hMaskTilecD56, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
del hMap, hMask


# In[ ]:


# PACT 220GHz
pathMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f220_night_map.fits"

tMap = enmap.read_map(pathMap)
hMap = enmap.to_healpix(tMap)
del tMap

# plot masked map
hp.mollview(hMap * hMask, title='PACT 220')

# compute power spectrum on various patches
# compute power spectrum on various patches
lCen, ClPACT220, sClPACT220 = powerSpectrum(hMap, mask=hMask, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
lCen, ClPACT220BND56, sClPACT220BND56 = powerSpectrum(hMap, mask=hMaskTilec, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
lCen, ClPACT220BN, sClPACT220BN = powerSpectrum(hMap, mask=hMaskTilecBN, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
lCen, ClPACT220D56, sClPACT220D56 = powerSpectrum(hMap, mask=hMaskTilecD56, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)
del hMap, hMask


# In[ ]:


# del hMask, hMaskTilec, hMaskTilecBN, hMaskTilecD56


# # Save all power spectra to file

# In[ ]:


data = np.zeros((len(lCen), 20))
data[:,0] = lCen
data[:,1] = ClTilec
data[:,2] = ClTilecBN
data[:,3] = ClTilecD56
#
data[:,4] = ClPact150
data[:,5] = ClPact150BND56
data[:,6] = ClPact150BN
data[:,7] = ClPact150D56
#
data[:,8] = ClPact90
data[:,9] = ClPact90BND56
data[:,10] = ClPact90BN
data[:,11] = ClPact90D56
#
data[:,12] = ClPact220
data[:,13] = ClPact220BND56
data[:,14] = ClPact220BN
data[:,15] = ClPact220D56
#
np.savetxt(pathOut + "compare_power_all.txt", data)


# # Compare power spectra

# ## Tilec VS 150GHz

# In[ ]:


# Tilec VS 150GHz
fig=plt.figure(0, figsize=(14, 10))
ax=fig.add_subplot(111)
#
f = lCen**2 / (2.*np.pi)
#
# BN + D56
ax.loglog(lCen, f * ClTilec, 'r', label=r'TileC BN+D56')
ax.loglog(lCen, f * ClPact150BND56, 'r--', label=r'PACT BN+D56')
#
# BN
ax.loglog(lCen, f * ClTilecBN, 'b', label=r'TileC BN')
ax.loglog(lCen, f * ClPact150BN, 'b--', label=r'PACT BN')
#
# D56
ax.loglog(lCen, f * ClTilecD56, 'g', label=r'TileC D56')
ax.loglog(lCen, f * ClPact150D56, 'g--', label=r'PACT D56')
#
ax.loglog(lCen, f * ClPact150, 'y', label=r'PACT')
#
ax.loglog(lCen, f * np.array(map(flensedTT, lCen)), 'k', label=r'lensed CMB')
#
ax.set_title(r'TileC VS PACT150')
ax.legend(loc=3)
ax.set_ylim((5.e1, 1.e4))
ax.set_xlim((10., 7.e3))
ax.set_xlabel(r'$\ell$')
ax.set_xlabel(r'$\ell^2 C_\ell$')
#
fig.savefig(pathFig + "compare_power_tilec_vs_pact220.pdf")

plt.show()


# ## Tilec VS 90GHz

# In[ ]:


# Tilec VS 90GHz
fig=plt.figure(0, figsize=(14, 10))
ax=fig.add_subplot(111)
#
f = lCen**2 / (2.*np.pi)
#
# BN + D56
ax.loglog(lCen, f * ClTilec, 'r', label=r'TileC BN+D56')
ax.loglog(lCen, f * ClPact90BND56, 'r--', label=r'PACT BN+D56')
#
# BN
ax.loglog(lCen, f * ClTilecBN, 'b', label=r'TileC BN')
ax.loglog(lCen, f * ClPact90BN, 'b--', label=r'PACT BN')
#
# D56
ax.loglog(lCen, f * ClTilecD56, 'g', label=r'TileC D56')
ax.loglog(lCen, f * ClPact90D56, 'g--', label=r'PACT D56')
#
ax.loglog(lCen, f * ClPact90, 'y', label=r'PACT')
#
ax.loglog(lCen, f * np.array(map(flensedTT, lCen)), 'k', label=r'lensed CMB')
#
ax.set_title(r'TileC VS PACT90')
ax.legend(loc=3)
ax.set_ylim((5.e1, 1.e4))
ax.set_xlim((10., 7.e3))
ax.set_xlabel(r'$\ell$')
ax.set_xlabel(r'$\ell^2 C_\ell$')
#
fig.savefig(pathFig + "compare_power_tilec_vs_pact90.pdf")

plt.show()


# ## Tilec VS 220GHz

# In[ ]:


# Tilec VS 90GHz
fig=plt.figure(0, figsize=(14, 10))
ax=fig.add_subplot(111)
#
f = lCen**2 / (2.*np.pi)
#
# BN + D56
ax.loglog(lCen, f * ClTilec, 'r', label=r'TileC BN+D56')
ax.loglog(lCen, f * ClPact90BND56, 'r--', label=r'PACT BN+D56')
#
# BN
ax.loglog(lCen, f * ClTilecBN, 'b', label=r'TileC BN')
ax.loglog(lCen, f * ClPact90BN, 'b--', label=r'PACT BN')
#
# D56
ax.loglog(lCen, f * ClTilecD56, 'g', label=r'TileC D56')
ax.loglog(lCen, f * ClPact90D56, 'g--', label=r'PACT D56')
#
ax.loglog(lCen, f * ClPact90, 'y', label=r'PACT')
#
ax.loglog(lCen, f * np.array(map(flensedTT, lCen)), 'k', label=r'lensed CMB')
#
ax.set_title(r'TileC VS PACT220')
ax.legend(loc=3)
ax.set_ylim((5.e1, 1.e4))
ax.set_xlim((10., 7.e3))
ax.set_xlabel(r'$\ell$')
ax.set_xlabel(r'$\ell^2 C_\ell$')
#
fig.savefig(pathFig + "compare_power_tilec_vs_pact220.pdf")

plt.show()


# In[ ]:




