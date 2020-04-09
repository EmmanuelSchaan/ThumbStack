#!/usr/bin/env python
# coding: utf-8

# In[16]:


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
#from matplotlib import rc
#rc('text', usetex=True)
#rc('font', family='serif')
#get_ipython().run_line_magic('matplotlib', 'inline')
#import matplotlib.pyplot as plt


# In[17]:


# Binned power spectrum wrapper for healpy

def powerSpectrum(iMap, mask=None, theory=[], fsCl=None, nBins=51, lRange=None, plot=False, path="./figures/tests/test_power.pdf", save=True):
   """Compute the power spectrum of a healpix map.
   """
   
   hMap = iMap.copy()
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


# In[18]:


# CMB power spectra for reference
cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)
flensedTT = lambda l: cmb1_4.flensedTT(l) * cmb1_4.fbeam(l)**2
fdetectorNoise = lambda l: cmb1_4.fdetectorNoise(l) * cmb1_4.fbeam(l)**2
ftotal = lambda l: cmb1_4.ftotal(l) * cmb1_4.fbeam(l)**2


# In[19]:


# Path to output and figures
pathOut = "./output/cmb_map/coadd"
pathFig = "./figures/cmb_map/coadd"


# In[20]:


# create output/figure folders if needed
if not os.path.exists(pathOut):
    os.makedirs(pathOut)
if not os.path.exists(pathFig):
    os.makedirs(pathFig)


# # Examine noise maps

# In[21]:


# hit count at 150
pathHit150 = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits"
hit150 = enmap.read_map(pathHit150)[0]

# hit count at 90
pathHit90 = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_ivar.fits"
hit90 = enmap.read_map(pathHit90)[0]


# Sigurd confirmed that these maps are inverse var maps in 1/muK^2. They are properly normalized to be used as weights in the coadd.

# In[22]:


#plot=enplot.plot(hit150, downgrade=16, min=1.e-10, max=1.e-3)
#enplot.show(plot)

#plot=enplot.plot(hit90, downgrade=16, min=1.e-10, max=1.e-3)
#enplot.show(plot)


# Check whether the ratio of hit counts changes a lot with position. If so, this would suggest we should coadd with a pixel-dependent relative weight.
# The ratio only changes for D56 for which there is more 150 data (may matter), and a tiny chunk of BOSS North (not important).

# In[23]:


#plot=enplot.plot(np.log10((hit150+1.e-20)/(hit90+1.e-20)), downgrade=16, min=-1, max=1.)
#enplot.show(plot)


# In[24]:


# del hit150, hit90


# # PACT 90/150 daynight: maps, mask

# In[25]:


# PACT 150GHz
path150 = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits"
pact150 = enmap.read_map(path150)[0]


# In[26]:


# PACT 90GHz
path90 = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_map.fits"
pact90 = enmap.read_map(path90)[0]


# In[27]:


# full mask: footprint + point sources + Milky Way
pathMask = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits"
fullMask = enmap.read_map(pathMask)[0]


# In[28]:


#plot=enplot.plot(pact150 * fullMask, downgrade=16)
#enplot.show(plot)

#plot=enplot.plot(pact90 * fullMask, downgrade=16)
#enplot.show(plot)


# # Reconvolve the 90GHz to the same beam as the 150GHz

# Map info: Coadd 20200228
# https://phy-wiki.princeton.edu/polwiki/pmwiki.php?n=Maps.FieldMapDirectory
# - 150GHz beam: 150GHz S16 PA2
# - 90GHz beam: 90GHz S16 PA3
# 
# Beam infos:
# https://phy-wiki.princeton.edu/polwiki/pmwiki.php?n=BeamAnalysis.BeamDistributionCenter
# - 150GHz S16 PA2: Omega = 191 nsr, ie fwhm = 1.41 arcmin
# - 90GHz S16 PA3: Omega = 516 nsr, ie fwhm = 2.32 arcmin

# In[29]:


# beam sigmas in radians (convert from fwhm to sigma)
s150 = 1.41 * np.pi/(180.*60.) / np.sqrt(8.*np.log(2.))
s90 = 2.32 * np.pi/(180.*60.) / np.sqrt(8.*np.log(2.))

# do the reconvolution in Fourier space
pact90ReconvF = enmap.fft(pact90)
l2 = np.sum(pact90.lmap()**2,0)
pact90ReconvF *= np.exp(-0.5*l2*(s150**2 - s90**2))
pact90 = enmap.ifft(pact90ReconvF).real


# # Measure power spectra

# In[ ]:


#hMask = enmap.to_healpix(fullMask)

#h150 = enmap.to_healpix(pact150)
#lCen, ClPACT150, sClPACT150 = powerSpectrum(h150, mask=hMask, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)

#h90 = enmap.to_healpix(pact90)
#lCen, ClPACT90, sClPACT90 = powerSpectrum(h90, mask=hMask, theory=[flensedTT, ftotal], fsCl=None, nBins=201, lRange=None, plot=False, path=None, save=False)


# In[ ]:


#fig=plt.figure(0, figsize=(14,10))
#ax=fig.add_subplot(111)
##
#f = lCen**2 / (2.*np.pi)
##
#ax.loglog(lCen, f * ClPACT150, 'b', label=r'150')
#ax.loglog(lCen, f * ClPACT90, 'r', label=r'90')
##
#ax.loglog(lCen, f * np.array(map(cmb1_4.flensedTT, lCen)), 'k', label=r'lensed CMB')
##
#ax.legend(loc=3)
#ax.set_xlim((20., 1.e4))
#ax.set_ylim((1.e2, 1.e4))
#ax.set_xlabel(r'$\ell$')
#ax.set_ylabel(r'$\ell^2 C_\ell / (2\pi)$')
##
#fig.savefig(pathFig + "/power_90_150.pdf", bbox_inches='tight')
#
#plt.show()


# ## Check the relative calibration

# In[ ]:


#fig=plt.figure(0, figsize=(14, 10))
#ax=fig.add_subplot(111)
##
#ax.semilogx(lCen, ClPACT150 / ClPACT90, 'g', label=r'150 / 90')
## ax.loglog(lCen, ClPACT150 / np.array(map(cmb1_4.flensedTT, lCen)), 'b', label=r'150')
## ax.loglog(lCen, ClPACT90 / np.array(map(cmb1_4.flensedTT, lCen)), 'r', label=r'90')
#ax.semilogx(lCen, np.ones_like(lCen), 'k', label=r'lensed CMB')
##
#ax.legend(loc=3)
#ax.set_xlim((20., 1.e4))
#ax.set_ylim((0.8, 1.6))
#ax.set_xlabel(r'$\ell$')
#ax.set_ylabel(r'$C_\ell / C_\ell^\text{lensed CMB}$')
##
#fig.savefig(pathFig + "/ratio_power_90_150.pdf", bbox_inches='tight')
#          
#plt.show()


# ## Equal weights coadd

# In principle, we would inverse-variance weight the maps.
# Here, the power spectra only differ by 50% at most. Let's try equal weights.

# In[ ]:


#hMean = 0.5 * (h150 + h90)
#lCen, ClMean, sClMean = powerSpectrum(hMean, mask=hMask, nBins=201, lRange=None, plot=False, path=None, save=False)
#del hMean


# In[ ]:


meanPact = 0.5 * (pact150 + pact90)
hitMeanPact = np.nan_to_num(1. / (0.25 * (1./hit150 + 1./hit90)))

enmap.write_map(pathOut + "/coadd_equal_weights.fits", meanPact)
enmap.write_map(pathOut + "/coadd_equal_weights_ivar.fits", hitMeanPact)
del meanPact, hitMeanPact


# ## Approximate inverse noise variance

# In[ ]:


ivar150 = 1./1.3
ivar90 = 1./1.
w150 = ivar150 / (ivar90 + ivar150)
w90 = ivar90 / (ivar90 + ivar150)


# In[ ]:


#hIvar = w150 * h150 + w90 * h90
#lCen, ClIvar, sClIvar = powerSpectrum(hIvar, mask=hMask, nBins=201, lRange=None, plot=False, path=None, save=False)
#del hIvar


# In[ ]:


ivarPact = w150 * pact150 + w90 * pact90
hitIvarPact = np.nan_to_num(1. / (w150**2 * 1./hit150 + w90**2 * 1./hit90))

enmap.write_map(pathOut + "/coadd_approx_inv_noise_weighted.fits", ivarPact)
enmap.write_map(pathOut + "/coadd_approx_inv_noise_weighted_ivar.fits", ivarPact)
del ivarPact


# ## Weight by hit map

# The hit maps are true inverse noise variance maps, in 1/muK^2.

# In[ ]:


coaddHit = np.nan_to_num((hit150*pact150 + hit90*pact90) / (hit150 + hit90))
hitHit = np.nan_to_num(hit150 + hit90)

enmap.write_map(pathOut + "/coadd_hitmap_weighted.fits", coaddHit)
enmap.write_map(pathOut + "/coadd_hitmap_weighted_ivar.fits", hitHit)
del hitHit


# In[ ]:

#
#hCoaddHit = enmap.to_healpix(coaddHit)
#del coaddHit

#hp.mollview(hCoaddHit * hMask)

#lCen, ClCoaddHit, sClCoaddHit = powerSpectrum(hCoaddHit, mask=hMask, nBins=201, lRange=None, plot=False, path=None, save=False)
#del hCoaddHit


# # Comparing power spectra

# - There seems to be a significant gain in noise for ell>2500 from coadding.
# - The choice of weights doesn't seem to matter much.
# - There is a small ~2% loss at ell<1000 compared to 90GHz, but this is probably negligible.
# 
# --> keep the equal weights coadd

# In[ ]:


#fig=plt.figure(0, figsize=(14,10))
#ax=fig.add_subplot(111)
##
#f = lCen**2 / (2.*np.pi)
##
#ax.loglog(lCen, f * ClPACT150, 'b', label=r'150')
#ax.loglog(lCen, f * ClPACT90, 'r', label=r'90')
#ax.loglog(lCen, f * ClMean, 'g', label=r'Mean')
#ax.loglog(lCen, f * ClIvar, 'y', label=r'cheap inv noise var')
#ax.loglog(lCen, f * ClCoaddHit, 'c', label=r'hit map weighted')
##
#ax.loglog(lCen, f * np.array(map(cmb1_4.flensedTT, lCen)), 'k', label=r'lensed CMB')
##
#ax.legend(loc=3)
#ax.set_xlim((20., 1.e4))
#ax.set_ylim((1.e2, 1.e4))
#ax.set_xlabel(r'$\ell$')
#ax.set_ylabel(r'$\ell^2 C_\ell / (2\pi)$')
##
#fig.savefig(pathFig + "/power_coadds.pdf", bbox_inches='tight')
#          
#plt.show()


# In[ ]:


#fig=plt.figure(0, figsize=(14,10))
#ax=fig.add_subplot(111)
##
#ax.semilogx(lCen, ClPACT150 / ClPACT90, 'b', label=r'150')
#ax.semilogx(lCen, ClPACT90 / ClPACT90, 'r', label=r'90')
#ax.semilogx(lCen, ClMean / ClPACT90, 'g', label=r'Mean')
#ax.semilogx(lCen, ClIvar / ClPACT90, 'y', label=r'cheap inv noise var')
#ax.semilogx(lCen, ClCoaddHit / ClPACT90, 'c', label=r'hit map weighted')
##
#ax.legend(loc=3)
#ax.set_xlim((20., 1.e4))
#ax.set_ylim((0.5, 1.5))
#ax.set_xlabel(r'$\ell$')
#ax.set_ylabel(r'$C_\ell / C_\ell^{90}$')
##
#fig.savefig(pathFig + "/ratio_power_coadds.pdf", bbox_inches='tight')
#    
#plt.show()


# In[ ]:




