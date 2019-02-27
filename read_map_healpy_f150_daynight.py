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

# path for figures
pathFig = "./figures/cmb_map/planck_act_coadd_2018_08_10/f150_daynight/"

# path for output
#pathOut = "./output/cmb_map/planck_act_coadd_2018_08_10/f150_daynight/"


#########################################################################
# CMB power spectra

cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)
cmb7_3 = StageIVCMB(beam=7.3, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)

#########################################################################
# load maps

# path for maps
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/"
pathMap = pathIn + "f150_daynight_all_map_mono.fits"
pathHit = pathIn + "f150_daynight_all_div_mono.fits"
pathPSMask = pathIn + "source_mask_s16_simonecat_sn5_cross_20171105.fits"

print "Read maps"
baseMap = enmap.read_map(pathMap)
hitMap = enmap.read_map(pathHit)
psMask = enmap.read_map(pathPSMask)

print("Map properties:")
print("nTQU, nY, nX = "+str(baseMap.shape))
print("WCS attributes: "+str(baseMap.wcs))


#########################################################################
# Histograms before masking

# T histogram
myHistogram(baseMap[0], nBins=71, lim=(-10.*np.std(baseMap[0]), 10.*np.std(baseMap[0])), sigma2Theory=110.**2, path=pathFig+"hist_T_nomask.pdf", nameLatex=r'$T$ [$\mu$K]', semilogx=False, semilogy=True, doGauss=True)

# Q histogram
myHistogram(baseMap[1], nBins=71, lim=(-10.*np.std(baseMap[1]), 10.*np.std(baseMap[1])), sigma2Theory=None, path=pathFig+"hist_Q_nomask.pdf", nameLatex=r'$Q$ [$\mu$K]', semilogx=False, semilogy=True, doGauss=True)

# U histogram
myHistogram(baseMap[2], nBins=71, lim=(-10.*np.std(baseMap[2]), 10.*np.std(baseMap[2])), sigma2Theory=None, path=pathFig+"hist_U_nomask.pdf", nameLatex=r'$U$ [$\mu$K]', semilogx=False, semilogy=True, doGauss=True)

# PS mask
myHistogram(psMask, nBins=71, lim=(np.min(psMask), np.max(psMask)), sigma2Theory=None, path=pathFig+"hist_psmask.pdf", nameLatex=r'Mask values', semilogx=False, semilogy=True, doGauss=False)

# PS mask
myHistogram(hitMap, nBins=71, lim=(np.min(hitMap), np.max(hitMap)), sigma2Theory=None, path=pathFig+"hist_hit_nomask.pdf", nameLatex=r'Hit map values', semilogx=False, semilogy=True, doGauss=False)



#########################################################################
#########################################################################
# convert map to healpix

print("Convert from CAR to healpix")
hMap = enmap.to_healpix(baseMap)
hHitMap = enmap.to_healpix(hitMap)
hPSMask = enmap.to_healpix(psMask)

# re-threshold the PS mask
hPSMask = (hPSMask>0.95).astype(np.float)

# Save healpix maps to files
hp.write_map(pathIn+"healpix_f150_daynight_all_map_mono.fits", hMap, overwrite=True)
hp.write_map(pathIn+"healpix_f150_daynight_all_div_mono.fits", hHitMap, overwrite=True)
hp.write_map(pathIn+"healpix_pointsourcemask.fits", hPSMask, overwrite=True)

nSide = hp.get_nside(hMap)
print("Nside = "+str(nSide))


#########################################################################
# For debugging: downgrade the map for speed
'''
hMap = hp.ud_grade(hMap, 1024)
hHitMap = hp.ud_grade(hHitMap, 1024)
hPSMap = hp.ud_grade(hPSMap, 1024)
nSide = hp.get_nside(hMap)
'''

#########################################################################
# plot the PS mask

fig=plt.figure(0)
hp.mollview(hPSMask, fig=0, title="Point source mask", max=np.max(hPSMask), coord=None, cbar=False, unit='')
fig.savefig(pathFig+"psmask.pdf")
fig.clf()


#########################################################################
# Plot the hit count map

fig=plt.figure(0)
hp.mollview(hHitMap, fig=0, title="Hit count", max=0.5*np.max(hHitMap), coord=None, cbar=True, unit='')
fig.savefig(pathFig+"f150_hit_0.5max.pdf")
fig.clf()

fig=plt.figure(0)
logHit = np.log10(np.abs(hHitMap)+1.e-5)
hp.mollview(logHit, fig=0, title="Hit count", min=np.min(logHit), max=np.max(logHit), coord=None, cbar=False, unit='')
fig.savefig(pathFig+"f150_hit_log.pdf")
fig.clf()


#########################################################################
# Generate a footprint mask from the hit map

# copy the masked temperature map
footMask = hHitMap.copy()

# generate a mask from all the null pixels,
footMask = (footMask<>0.).astype(np.float)

fig=plt.figure(0)
hp.mollview(footMask, fig=0, title="Footprint", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"f150_footmask.pdf")
fig.clf()

# save the footprint to file
hp.write_map(pathIn+"f150_footprint_mask.fits", footMask, overwrite=True)

# read the footprint mask
footMask = hp.read_map(pathIn+"f150_footprint_mask.fits")


#########################################################################
# Mask the Milky Way with a Planck mask

# Read the Planck lensing mask
#pathPlanckLensingMask = pathIn+"/COM_Mask_Lensing_2048_R2.00.fits"
#planckMask = hp.read_map(pathPlanckMask)

# Read the Planck Galactic mask
pathPlanckMask = pathIn+"/HFI_Mask_GalPlane-apo0_2048_R2.00.fits"
# field, fsky:
#0, 0.20
#1, 0.40
#2, 0.60
#3, 0.70
#4, 0.80
#5, 0.90
#6, 0.97
#7, 0.99
planckMask = hp.read_map(pathPlanckMask, field=3)

# upgrade to Nside=4096
planckMask = hp.ud_grade(planckMask, nSide)

# plot the Planck mask in Galactic coordinates
fig=plt.figure(0)
hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"planck_mask_G.pdf")
fig.clf()

# rotate it to equatorial coordinates
rot = hp.Rotator(coord=['G','C'])  # Transforms galactic to equatorial coordinates
planckMask = rot.rotate_map(planckMask)

fig=plt.figure(0)
hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"planck_mask_C.pdf")
fig.clf()

# rethreshold the mask, to make it zeros and ones
planckMask = (planckMask>0.999).astype(np.float)

fig=plt.figure(0)
hp.mollview(planckMask, fig=0, title="Planck Galactic mask", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"planck_mask_C_thresh.pdf")
fig.clf()

# save the Planck mask, rotated
hp.write_map(pathIn+"Planck_mask_4096_equatorial_coord.fits", planckMask, overwrite=True)

# read the Planck mask
planckMask = hp.read_map(pathIn+"Planck_mask_4096_equatorial_coord.fits")


#########################################################################
# Combine the masks into one, and extend it around the edges

fullMask = footMask * planckMask

print "fsky =", np.sum(fullMask) / len(fullMask)

fig=plt.figure(0)
hp.mollview(fullMask, fig=0, title="Full mask", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"f150_mask_foot_planck.pdf")
fig.clf()

# Gaussian smooth the mask
sigma = 4. * np.pi/180. # sigma in rad. We will truncate at around 2 sigmas
fullMask = hp.smoothing(fullMask, sigma=sigma)

fig=plt.figure(0)
hp.mollview(fullMask, fig=0, title="Footprint", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"f150_mask_foot_planck_smoothed.pdf")
fig.clf()

# threshold the mask
fullMask = (fullMask>0.95).astype(np.float)

# combine it with the point source mask
fullMask *= hPSMask

fig=plt.figure(0)
hp.mollview(fullMask, fig=0, title="Footprint", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"f150_mask_foot_planck_thresh.pdf")
fig.clf()

# save the full mask
hp.write_map(pathIn+"f150_mask_foot_planck_ps.fits", fullMask, overwrite=True)

# read the full mask
fullMask = hp.read_map(pathIn+"f150_mask_foot_planck_ps.fits")

fSky = np.sum(fullMask) / len(fullMask)
print "unmasked fsky =", fSky


# convert to CAR
carFullMask = reproject.enmap_from_healpix(fullMask, hitMap.shape, hitMap.wcs, rot=None)
# re-threshold
carFullMask = (carFullMask>0.95).astype(np.float)
# save
enmap.write_map(pathIn+"f150_mask_foot_planck_ps_car.fits", carFullMask)
# check that fsky is the same
checkFSky = np.sum(carFullMask) / (carFullMask.shape[1] * carFullMask.shape[2])
print "checking unmasked fSky =", checkFSky


#########################################################################
# Plot the masked hit count map

fig=plt.figure(0)
logHit = np.log10(np.abs(fullMask * hHitMap)+1.e-5)
hp.mollview(logHit, fig=0, title="Hit count", min=np.min(logHit), max=np.max(logHit), coord=None, cbar=False, unit='')
fig.savefig(pathFig+"f150_masked_hit_log.pdf")
fig.clf()

# save this map to file, for later illustrations
hp.write_map(pathIn+"f150_masked_hit_log.fits", logHit, overwrite=True)


#########################################################################
# Mask the temperature map

# mask T, Q, U
hMap[0] *= fullMask
hMap[1] *= fullMask
hMap[2] *= fullMask


#########################################################################
# Histograms: Healpy maps

# T histogram
x = hMap[0][fullMask.astype('bool')]
#x -= np.mean(x)
myHistogram(x, nBins=71, lim=(-10.*np.std(x), 10.*np.std(x)), sigma2Theory=110.**2, path=pathFig+"hist_T_masked.pdf", nameLatex=r'$T$ [$\mu$K]', semilogx=False, semilogy=True, doGauss=True)

# Q histogram
x = hMap[1][fullMask.astype('bool')]
#x -= np.mean(x)
myHistogram(x, nBins=71, lim=(-10.*np.std(x), 10.*np.std(x)), sigma2Theory=None, path=pathFig+"hist_Q_masked.pdf", nameLatex=r'$Q$ [$\mu$K]', semilogx=False, semilogy=True, doGauss=True)

# U histogram
x = hMap[2][fullMask.astype('bool')]
#x -= np.mean(x)
myHistogram(x, nBins=71, lim=(-10.*np.std(x), 10.*np.std(x)), sigma2Theory=None, path=pathFig+"hist_U_masked.pdf", nameLatex=r'$U$ [$\mu$K]', semilogx=False, semilogy=True, doGauss=True)


#########################################################################
# Histograms: CAR map, without pre-filtering

# T histogram
x = baseMap[0][carFullMask[0].astype('bool')]
myHistogram(x, nBins=71, lim=(-10.*np.std(x), 10.*np.std(x)), sigma2Theory=110.**2, path=pathFig+"hist_T_masked_car.pdf", nameLatex=r'$T$ [$\mu$K]', semilogx=False, semilogy=True, doGauss=True)

#########################################################################
# Histograms: CAR map, after pre-filtering
# --> prefiltering affects the pixel histograms

# do the pre-filtering
prefiltMap = utils.interpol_prefilter(baseMap[0], inplace=False)
prefiltCarFullMask = utils.interpol_prefilter(carFullMask[0], inplace=False)

# T histogram
x = prefiltMap[prefiltCarFullMask.astype('bool')]
myHistogram(x, nBins=71, lim=(-10.*np.std(x), 10.*np.std(x)), sigma2Theory=110.**2, path=pathFig+"hist_T_masked_car_prefiltered_cubic.pdf", nameLatex=r'$T$ [$\mu$K]', semilogx=False, semilogy=True, doGauss=True)



# do the pre-filtering
prefiltMap = utils.interpol_prefilter(baseMap[0], inplace=False, order=2)
prefiltCarFullMask = utils.interpol_prefilter(carFullMask[0], inplace=False, order=2)

# T histogram
x = prefiltMap[prefiltCarFullMask.astype('bool')]
myHistogram(x, nBins=71, lim=(-10.*np.std(x), 10.*np.std(x)), sigma2Theory=110.**2, path=pathFig+"hist_T_masked_car_prefiltered_quadratic.pdf", nameLatex=r'$T$ [$\mu$K]', semilogx=False, semilogy=True, doGauss=True)


#########################################################################

# plot masked T
mean = np.mean(hMap[0]) / fSky
sigma = np.std(hMap[0]) / fSky
#
fig=plt.figure(0)
hp.mollview(hMap[0], fig=0, min=mean-3.*sigma, max=mean+3.*sigma, title="T", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"f150_T_masked.pdf")
fig.clf()

# plot masked Q
mean = np.mean(hMap[1]) / fSky
sigma = np.std(hMap[1]) / fSky
#
fig=plt.figure(0)
hp.mollview(hMap[1], fig=0, min=mean-3.*sigma, max=mean+3.*sigma, title="Q", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"f150_Q_masked.pdf")
fig.clf()

# plot masked U
mean = np.mean(hMap[2]) / fSky
sigma = np.std(hMap[2]) / fSky
#
fig=plt.figure(0)
hp.mollview(hMap[2], fig=0, min=mean-3.*sigma, max=mean+3.*sigma, title="U", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"f150_U_masked.pdf")
fig.clf()


#########################################################################
# Power spectrum wrapper for healpy


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
         fig.clf()
      else:
         plt.show()

   return lCen, Cl, sCl


#########################################################################
# test my power spectrum code: generate a GRF, mask it, then take its power spectrum

# generate GRF
cl = np.array(map(cmb1_4.flensedTT, np.arange(3*nSide-1)))
cl = np.nan_to_num(cl)
testMap = hp.synfast(cl, nside=nSide)

# mask the map
testMap *= planckMask

# measure power spectrum, correcting for fsky
lCen, Cl, sCl = powerSpectrum(testMap, mask=planckMask, theory=[cmb1_4.flensedTT], fsCl=None, nBins=101, lRange=None, plot=True, path=pathFig+"check_synfast_mask_anafast.pdf", save=True)


#########################################################################
# Measure the power spectrum of the masked map

# CMB power spectra
flensedTT = lambda l: cmb1_4.flensedTT(l) * cmb1_4.fbeam(l)**2
fdetectorNoise = lambda l: cmb1_4.fdetectorNoise(l) * cmb1_4.fbeam(l)**2
ftotal = lambda l: cmb1_4.ftotal(l) * cmb1_4.fbeam(l)**2

# Plot the power spectrum
lCen, Cl, sCl = powerSpectrum(hMap[0], mask=fullMask, theory=[ftotal], fsCl=None, nBins=101, lRange=None, plot=True, path=pathFig+"f150_power_T_masked.pdf", save=True)

# Save it to file
data = np.zeros((len(lCen), 3))
data[:,0] = lCen
data[:,1] = np.nan_to_num(Cl)
data[:,2] = np.nan_to_num(sCl)
np.savetxt(pathIn+"f150_power_T_masked.txt", data)


#########################################################################
# Planck 143GHz map: measure its power spectrum

# load the Planck 143GHz
hP143Map = hp.read_map(pathIn+"HFI_SkyMap_143_2048_R2.02_full.fits")

# rotate the map to equatorial coordinates
rot = hp.Rotator(coord=['G','C'])  # Transforms galactic to equatorial coordinates
hP143Map = rot.rotate_map(hP143Map)

# convert it to muK^2
hP143Map *= 1.e6

# convert it to the right nside
hP143Map = hp.ud_grade(hP143Map, nSide)

# make it a fake T, Q, U map
#hMap = [hMap, hMap, hMap]

# CMB power spectra
flensedTT = lambda l: cmb7_3.flensedTT(l) * cmb7_3.fbeam(l)**2
fdetectorNoise = lambda l: cmb7_3.fdetectorNoise(l) * cmb7_3.fbeam(l)**2
ftotal = lambda l: cmb7_3.ftotal(l) * cmb7_3.fbeam(l)**2

# Measure its power spectrum
lCenP143, ClP143, sClP143 = powerSpectrum(hP143Map, mask=fullMask, theory=[ftotal], fsCl=None, nBins=101, lRange=None, plot=True, path=pathFig+"power_P143_masked.pdf", save=True)




#########################################################################
# Compare ACT+Planck and Planck, after converting to the same beam (1.4')

fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# theory
ClTheory = np.array(map(cmb1_4.flensedTT, lCen))
ax.loglog(lCen, lCen**2 * ClTheory, 'k-', label=r'CMB, $1.4^\prime$ beam')
ClTheory = np.array(map(cmb1_4.ftotalTT, lCen))
ax.loglog(lCen, lCen**2 * ClTheory, 'k--', label=r'+noise')
#
# Planck: convert from 5' beam to 1.4' beam
f = lambda l: (cmb1_4.fbeam(l) / cmb7_3.fbeam(l))**2
reBeam = np.array(map(f, lCen))
ax.errorbar(lCenP143, lCenP143**2 * ClP143 * reBeam, yerr=lCenP143**2 * sClP143 * reBeam, c='r', fmt='.', label=r'Planck 143, conv. to $1.4^\prime$ beam')
#
# ACT + Planck
ax.errorbar(lCen, lCen**2 * Cl, yerr=lCen**2 * sCl, c='b', fmt='.', label=r'ACT+Planck, $1.4^\prime$ beam')
#
ax.set_ylim((1.e2, 1.e5))
ax.legend(loc=3)
#
fig.savefig(pathFig+"f150_power_ACTPlanck_vs_P143.pdf")
fig.clf()


