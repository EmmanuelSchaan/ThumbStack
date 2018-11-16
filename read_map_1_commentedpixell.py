import numpy as np, time
import matplotlib.pyplot as plt

import flat_map
reload(flat_map)
from flat_map import *

import cmb
reload(cmb)
from cmb import *

#from enlib import enmap, utils, powspec
from pixell import enmap, utils, powspec, enplot

import healpy as hp

# copy rotfuncs.py somewhere on your python path,
# so you can import it
import rotfuncs


# for cori
plt.switch_backend('Agg')


#########################################################################
#########################################################################

# path for figures
pathFig = "./figures/cmb_map/planck_act_coadd_2018_08_10/"

# path for output
pathOut = "./output/cmb_map/planck_act_coadd_2018_08_10/"


#########################################################################
# CMB power spectra

cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=True)
cmb7_3 = StageIVCMB(beam=7.3, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=True)



#########################################################################
# load maps

# path for maps
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/"
pathMap = pathIn + "f150_daynight_all_map_mono.fits"
pathHit = pathIn + "f150_daynight_all_div_mono.fits"

print "Read map"
baseMap = enmap.read_map(pathMap)
hitMap = enmap.read_map(pathHit)

print("Map properties:")
print("nTQU, nY, nX = "+str(baseMap.shape))
print("WCS attributes: "+str(baseMap.wcs))


#########################################################################
#########################################################################
#########################################################################
# convert map to healpix, to plot and check the power spectrum

print("Convert from CAR to healpix")
hMap = enmap.to_healpix(baseMap)
hHitMap = enmap.to_healpix(hitMap)

# save the hit count map to file, for plotting purposes
#hp.write_map(pathIn+"healpix_f150_daynight_all_div_mono.fits", hHitMap, overwrite=True)

nSide = hp.get_nside(hMap)
print("Nside = "+str(nSide))



#########################################################################
#!!!!!!!! Temporarily, while debugging
'''
#!!!!!!!! load the Planck 143GHz instead of the Planck+ACT map
hMap = hp.read_map(pathIn+"HFI_SkyMap_143_2048_R2.02_full.fits")

# rotate the map to equatorial coordinates
rot = hp.Rotator(coord=['G','C'])  # Transforms galactic to equatorial coordinates
hMap = rot.rotate_map(hMap)

# convert it into muK^2
hMap *= 1.e6

# convert it to the right nside
hMap = hp.ud_grade(hMap, nSide)

# make it a fake T, Q, U map
hMap = [hMap, hMap, hMap]
'''

#########################################################################

'''
#!!!!!!!! Temporarily, while debugging: downgrade the map for speed
hMap = hp.ud_grade(hMap, 1024)
hHitMap = hp.ud_grade(hHitMap, 1024)
nSide = hp.get_nside(hMap)
'''

#########################################################################
# Plot the hit count map


fig=plt.figure(0)
hp.mollview(hHitMap, fig=0, title="Hit count", max=0.5*np.max(hHitMap), coord=None, cbar=True, unit='')
fig.savefig(pathFig+"hit_0.5max.pdf")
fig.clf()

fig=plt.figure(0)
logHit = np.log10(np.abs(hHitMap)+1.e-5)
hp.mollview(logHit, fig=0, title="Hit count", min=np.min(logHit), max=np.max(logHit), coord=None, cbar=False, unit='')
fig.savefig(pathFig+"hit_log.pdf")
fig.clf()


#########################################################################
# Generate a footprint mask from the hit map


# copy the masked temperature map
footMask = hHitMap.copy()

# generate a mask from all the null pixels,
footMask = (footMask<>0.).astype(np.float)

fig=plt.figure(0)
hp.mollview(footMask, fig=0, title="Footprint", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"footmask.pdf")
fig.clf()

# save the footprint to file
hp.write_map(pathIn+"footprint_mask.fits", footMask, overwrite=True)


# read the footprint mask
footMask = hp.read_map(pathIn+"footprint_mask.fits")


#########################################################################
# Mask the Milky Way with a Planck mask


# Read the Planck mask
#pathPlanckLensingMask = pathIn+"/COM_Mask_Lensing_2048_R2.00.fits"
#planckMask = hp.read_map(pathPlanckMask)

pathPlanckMask = pathIn+"/HFI_Mask_GalPlane-apo0_2048_R2.00.fits"
# fields, unmasked fsky:
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
fig.savefig(pathFig+"mask_foot_planck.pdf")
fig.clf()


# Gaussian smooth the mask
sigma = 4. * np.pi/180. # sigma in rad. We will truncate at around 2 sigmas
fullMask = hp.smoothing(fullMask, sigma=sigma)

fig=plt.figure(0)
hp.mollview(fullMask, fig=0, title="Footprint", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"mask_foot_planck_smoothed.pdf")
fig.clf()

# threshold the mask
fullMask = (fullMask>0.95).astype(np.float)

fig=plt.figure(0)
hp.mollview(fullMask, fig=0, title="Footprint", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"mask_foot_planck_thresh.pdf")
fig.clf()

# save the full mask
hp.write_map(pathIn+"mask_foot_planck.fits", fullMask, overwrite=True)


# read the full mask
fullMask = hp.read_map(pathIn+"mask_foot_planck.fits")


fSky = np.sum(fullMask) / len(fullMask)
print "unmasked fsky =", fSky


#########################################################################
# Plot the masked hit count map

fig=plt.figure(0)
logHit = np.log10(np.abs(fullMask * hHitMap)+1.e-5)
hp.mollview(logHit, fig=0, title="Hit count", min=np.min(logHit), max=np.max(logHit), coord=None, cbar=False, unit='')
fig.savefig(pathFig+"masked_hit_log.pdf")
fig.clf()

# save this map to file, for later illustrations
hp.write_map(pathIn+"masked_hit_log.fits", logHit, overwrite=True)


#########################################################################
# Mask the temperature map

# mask T, Q, U
hMap[0] *= fullMask
hMap[1] *= fullMask
hMap[2] *= fullMask

# plot masked T
mean = np.mean(hMap[0]) / fSky
sigma = np.std(hMap[0]) / fSky
#
fig=plt.figure(0)
hp.mollview(hMap[0], fig=0, min=mean-3.*sigma, max=mean+3.*sigma, title="T", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"T_masked.pdf")
fig.clf()

# plot masked Q
mean = np.mean(hMap[1]) / fSky
sigma = np.std(hMap[1]) / fSky
#
fig=plt.figure(0)
hp.mollview(hMap[1], fig=0, min=mean-3.*sigma, max=mean+3.*sigma, title="Q", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"Q_masked.pdf")
fig.clf()

# plot masked U
mean = np.mean(hMap[2]) / fSky
sigma = np.std(hMap[2]) / fSky
#
fig=plt.figure(0)
hp.mollview(hMap[2], fig=0, min=mean-3.*sigma, max=mean+3.*sigma, title="U", coord=None, cbar=True, unit='')
fig.savefig(pathFig+"U_masked.pdf")
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
#         ax.axhline(0.)
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
lCen, Cl, sCl = powerSpectrum(hMap[0], mask=fullMask, theory=[ftotal], fsCl=None, nBins=101, lRange=None, plot=True, path=pathFig+"power_T_masked.pdf", save=True)


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
# Compare ACT+Planck and Planck, after converting to the same beam

fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# theory
ClTheory = np.array(map(cmb1_4.flensedTT, lCen))
ax.loglog(lCen, lCen**2 * ClTheory, 'k-', label=r'CMB, $1.4^\prime$ beam')
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
fig.savefig(pathFig+"power_ACTPlanck_vs_P143.pdf")
fig.clf()










































#########################################################################
#########################################################################
# Do I understand the map? i.e. mean, std dev, units, plot?




'''
# get the map of coordinate values
pos = baseMap.posmap()
dec = pos[0]
ra = pos[1]
'''

'''
def quickPowerEnlib(map1, map2 = None, window = None, binningParams = (200, 10000, 40), return2dPower = False):
    """ Code snippet from Alex van Engelen
    """

    import orphics.maps, orphics.stats
    import numpy as np
    map1Local = map1.copy()
    

    if map2 is  None:
        map2Local = map1.copy()
    else:
        map2Local = map2.copy()
    

    if window is not  None:
        map1Local *= window / np.sqrt(np.mean(window**2))
        map2Local *= window / np.sqrt(np.mean(window**2))
    

    fc = orphics.maps.FourierCalc(map1Local.shape,map1Local.wcs)
    power2d, kmap1, kmap2 =  fc.power2d(np.nan_to_num(map1Local),np.nan_to_num( map2Local))


    lbin_edges = np.arange(binningParams[0], binningParams[1], binningParams[2])
    lbinner = orphics.stats.bin2D( power2d.modlmap(), lbin_edges)


    binnedPower =    lbinner.bin(power2d)


    return {'lbin':  binnedPower[0], \
            'clbin' : binnedPower[1], \
            'dlbin': binnedPower[1] * binnedPower[0] * (binnedPower[0] + 1) / (2. * np.pi) ,
            'power2d' : power2d if return2dPower else None}
'''

#########################################################################
# Plot the map using enplot
#!!! Not working, probably because of my PIL version

'''
enplot.plot(baseMap, oname="./figures/tests/planck_act_coadd_2018_08_10_f150.pdf", quantile=0.16)
# range = (0.1, 0.4, 0.2)
# min = 0
# max = 1
'''

'''
plots = enplot.get_plots(baseMap)
enplot.write("./figures/tests/planck_act_coadd_2018_08_10_f150.png",plots)
'''

#########################################################################
# Naive imshow plots of the maps
'''
# temperature
fig=plt.figure(0)
ax=fig.add_subplot(111)
im=ax.imshow(baseMap[0], vmin=-1.*np.std(baseMap[0].flatten()), vmax=1.*np.std(baseMap[0].flatten()))   # T
fig.savefig("./figures/tests/imshow_T.pdf")
fig.clf()

# Q polarization
fig=plt.figure(0)
ax=fig.add_subplot(111)
im=ax.imshow(baseMap[1], vmin=-1.*np.std(baseMap[1].flatten()), vmax=1.*np.std(baseMap[1].flatten()))   # T
fig.savefig("./figures/tests/imshow_Q.pdf")
fig.clf()

# U polarization
fig=plt.figure(0)
ax=fig.add_subplot(111)
im=ax.imshow(baseMap[2], vmin=-1.*np.std(baseMap[2].flatten()), vmax=1.*np.std(baseMap[2].flatten()))   # T
fig.savefig("./figures/tests/imshow_U.pdf")
fig.clf()
'''


#########################################################################
# use my map class to get nice plots
'''
nPol, nX, nY = np.shape(baseMap)
sizeX = np.max(dec) - np.min(dec)
sizeY = np.max(ra) - np.min(ra)
print "map size=", sizeX*180./np.pi, "deg *", sizeX*180./np.pi, "deg"
print "pixel numbers =", nX, "*", nY

myMap = FlatMap(nX=nX, nY=nY, sizeX=sizeX, sizeY=sizeY)



print "Plot T, Q, U maps"
for i in range(3):
   path = "./figures/tests/act_planck_coadd_"+str(i)+".pdf"
   myMap.plot(baseMap[i], save=True, path=path)
#   fig=plt.figure(0)
#   ax=fig.add_subplot(111)
#   #
#   ax.imshow(baseMap[i])
#   #
#   fig.savefig("./figures/tests/act_planck_coadd_"+str(i)+".pdf")
#   fig.clf()
   print "- done with i="+str(i+1)+" of 3"
'''


#########################################################################
# Play with the powspec class
'''
# generate an empty map
deg = utils.degree
# Build an empty map from scratch (usually you don't
# have to do this, as you'll be reading in an existing map)
shape, wcs = enmap.geometry(pos=[[-5*deg,-5*deg],[5*deg,5*deg]], res=0.5*utils.arcmin, proj="car")
map = enmap.zeros(shape, wcs)
print "Made a map with shape %s and wcs %s" % (str(shape),str(wcs))

# We didn't get Q and U because we didn't specify a 3-component
# shape. Let's add this to the shape:
shape = (3,)+shape

# Let's try to make a map with a CMB power spectrum. We first
# read in a spectrum, and then make a random realization using it.
ps = powspec.read_spectrum("cl_lensed.dat")
map = enmap.rand_map(shape, wcs, ps)

print "Displaying full T,Q,U map"
for m in map: plt.matshow(m)
show()

# Let's apply a lowpass filter. We could do that
# simply by using enmap.smooth_gauss, but we will
# do it manually to show how it's done. The lowpass
# filter is defined in fourier space. map2harm takes
# you from real-space TQU to fourier-space TEB.
# It's the flat-sky equivalent of map2alm.
# The result is ({t,e,b},nly,nlx).
fmap = enmap.map2harm(map)
# Get the 2d fourier mode (ly,lx) corresponding to
# each fourier-space pixel
lmap = fmap.lmap()
# I only want the 1d multipole (this is just pythagoras)
l    = np.sum(lmap**2,0)**0.5
# Apply a gaussian filter with standard deviation of l=1000
fmap *= np.exp(-0.5*(l/1000)**2)
# And transform back
smooth = enmap.harm2map(fmap)
print "Displaying smoothed map"
for m in smooth: matshow(m)
show()

# We can also use this interpolation to zoom. Let's
# create a very high res coordinate system in a small
# area.
shape3, wcs3 = enmap.geometry(pos=[[-0.2*deg,-0.2*deg],[0.2*deg,0.2*deg]],res=0.05*utils.arcmin)
map_zoom = map.project(shape3, wcs3)
print "Displaying zoomed map"
matshow(map_zoom[0]); show()

# The interpolation uses 3rd order spline interpolation by default.
# If we want to see how the old pixels map to the new ones we can
# use 0th order interpolation (nearest neighbor). This is very fast,
# but not what you usually want to do.
map_zoom_nn = map.project(shape3, wcs3, order=0)
print "Displaying 0th order zoomed map"
matshow(map_zoom_nn[0]); show()

'''


#########################################################################
#########################################################################
# Extract postage stamp map at desired location
'''
# setup the interpolation algorithm,
# done once for all, to speed up subsequent calls
print "Set up interpolation"
tStart = time()
baseMap = utils.interpol_prefilter(baseMap, inplace=True)
tStop = time()
print "took", tStop-tStart, "sec"


# define geometry of small square maps to be extracted
# here 1deg * 1deg, with 0.25arcmin pixel
# the center of the small square map will be specified later
# car: not equal area, but equally-spaced coordinates
# cea: equal area pixels, but coordinates are not equally spaced
# previous kSZ paper: went out to 5arcmin
# Planck paper CV: went out to 18 arcmin
# probably good to extract bigger stamp than needed for now
shape, wcs = enmap.geometry(np.array([[-1,-1],[1,1]])*utils.degree, res=0.25*utils.arcmin, proj='cea')
stampMap = enmap.zeros(shape, wcs)

# coordinates of the square map (between -1 and 1 deg)
# output map position [{dec,ra},ny,nx]
opos = stampMap.posmap()

# coordinate of the center of the square map we want to extract
# (known point source)
# ra, dec in this order
sourcecoord = np.array([-8.14, -2.08])*utils.degree

# corresponding true coordinates on the big healpy map
ipos = rotfuncs.recenter(opos[::-1], [0,0,sourcecoord[0],sourcecoord[1]])[::-1]

# extract the small square map by interpolating the big healpy map
stampMap = baseMap.at(ipos, prefilter=False, mask_nan=False)

# save the map
enmap.write_map("./output/tests/stamp_map.fits", stampMap)


# take a quick look in temperature,
# to check that the source is there
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.imshow(stampMap[0])
#
fig.savefig("./figures/tests/stamp.pdf")
'''


#########################################################################

'''
# mpi4py: parallel loop over sources
# without load balancing
sources = np.loadtxt()
nsrc = len(sources)
for ind in range(comm.rank, nsrc, comm.size):
    source = sources[ind]

# useful functions from enlib
# the fft assumes that the pixels are equally spaced
enmap.fft
enmap.ifft
enmap.lmap
enmap.apod
l = np.sum(m.lmap()**2,0)**0.5
l = m.modlmap()
'''

