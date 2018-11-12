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
# Look at full map

# load a healpy map
path = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/"
pathMap = path + "f150_daynight_all_map_mono.fits"
pathHit = path + "f150_daynight_all_div_mono.fits"


print "Read map"
baseMap = enmap.read_map(pathMap)
hitMap = enmap.read_map(pathHit)

print("Map properties:")
print("nTQU, nY, nX = "+str(baseMap.shape))
print("WCS attributes: "+str(baseMap.wcs))

print("Map properties:")
print("nTQU, nY, nX = "+str(baseMap.shape))
print("WCS attributes: "+str(baseMap.wcs))

# setup the interpolation algorithm,
# done once for all, to speed up subsequent calls
print "Set up interpolation"
tStart = time()
baseMap = utils.interpol_prefilter(baseMap, inplace=True)
tStop = time()
print "took", tStop-tStart, "sec"

#########################################################################
# Plot the map using enplot

'''
enplot.plot(baseMap, oname="./figures/tests/planck_act_coadd_2018_08_10_f150.pdf", quantile=0.16)
# range = (0.1, 0.4, 0.2)
# min = 0
# max = 1
'''

#!!! Not working, probably because of my PIL version
plots = enplot.get_plots(baseMap)
enplot.write("./figures/tests/planck_act_coadd_2018_08_10_f150.png",plots)









#########################################################################
# Naive imshow plots of the maps

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


#########################################################################
#########################################################################
# convert map to healpix, to plot and check the power spectrum


print("Convert CAR to healpix, to plot")
hMap = enmap.to_healpix(baseMap)
hHitMap = enmap.to_healpix(hitMap)

# save the hit count map to file, for plotting purposes
hp.write_map(path+"healpix_f150_daynight_all_div_mono.fits", hHitMap)

nSide = hp.get_nside(hMap)
print("Nside = "+str(nSide))



# Hit map
fig=plt.figure(0)
#
#hp.mollview(hHitMap, fig=0, min=0., max=0.2*np.max(hHitMap), title="hit", coord=None, cbar=True, unit='')
hp.mollview(np.log(np.abs(hHitMap)+1.e-5), fig=0, title="", coord=None, cbar=False, unit='')
fig.savefig("./figures/tests/mollweide_hit_log.pdf")
fig.clf()


# T
mean = np.mean(hMap[0])
sigma = np.std(hMap[0])
#
fig=plt.figure(0)
hp.mollview(hMap[0], fig=0, min=mean-sigma, max=mean+sigma, title="T", coord=None, cbar=True, unit='')
fig.savefig("./figures/tests/mollweide_T.pdf")
fig.clf()

# Q
mean = np.mean(hMap[1])
sigma = np.std(hMap[1])
#
fig=plt.figure(0)
hp.mollview(hMap[1], fig=0, min=mean-sigma, max=mean+sigma, title="Q", coord=None, cbar=True, unit='')
fig.savefig("./figures/tests/mollweide_Q.pdf")
fig.clf()

# U
mean = np.mean(hMap[2])
sigma = np.std(hMap[2])
#
fig=plt.figure(0)
hp.mollview(hMap[2], fig=0, min=mean-sigma, max=mean+sigma, title="U", coord=None, cbar=True, unit='')
fig.savefig("./figures/tests/mollweide_U.pdf")
fig.clf()


## Other projections
#
## T
#mean = np.mean(hMap[0])
#sigma = np.std(hMap[0])
##
#fig=plt.figure(0)
#hp.cartview(hMap[0], fig=0, min=mean-sigma, max=mean+sigma, title="T", coord=None, cbar=True, unit='')
#fig.savefig("./figures/tests/cartview_T.pdf")
#fig.clf()
#
## T
#mean = np.mean(hMap[0])
#sigma = np.std(hMap[0])
##
#fig=plt.figure(0)
#hp.orthview(hMap[0], fig=0, min=mean-sigma, max=mean+sigma, title="T", coord=None, cbar=True, unit='')
#fig.savefig("./figures/tests/orthview_T.pdf")
#fig.clf()
#
## T
#mean = np.mean(hMap[0])
#sigma = np.std(hMap[0])
##
#fig=plt.figure(0)
#hp.gnomview(hMap[0], fig=0, min=mean-sigma, max=mean+sigma, title="T", coord=None, cbar=True, unit='')
#fig.savefig("./figures/tests/gnomview_T.pdf")
#fig.clf()
#
## T
#mean = np.mean(hMap[0])
#sigma = np.std(hMap[0])
##
#fig=plt.figure(0)
#hp.azeqview(hMap[0], fig=0, min=mean-sigma, max=mean+sigma, title="T", coord=None, cbar=True, unit='')
#fig.savefig("./figures/tests/azeqview_T.pdf")
#fig.clf()



#########################################################################
#########################################################################
# Do I understand the map? i.e. mean, std dev, units, plot?


# print map size in deg, nb of pixels, angular resolution in arcmin
# get units, rescale the color plot, get a power spectrum
print "T map has mean", np.mean(baseMap[0]), " st.d dev.", np.std(baseMap[0])

# convert map units from Jy to K?
cmb = CMB()
nu = 150.e9 # freq in Hz
#cmb.dBdT(nu, cmb.Tcmb)


# get the map of coordinate values
pos = baseMap.posmap()
dec = pos[0]
ra = pos[1]


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
# use healpix to get power spectrum and everything?
'''
# convert map to healpix, so that I can check the power spectrum
healpMap = enmap.to_healpix(baseMap)

## get power spectrum
#Cl = hp.anafast(masked_pl, lmax = 1500)
#pixwinf = hp.pixwin(N_side)[0:len(Cl)]
#Cl = Cl / (pixwinf **2)
#
## Number of bins and range
#Nbins = 8
#lmin = 30
#lmax = 1200
#
## binning
#bins = np.round(np.linspace(lmin, lmax, Nbins+1)) bins = bins.astype(int)
#lcenterbin = np.zeros(len(bins)-1)
#binnedCl = np.zeros(len(bins)-1)
#
## binning
#for k in range(0, len(bins)-1):
#lmaxvec = np.arange(bins[k], bins[k+1], 1) lcenterbin[k] = np.round(0.5 * (bins[k] + bins[k+1])) for l in lmaxvec:
#binnedCl[k] += Cl[l]
#binnedCl[k] = binnedCl[k] / len(lmaxvec)
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

