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
#plt.switcsh_backend('Agg')


#########################################################################
#########################################################################

# path for figures
pathFig = "./figures/cmb_map/planck_act_coadd_2018_08_10/f150_daynight/"

# path for output
pathOut = "./output/cmb_map/planck_act_coadd_2018_08_10/f150_daynight/"


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
pathFullMask = pathIn + "f150_mask_foot_planck_ps_car.fits"

print "Read maps"
baseMap = enmap.read_map(pathMap)
hitMap = enmap.read_map(pathHit)
fullMask = enmap.read_map(pathFullMask)

print("Map properties:")
print("nTQU, nY, nX = "+str(baseMap.shape))
print("WCS attributes: "+str(baseMap.wcs))


#########################################################################
#########################################################################
# Extract postage stamp map at desired location

# setup the interpolation algorithm,
# done once for all, to speed up subsequent calls
print "Set up interpolation"
tStart = time()
baseMap = utils.interpol_prefilter(baseMap, inplace=True)
tStop = time()
print "took", tStop-tStart, "sec"


# define geometry of small square maps to be extracted
# here 0.5deg * 0.5deg, with 0.25arcmin pixel
# the center of the small square map will be specified later
# car: not equal area, but equally-spaced coordinates
# cea: equal area pixels, but coordinates are not equally spaced
# previous kSZ paper: went out to 5arcmin
# Planck paper CV: went out to 18 arcmin
# --> let's go to 18'*sqrt(2) = 25' \simeq 0.5deg
# probably good to extract bigger stamp than needed for now
shape, wcs = enmap.geometry(np.array([[-0.5,-0.5],[0.5,0.5]])*utils.degree, res=0.25*utils.arcmin, proj='cea')
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
#enmap.write_map("./output/tests/stamp_map.fits", stampMap)


# take a quick look in temperature,
# to check that the source is there
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.imshow(stampMap[0])
#
fig.savefig("./figures/tests/stamp.pdf")



#########################################################################

# implement the disk - ring filter

# local coordinates in rad.
# zero is at the center of the map
dec = opos[0,:,:]
ra = opos[1,:,:]
radius = np.sqrt(ra**2 + dec**2)

theta0 = 1./60. * np.pi/180.
theta1 = theta0 * np.sqrt(2.)

# filter map
inDisk = 1.*(radius<=theta0)
inRing = 1.*(radius>theta0)*(radius<=theta1)
filter = inDisk / np.sum(inDisk)
filter -= inRing / np.sum(inRing)

# check that the filter integrates to zero
print "- disk-ring filter sums over pixels to "+str(np.sum(filter))
print "  (should be 0; to be compared with "+str(len(filter.flatten()))+")"

# output of the filter
filtOutput = np.sum(filter * stampMap)
print "- output of the matched filtering:"+str(filtOutput)+" muK"


# mean variance of the filter
varFilter = np.sum(filter**2 * 1./stampHit)
print "- mean var from weight map:"+str(varFilter)+" (arbitrary unit)"

# count nb of pixels where filter is strictly positive
nbPix = len(np.where(filter>0.)[0])
print "- nb of pixels where filter>0: "+str(nbPix)

# estimate area of strictly positive part of filter
pixArea = ra.area / len(ra.flatten()) # in sr
diskArea = np.sum(inDisk) * pixArea  # disk area in sr
print "  ie area of "+str(diskArea)+" sr"
#filtArea *= (self.U.ComovDist(1./(1.+z0), 1.) /(1.+z0) )**2
#print "  ie "+str(diskArea)+"(Mpc/h)^2 physical at z="+str(z0)

#print "- stellar mass: "+str(self.Catalog.Mstellar[iobj])+" Msun"
#print "- halo mass: "+str(self.Catalog.Mvir[iobj])+" Msun"
#print "- reconstructed velocity: "+str(self.Catalog.velR[iobj])+" km/s"

## just for testing purpose
#expectedTSZ = self.TSZ[iobj] / filtArea
#expectedKSZ = self.Tau[iobj] * (self.Catalog.velR[iobj]/3.e5) * 2.726e6 / filtArea
#print "- expected tSZ signal: "+str(expectedTSZ)+" muK"
#print "- expected kSZ signal: "+str(expectedKSZ)+" muK"
#print "- output of the matched filtering:"+str(filtOutput)+" muK"







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

