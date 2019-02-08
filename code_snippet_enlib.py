'''
import numpy as np, time
import matplotlib.pyplot as pl

# copy rotfuncs.py somewhere on your python path,
# so you can import it
import rotfuncs

# git clone https://github.com/amaurea/enlib
# then some synlink magic...
# it all happened in seconds when Sigurd did it :)
from enlib import enmap, utils
'''
from headers import *

# Planck + ACT 150GHz day and night
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/"
pathMap = pathIn + "f150_daynight_all_map_mono.fits"
pathHit = pathIn + "f150_daynight_all_div_mono.fits"
pathMask = pathIn + "f150_mask_foot_planck_ps_car.fits"

# load a healpy map
bigmap = enmap.read_map(pathHit)


#######################################################################
#######################################################################
# Get map value at sky location: interpolated

# setup the interpolation algorithm,
# done once for all, to speed up subsequent calls
bigmap = utils.interpol_prefilter(bigmap, inplace=True)

# specify sky coordinates
ra = 22. # deg
dec = 0. # deg

# interpolate the map to the given sky coordinates
sourcecoord = np.array([dec, ra])*utils.degree
bigmap.at(sourcecoord, prefilter=False, mask_nan=False)


#######################################################################
# Get map value at sky location: closest pixel

# specify sky coordinates
ra = 22. # deg
dec = -90. # deg

yPix,xPix = enmap.sky2pix(bigmap.shape, bigmap.wcs, np.array([dec,ra])*utils.degree)

# round to closest pixel
yPix = np.int(np.round(yPix))
xPix = np.int(np.round(xPix))

print "ra, dec =", ra, dec
print "x, y =", xPix, yPix
print "nX, nY =", bigmap.shape[1], bigmap.shape[0]

# evaluate the map at the desired pixel
bigmap[yPix, xPix]



#######################################################################
# Testing them on the hit count map, to see typical hit values


# check the hit map as a function of dec, at fixed ra = 0
# this will return an error if the ra,dec value is outside the map
def sky2map(dec, ra=45.):
   yPix,xPix = enmap.sky2pix(bigmap.shape, bigmap.wcs, np.array([dec,ra])*utils.degree)
   # round to closest pixel
   yPix = np.int(np.round(yPix))
   xPix = np.int(np.round(xPix))
   # evaluate the map at the desired pixel
   return bigmap[yPix, xPix]


# check the hit map as a function of dec, at fixed ra = 0
def sky2map(dec, ra=45.):
   # interpolate the map to the given sky coordinates
   sourcecoord = np.array([dec, ra])*utils.degree
   return bigmap.at(sourcecoord, prefilter=False, mask_nan=False)



Dec = np.linspace(-90., 90., 201)

Hit = np.array(map(sky2map, Dec))

fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.plot(Dec, Hit, ls='-')
ax.plot(Dec, -Hit, ls='--')
#
ax.set_yscale('log', nonposx='clip')
ax.set_xlabel(r'Dec')
ax.set_ylabel(r'Hit')
ax.set_xlim((-90., 90.))
ax.set_ylim((1.e-5, 1.e-2))
#
fig.savefig("./test.pdf", bbox_inches='tight')
fig.clf()




#It should be (dec,ra). Note that enmap.at() uses interpolation to get that value. If you want to get the value at the nearest pixel, I would do
#y,x = enmap.sky2pix(np.array([dec,ra])*utils.degree)
#then round it however you like, and then access it with imap[...,y,x]. These are all approximate of course. If you want something that has been written a bit more carefully, you should consider
#pixell.reproject.postage_stamp (https://github.com/simonsobs/pixell/blob/master/pixell/reproject.py#L13 )
#whose central pixel should have what you need. One still needs to worry about pixel window from interpolation.



#######################################################################
##########################################################################
# Extract a square cutout of the map

# setup the interpolation algorithm,
# done once for all, to speed up subsequent calls
bigmap = utils.interpol_prefilter(bigmap, inplace=True)

# define geometry of small square maps to be extracted
# here 1deg * 1deg, with 0.25arcmin pixel
# car: not equal area, but equally-spaced coordinates
# cea: equal area pixels, but coordinates are not equally spaced
shape, wcs = enmap.geometry(np.array([[-1,-1],[1,1]])*utils.degree, res=0.25*utils.arcmin, proj='cea')
sqmap = enmap.zeros(shape, wcs)

# coordinate of the center of the square map we want to extract
# (known point source)
# ra, dec in this order
sourcecoord = np.array([-8.14, -2.08])*utils.degree

# coordinates of the square map (between -1 and 1 deg)
# output map position [{dec,ra},ny,nx]
opos = sqmap.posmap()

# corresponding true coordinates on the big healpy map
ipos = rotfuncs.recenter(opos[::-1], [0,0,sourcecoord[0],sourcecoord[1]])[::-1]

# extract the small square map by interpolating the big healpy map
sqmap[:] = bigmap.at(ipos, prefilter=False, mask_nan=False)

# save the map
enmap.write_map("out_map.fits", sqmap)

# take a quick look,
# to check that the source is there
pl.imshow(sqmap[:])
pl.show()

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
