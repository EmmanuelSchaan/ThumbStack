import numpy as np, time
import matplotlib.pyplot as plt

import flat_map
reload(flat_map)
from flat_map import *

import cmb
reload(cmb)
from cmb import *

from enlib import enmap, utils, powspec
# copy rotfuncs.py somewhere on your python path,
# so you can import it
import rotfuncs


#########################################################################
# Look at full map




# load a healpy map
path = "/global/cscratch1/sd/eschaan/ksz_act_planck/data/maps/planck_act_coadd_2018_08_10/"
path += "f150_daynight_all_map_mono.fits"

print "Read map"
baseMap = enmap.read_map(path)

# setup the interpolation algorithm,
# done once for all, to speed up subsequent calls
print "Set up interpolation"
baseMap = utils.interpol_prefilter(baseMap, inplace=True)

#########################################################################
#########################################################################
# Do I understand the map? ie mean, std dev, units, plot?


# print map size in deg, nb of pixels, angular resolution in arcmin
# get units, rescale the color plot, get a power spectrum
print "T map has mean", np.mean(baseMap[0]), " st.d dev.", np.std(baseMap[0])

# convert map units from Jy to K?
cmb = CMB()
nu = 150.e9 # freq in Hz
cmb.dBdT(nu, self.Tcmb)


# get the map of coordinate values
pos = baseMap.posmap()
dec = pos[0]
ra = pos[1]

# use my map class to get nice plots
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

'''
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

