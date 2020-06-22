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

pathMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits"

pathFig = "./figures/cmb_map/"


##########################################################################
# Fourier beams
print("Reading the Fourier beams")

# the coadds have complex beams
pathBeam150 = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28_r2/" + "beam_f150_daynight.txt"
data = np.loadtxt(pathBeam150)
ell = data[:,0]
beam150F = data[:,1]
fBeam150F = interp1d(ell, beam150F, kind='linear', bounds_error=False, fill_value=(beam150F[0], beam150F[-1]))
#
pathBeam90 = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28_r2/" + "beam_f090_daynight.txt"
data = np.loadtxt(pathBeam90)
ell = data[:,0]
beam90F = data[:,1]
fBeam90F = interp1d(ell, beam90F, kind='linear', bounds_error=False, fill_value=(beam90F[0], beam90F[-1]))

# TileC maps have Gussian beams
s = 1.6 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
beamTilecF = np.exp(-0.5 * s**2 * ell**2)
fBeamTilecF = interp1d(ell, beamTilecF, kind='linear', bounds_error=False, fill_value=(beamTilecF[0], beamTilecF[-1]))
#
s = 2.4 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
beamTilecDeprojF = np.exp(-0.5 * s**2 * ell**2)
fBeamTilecDeprojF = interp1d(ell, beamTilecDeprojF, kind='linear', bounds_error=False, fill_value=(beamTilecDeprojF[0], beamTilecDeprojF[-1]))



# Plot Fourier space beams
#fig=plt.figure(0)
#ax=fig.add_subplot(111)
##
#ax.plot(ell, beam150F, 'b-', label=r'150GHz')
#s = 1.3 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
#y = np.exp(-0.5 * s**2 * ell**2)
#ax.plot(ell, y, 'b--', label=r'Gaussian with fwhm=1.3')
##
#ax.plot(ell, beam90F, 'r-', label=r'90GHz')
#s = 2.1 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
#y = np.exp(-0.5 * s**2 * ell**2)
#ax.plot(ell, y, 'r--', label=r'Gaussian with fwhm=2.1')
##
#ax.legend(loc=1)
#ax.set_xscale('log', nonposx='clip')
#ax.set_yscale('log', nonposy='clip')
##
#plt.show()


##########################################################################
# Real-space beams
print("Checking the real space beams")

# compute real space beams
theta = np.linspace(0., 10., 101) * np.pi / 180. / 60. # [rad]
beam150 = utils.beam_transform_to_profile(beam150F, theta)
#beam150 /= np.trapz(2. * np.pi * theta * beam150, x=theta)
beam150 /= beam150[0]

beam90 = utils.beam_transform_to_profile(beam90F, theta)
#beam90 /= np.trapz(2. * np.pi * theta * beam90, x=theta)
beam90 /= beam90[0]


# plot real space beams
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# Coadd at 150
ax.plot(theta*180.*60./np.pi, beam150, 'b-', label=r'150GHz')
s = 1.3 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
y = np.exp(-0.5 * theta**2 / s**2) #/ (2. * np.pi * s**2)
ax.plot(theta*180.*60./np.pi, y, 'b--', label=r'G fwhm=$1.3^\prime$')
#
# Coadd at 90
ax.plot(theta*180.*60./np.pi, beam90, 'r-', label=r'90GHz')
s = 2.1 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
y = np.exp(-0.5 * theta**2 / s**2) #/ (2. * np.pi * s**2)
ax.plot(theta*180.*60./np.pi, y, 'r--', label=r'G fwhm=$2.1^\prime$')
#
# TileC no deproj
s = 1.6 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
y = np.exp(-0.5 * theta**2 / s**2) #/ (2. * np.pi * s**2)
ax.plot(theta*180.*60./np.pi, y, 'g-', label=r'TileC: G fwhm=$1.6^\prime$')
#
# TileC  deproj
s = 2.4 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
y = np.exp(-0.5 * theta**2 / s**2) #/ (2. * np.pi * s**2)
ax.plot(theta*180.*60./np.pi, y, 'c-', label=r'TileC deproj: G fwhm=$2.4^\prime$')
#
ax.legend(loc=3, fontsize='x-small', labelspacing=0.1, frameon=False)
ax.set_xscale('log', nonposx='clip')
ax.set_yscale('log', nonposy='clip')
ax.set_ylim((1.e-3, 2.))
ax.set_xlabel(r'$\theta$ [arcmin]')
ax.set_ylabel(r'$B\left(\theta \right) / B(0)$ ')
#
fig.savefig(pathFig + "beams.pdf", bbox_inches='tight')
#plt.show()
fig.clf()


##########################################################################
print("Reconvolve 150 to 90")

pathOutMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvto90.fits"


# read maps
iMap = enmap.read_map(pathMap)
# do the reconvolution in Fourier space
mapF = enmap.fft(iMap)
lMap = np.sum(iMap.lmap()**2,0)
mapF *= fBeam90F(lMap) / fBeam150F(lMap)
oMap = enmap.ifft(mapF).real
# save it to file
enmap.write_map(pathOutMap, oMap)


##########################################################################
print("Reconvolve 150 to TileC deproj")

pathOutMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvtotilecdeproj.fits"


# read maps
iMap = enmap.read_map(pathMap)
# do the reconvolution in Fourier space
mapF = enmap.fft(iMap)
lMap = np.sum(iMap.lmap()**2,0)
mapF *= fBeamTilecDeprojF(lMap) / fBeam150F(lMap)
oMap = enmap.ifft(mapF).real
# save it to file
enmap.write_map(pathOutMap, oMap)


##########################################################################
print("Reconvolve 150 to TileC")

pathOutMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvtotilec.fits"


# read maps
iMap = enmap.read_map(pathMap)
# do the reconvolution in Fourier space
mapF = enmap.fft(iMap)
lMap = np.sum(iMap.lmap()**2,0)
mapF *= fBeamTilecF(lMap) / fBeam150F(lMap)
oMap = enmap.ifft(mapF).real
# save it to file
enmap.write_map(pathOutMap, oMap)



