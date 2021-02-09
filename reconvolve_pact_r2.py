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


pathFig = "./figures/cmb_map/"


##########################################################################
# Fourier beams
print("Reading the Fourier beams")

# the coadds have complicated beams
pathBeam150 = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "beam_f150_daynight.txt"
data = np.loadtxt(pathBeam150)
ell = data[:,0]
beam150F = data[:,1]
# regularization, as in Naess+20
v_cut = 1.e-2
ell_cut = np.argmin((np.abs(beam150F)-v_cut)**2)
beam150F = (beam150F>=v_cut) * beam150F + (beam150F<v_cut) * v_cut * beam150F[0] * (ell/ell_cut)**(2. * np.log(v_cut))
beam150F[0] = 1.
fBeam150F = interp1d(ell, beam150F, kind='linear', bounds_error=False, fill_value=(beam150F[0], beam150F[-1]))

pathBeam90 = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "beam_f090_daynight.txt"
data = np.loadtxt(pathBeam90)
ell = data[:,0]
beam90F = data[:,1]
# regularization, as in Naess+20
v_cut = 1.e-2
ell_cut = np.argmin((np.abs(beam90F)-v_cut)**2)
beam90F = (beam90F>=v_cut) * beam90F + (beam90F<v_cut) * v_cut * beam90F[0] * (ell/ell_cut)**(2. * np.log(v_cut))
beam90F[0] = 1.
fBeam90F = interp1d(ell, beam90F, kind='linear', bounds_error=False, fill_value=(beam90F[0], beam90F[-1]))

# night-only data
pathBeam150Night = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "beam_f150_night.txt"
data = np.loadtxt(pathBeam150)
ell = data[:,0]
beam150NightF = data[:,1]
# regularization, as in Naess+20
v_cut = 1.e-2
ell_cut = np.argmin((np.abs(beam150NightF)-v_cut)**2)
beam150NightF = (beam150F>=v_cut) * beam150F + (beam150F<v_cut) * v_cut * beam150F[0] * (ell/ell_cut)**(2. * np.log(v_cut))
beam150NightF[0] = 1.
fBeam150NightF = interp1d(ell, beam150NightF, kind='linear', bounds_error=False, fill_value=(beam150F[0], beam150F[-1]))



# TileC maps have Gussian beams
s = 1.6 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
beamTilecF = np.exp(-0.5 * s**2 * ell**2)
fBeamTilecF = interp1d(ell, beamTilecF, kind='linear', bounds_error=False, fill_value=(beamTilecF[0], beamTilecF[-1]))
#
s = 2.4 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
beamTilecDeprojF = np.exp(-0.5 * s**2 * ell**2)
fBeamTilecDeprojF = interp1d(ell, beamTilecDeprojF, kind='linear', bounds_error=False, fill_value=(beamTilecDeprojF[0], beamTilecDeprojF[-1]))



# Plot Fourier space beams
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
L = np.linspace(0., 5.e5, 10001)
#
ax.plot(L, fBeam150F(L), 'b-', label=r'150 GHz DR5')
ax.plot(L, -fBeam150F(L), 'b--')
#
#ax.plot(L, fBeam150NightF(L), 'y:', label=r'150 GHz DR5 night-only')
#ax.plot(L, -fBeam150NightF(L), 'y:')
#
ax.plot(L, fBeam90F(L), 'r-', label=r'98 GHz DR5')
ax.plot(L, -fBeam90F(L), 'r--')
#
ax.plot(L, fBeamTilecF(L), 'g-', label=r'ILC DR4')
ax.plot(L, fBeamTilecDeprojF(L), 'c-', label=r'ILC deproj DR4')
#
ax.legend(loc=3)
ax.set_xscale('log', nonposx='clip')
ax.set_yscale('log', nonposy='clip')
#
fig.savefig(pathFig + "beams_fourier.pdf", bbox_inches='tight')
plt.show()
fig.clf()






# Plot Fourier deconvolution kernels
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
L = np.linspace(0., 5.e5, 10001)
#
ax.plot(L, fBeam90F(L) / fBeam150F(L), 'b-', label=r'150 GHz to 98 GHz')
ax.plot(L, -fBeam90F(L) / fBeam150F(L), 'b--')
#
ax.plot(L, fBeamTilecF(L) / fBeam150F(L), 'g-', label=r'150 GHz to ILC')
ax.plot(L, -fBeamTilecF(L) / fBeam150F(L), 'g--')
#
ax.plot(L, fBeamTilecDeprojF(L) / fBeam150F(L), 'c-', label=r'150 GHz to ILC deproj')
ax.plot(L, -fBeamTilecDeprojF(L) / fBeam150F(L), 'c--')
#
ax.legend(loc=3)
ax.set_xscale('log', nonposx='clip')
ax.set_yscale('log', nonposy='clip')
#
fig.savefig(pathFig + "beams_deconvolution_kernels.pdf", bbox_inches='tight')
#plt.show()
fig.clf()






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
ax.plot(theta*180.*60./np.pi, beam150, 'b-', label=r'150 GHz DR5')
s = 1.3 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
y = np.exp(-0.5 * theta**2 / s**2) #/ (2. * np.pi * s**2)
ax.plot(theta*180.*60./np.pi, y, 'b--', label=r'G fwhm=$1.3^\prime$')
#
# Coadd at 90
ax.plot(theta*180.*60./np.pi, beam90, 'r-', label=r'98 GHz DR5')
s = 2.1 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
y = np.exp(-0.5 * theta**2 / s**2) #/ (2. * np.pi * s**2)
ax.plot(theta*180.*60./np.pi, y, 'r--', label=r'G fwhm=$2.1^\prime$')
#
# TileC no deproj
s = 1.6 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
y = np.exp(-0.5 * theta**2 / s**2) #/ (2. * np.pi * s**2)
ax.plot(theta*180.*60./np.pi, y, 'g-', label=r'ILC DR4: G fwhm=$1.6^\prime$')
#
# TileC  deproj
s = 2.4 * np.pi / (180. * 60.) / np.sqrt(8.*np.log(2.))
y = np.exp(-0.5 * theta**2 / s**2) #/ (2. * np.pi * s**2)
ax.plot(theta*180.*60./np.pi, y, 'c-', label=r'ILC deproj DR4: G fwhm=$2.4^\prime$')
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
##########################################################################
print("Read PACT 150 map")

pathMap = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits"

# read maps
iMap = enmap.read_map(pathMap)
# do the reconvolution in Fourier space
mapF = enmap.fft(iMap)
lMap = np.sqrt(np.sum(iMap.lmap()**2,0))


##########################################################################
# Reconvolve the map

print("Reconvolve 150 to 90")
pathOutMap = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvto90.fits"
oMap = enmap.ifft(mapF * fBeam90F(lMap) / fBeam150F(lMap)).real
enmap.write_map(pathOutMap, oMap)
print("check finite sum "+str(np.sum(oMap)))


print("Reconvolve 150 to TileC deproj")
pathOutMap = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvtotilecdeproj.fits"
oMap = enmap.ifft(mapF * fBeamTilecDeprojF(lMap) / fBeam150F(lMap)).real
enmap.write_map(pathOutMap, oMap)
print("check finite sum "+str(np.sum(oMap)))


print("Reconvolve 150 to TileC")
pathOutMap = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvtotilec.fits"
oMap = enmap.ifft(mapF * fBeamTilecF(lMap) / fBeam150F(lMap)).real
enmap.write_map(pathOutMap, oMap)
print("check finite sum "+str(np.sum(oMap)))



##########################################################################
##########################################################################
print("Read PACT 90 map")

pathMap = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f090_daynight_map.fits"

# read maps
iMap = enmap.read_map(pathMap)
# do the reconvolution in Fourier space
mapF = enmap.fft(iMap)
lMap = np.sqrt(np.sum(iMap.lmap()**2,0))


##########################################################################
# Reconvolve the map

print("Reconvolve 90 to TileC deproj")
pathOutMap = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f090_daynight_map_reconvtotilecdeproj.fits"
oMap = enmap.ifft(mapF * fBeamTilecDeprojF(lMap) / fBeam90F(lMap)).real
enmap.write_map(pathOutMap, oMap)
print("check finite sum "+str(np.sum(oMap)))


