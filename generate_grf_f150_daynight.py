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

pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/"


#########################################################################

# theory curve for the CMB power spectrum
cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)

# Load the actual power spectrum of the Planck+ACT+ACTPol coadd
data = np.genfromtxt(pathIn+"f150_power_T_masked.txt")
lCen = data[:,0]
ClMeasured = data[:,1]
sClMeasured = data[:,2]

# Stitch the two spectra at ell=100, to avoid high noise, and interpolate
ClTheory = np.array(map(cmb1_4.flensedTT, lCen))
ClStitched = ClTheory * (lCen<=100.) + ClMeasured * (lCen>100.)
fClStitched = interp1d(lCen, ClStitched, kind='linear', bounds_error=False, fill_value=0.)

# generate power spectrum array for every ell
L = np.arange(1., np.max(lCen))
Cl = np.array(map(fClStitched, L))

# Generate pixell map, GRF with the desired power spectrum
nSide = 4096   # this is what pixell chooses when converting our CAR map to healpix
hpGrfMap = hp.sphtfunc.synfast(Cl, nSide, lmax=None, mmax=None, alm=False, pol=False, pixwin=False, fwhm=0.0, sigma=None, new=False, verbose=False)


# Check the power spectrum of the GRF matches the input



# read CAR hit map to get the desired pixellation properties
pathHit = pathIn + "f150_daynight_all_div_mono.fits"
hitMap = enmap.read_map(pathHit)

# Convert to pixell CAR map
grfMap = reproject.enmap_from_healpix(hpGrfMap, hitMap.shape, hitMap.wcs, rot=None)

# save the new map
enmap.write_map(pathIn+"grf_f150_daynight.fits", grfMap)








