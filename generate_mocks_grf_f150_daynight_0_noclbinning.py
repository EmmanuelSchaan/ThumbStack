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

nProc = 10  # 32 cores on cori, but would run out of memory. 10 works. 15 runs out of memory.
nMocks = 10

pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/"
pathOut = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/output/cmb_map/mocks_grf_planck_act_coadd_2018_08_10/"
pathFig = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/cmb_map/mocks_grf_planck_act_coadd_2018_08_10/"


#########################################################################

# theory curve for the CMB power spectrum
cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)

# Load the actual power spectrum of the Planck+ACT+ACTPol coadd
data = np.genfromtxt(pathIn+"f150_power_T_masked.txt")
lCen = data[:,0]
ClMeasured = np.nan_to_num(data[:,1])
sClMeasured = np.nan_to_num(data[:,2])

# Stitch the two spectra at ell=100, to avoid high noise, and interpolate
ClTheory = np.array(map(cmb1_4.flensedTT, lCen))
ClStitched = ClTheory * (lCen<=100.) + ClMeasured * (lCen>100.)
fClStitched = interp1d(lCen, ClStitched, kind='linear', bounds_error=False, fill_value=0.)

# generate power spectrum array for every ell
L = np.arange(0., np.max(lCen))
ClInput = np.array(map(fClStitched, L))

# Map properties needed to generate the GRF
nSide = 4096   # this is what pixell chooses when converting our CAR map to healpix
# read CAR hit map to get the desired pixellation properties
pathHit = pathIn + "f150_daynight_all_div_mono.fits"
hitMap = enmap.read_map(pathHit)
hitShape = hitMap.shape
hitWcs = hitMap.wcs

#########################################################################
# Create the mocks

def genGRF(iMock):
   # set the random seed
   np.random.seed(iMock)
   # Generate GRF healpix map
   hpGrfMap = hp.synfast(ClInput, nSide, lmax=None, mmax=None, alm=False, pol=False, pixwin=False, fwhm=0.0, sigma=None, new=False, verbose=False)
   # save healpix map
   hp.write_map(pathOut+"hp_mock_"+str(iMock)+"_grf_f150_daynight.fits", hpGrfMap, overwrite=True)
   # Convert to pixell CAR map
   grfMap = reproject.enmap_from_healpix(hpGrfMap, hitShape, hitWcs, rot=None)
   # save CAR map
   enmap.write_map(pathOut+"mock_"+str(iMock)+"_grf_f150_daynight.fits", grfMap)


print "Generating mocks"
tStart = time()
with sharedmem.MapReduce(np=nProc) as pool:
   np.array(pool.map(genGRF, range(nMocks)))
tStop = time()
print "Took", (tStop-tStart)/60., "min"


#########################################################################
# Check that the power spectra match the input


def measurePower(iMock):
   # read the map
   hpGrfMap = hp.read_map(pathOut+"hp_mock_"+str(iMock)+"_grf_f150_daynight.fits")
   # measure power spectrum
   Cl = hp.anafast(hpGrfMap)
   Cl = np.nan_to_num(Cl)
   return Cl

print "Checking power spectra"
tStart = time()
with sharedmem.MapReduce(np=nProc) as pool:
   result = np.array(pool.map(measurePower, range(nMocks)))
tStop = time()
print "Took", (tStop-tStart)/60., "min"


# get mean power from mocks
ClMocks = np.mean(result, axis=0)
sClMocks = np.std(result, axis=0) / np.sqrt(nMocks)


fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ell = range(len(ClMocks))
factor = ell*(ell+1.)/(2.*np.pi)
ax.errorbar(ell, factor * ClMocks, yerr=factor * sClMocks, fmt='.', label=r'Mocks')
ax.plot(L, L*(L+1.)/(2.*np.pi) * ClInput, 'k', label=r'Input')
#
for iMock in range(nMocks):
   ax.plot(ell, factor * result[iMock,:], 'g', alpha=0.2)
#
ax.set_ylim((1.e2, 1.e5))
ax.set_xscale('log', nonposx='clip')
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$\ell(\ell+1) C_\ell / (2\pi)$')
#
fig.savefig(pathFig+"power_"+str(nMocks)+"mocks_grf_f150_daynight.pdf", bbox_inches='tight')
fig.clf()
