import numpy as np
from time import time
import matplotlib.pyplot as plt

import universe
reload(universe)
from universe import *

import mass_conversion
reload(mass_conversion)
from mass_conversion import *

import catalog
reload(catalog)
from catalog import *

import flat_map
reload(flat_map)
from flat_map import *

import cmb
reload(cmb)
from cmb import *

#from enlib import enmap, utils, powspec
from pixell import enmap, utils, powspec, enplot, reproject

import healpy as hp

# for cori
plt.switch_backend('Agg')


#########################################################################
#########################################################################
# Mock numbers to generate

# First mock to generate
iMock0 = np.int(sys.argv[1])

# Number of mocks to generate
nMocks = np.int(sys.argv[2])

print("Generating and analyzing  mocks "+str(iMock0)+" to "+str(iMock0+nMocks))

#########################################################################

nProc = 10  # 32 cores on cori, but would run out of memory. 10 works. 15 runs out of memory.

pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/"
pathOut = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/output/cmb_map/mocks_grf_planck_act_coadd_2019_03_11/"
pathFig = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/cmb_map/mocks_grf_planck_act_coadd_2019_03_11/"

if not os.path.exists(pathOut):
   os.makedirs(pathOut)
if not os.path.exists(pathFig):
   os.makedirs(pathFig)

#########################################################################

# theory curve for the CMB power spectrum
cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)

# Load the actual power spectrum of the Planck+ACT+ACTPol coadd
data = np.genfromtxt(pathIn+"f150_power_T_masked.txt")
lCen = data[:,0]
ClMeasured = np.nan_to_num(data[:,1])
sClMeasured = np.nan_to_num(data[:,2])

# Stitch the two spectra at ell=100, to avoid high noise at low ell, and interpolate
ClTheory = np.array(map(cmb1_4.flensedTT, lCen))
ClStitched = ClTheory * (lCen<=100.) + ClMeasured * (lCen>100.)
fClStitched = interp1d(lCen, ClStitched, kind='linear', bounds_error=False, fill_value=0.)

# generate power spectrum array for every ell
L = np.arange(0., np.max(lCen))
Cl = np.array(map(fClStitched, L))

# Map properties needed to generate the GRF
nSide = 4096   # this is what pixell chooses when converting our CAR map to healpix
# read CAR hit map to get the desired pixellation properties
pathHit = pathIn + "act_planck_f150_prelim_div_mono.fits"
hitMap = enmap.read_map(pathHit)
hitShape = hitMap.shape
hitWcs = hitMap.wcs

#########################################################################
# Create the mocks

# Generate healpix map, then convert to CAR
def genGRF(iMock):
   print("- generating mock "+str(iMock))
   # set the random seed
   np.random.seed(iMock)
   # Generate GRF healpix map
   hpGrfMap = hp.synfast(Cl, nSide, lmax=None, mmax=None, alm=False, pol=False, pixwin=False, fwhm=0.0, sigma=None, new=False, verbose=False)
   # save healpix map
   hp.write_map(pathOut+"hp_mock_"+str(iMock)+"_grf_f150_daynight.fits", hpGrfMap, overwrite=True)
   # Convert to pixell CAR map
   grfMap = reproject.enmap_from_healpix(hpGrfMap, hitShape, hitWcs, rot=None)
   # save CAR map
   enmap.write_map(pathOut+"mock_"+str(iMock)+"_grf_f150_daynight.fits", grfMap)


## generate directly a CAR map
#def genGRF(iMock):
   ## set the random seed
   #np.random.seed(iMock)
   ## Set up the cov argument
   ##cov = np.column_stack((L, Cl))
   #cov = Cl.copy()
   ##cov = np.zeros((1,1,len(Cl)))
   ##cov[0,0,:] = Cl.copy()
   ## generate the CAR map
   ##grfMap = enmap.rand_map(hitMap.shape, hitMap.wcs, cov[None,None], scalar=True, pixel_units=False)
   #grfMap = enmap.rand_map(hitMap.shape, hitMap.wcs, cov[None,None])
   ## save CAR map
   #enmap.write_map(pathOut+"mock_"+str(iMock)+"_grf_f150_daynight.fits", grfMap)

'''
print "Generating mocks"
tStart = time()
with sharedmem.MapReduce(np=nProc) as pool:
   np.array(pool.map(genGRF, range(iMock0, iMock0+nMocks)))
#np.array(map(genGRF, range(nMocks)))
tStop = time()
print "Took", (tStop-tStart)/60., "min"
'''

#########################################################################
# Check that the power spectra match the input


def powerSpectrum(hMap, mask=None, theory=[], fsCl=None, nBins=51, lRange=None, plot=False, path="./figures/tests/test_power.pdf", save=False):
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



def measurePower(iMock):
   print "- measuring power for mock", iMock
   # read the map
   hpGrfMap = hp.read_map(pathOut+"hp_mock_"+str(iMock)+"_grf_f150_daynight.fits")
   # measure its power spectrum
   lCen, Cl, sCl = powerSpectrum(hpGrfMap, nBins=101, lRange=None)
   return lCen, Cl, sCl

#def measurePower(iMock):
   #print "- measuring power for mock", iMock
   ## read the map
   #grfMap = enmap.read_map(pathOut+"mock_"+str(iMock)+"_grf_f150_daynight.fits")
   #hpGrfMap = enmap.to_healpix(grfMap)
   ## measure its power spectrum
   #lCen, Cl, sCl = powerSpectrum(hpGrfMap, nBins=101, lRange=None)
   #return lCen, Cl, sCl

'''
print "Checking power spectra"
tStart = time()
with sharedmem.MapReduce(np=nProc) as pool:
   result = np.array(pool.map(measurePower, range(iMock0, iMock0+nMocks)))
#result = np.array(map(measurePower, range(nMocks)))
tStop = time()
print "Took", (tStop-tStart)/60., "min"

# get mean power from mocks
lCen = result[0,0,:]
ClMocks = np.mean(result[:,1,:], axis=0)
sClMocks = np.std(result[:,2,:], axis=0) / np.sqrt(nMocks)

# save to file
data = np.zeros((len(lCen), 3))
data[:,0] = lCen
data[:,1] = ClMocks
data[:,2] = sClMocks
np.savetxt(pathOut + "mean_cl_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".txt", data)
'''

'''
# read mean power spectrum from file
data = np.genfromtxt(pathOut + "mean_cl_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".txt")
lCen = data[:,0]
ClMocks = data[:,1]
sClMocks = data[:,2]

# plot power spectrum
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
factor = lCen*(lCen+1.)/(2.*np.pi)
ax.errorbar(lCen, factor * ClMocks, yerr=factor * sClMocks, fmt='.', label=r'Mocks')
ax.plot(lCen, lCen*(lCen+1.)/(2.*np.pi) * ClStitched, 'k--', label=r'Input')
ax.plot(L, L*(L+1.)/(2.*np.pi) * Cl, 'k', label=r'Input')
#
#for iMock in range(nMocks):
#   ax.plot(lCen, factor * result[iMock,1,:], 'g', alpha=0.2)
#
ax.legend(loc=2)
ax.set_ylim((1.e2, 1.e5))
ax.set_xscale('log', nonposx='clip')
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel(r'$\ell$')
ax.set_ylabel(r'$\ell(\ell+1) C_\ell / (2\pi)$')
#
fig.savefig(pathFig+"power_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+"_grf_f150_daynight.pdf", bbox_inches='tight')
fig.clf()

#plt.show()
'''



##################################################################################
##################################################################################
##################################################################################
# do the stacking on the mock maps

nProc = 32 # 1 haswell node on cori

# cosmological parameters
u = UnivMariana()

# M*-Mh relation
massConversion = MassConversionKravtsov14()
#massConversion.plot()

###################################################################################
# Galaxy catalogs

###################################################################################
# Mariana

# CMASS
#cmassSMariana = Catalog(u, massConversion, name="cmass_s_mariana", nameLong="CMASS S M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt", save=False)
#cmassNMariana = Catalog(u, massConversion, name="cmass_n_mariana", nameLong="CMASS N M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt", save=False)
# combined catalog
#cmassMariana = cmassSMariana.copy(name="cmass_mariana", nameLong="CMASS M")
#cmassMariana.addCatalog(cmassNMariana, save=True)
cmassMariana = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False)


###################################################################################
###################################################################################
# Read CMB mask and hit count

# path to true hit count map and mask: Planck + ACT
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/"
pathHit = pathIn + "act_planck_f150_prelim_div_mono.fits"
pathMask = pathIn + "f150_mask_foot_planck_ps_car.fits"
pathPower = pathIn + "f150_power_T_masked.txt"

# read maps in common for all mocks
pactMask = enmap.read_map(pathMask)
pactHit = enmap.read_map(pathHit)

# theory power spectrum
cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)


###################################################################################
###################################################################################
# Do the stacking on all the mock GRF CMB maps

import thumbstack
reload(thumbstack)
from thumbstack import *

# recompute or not the stacked profiles
save = True


def doStacking(iMock):
   print("Stacking on mock "+str(iMock))
   pathMap = pathOut+"mock_"+str(iMock)+"_grf_f150_daynight.fits"
   pactMap = enmap.read_map(pathMap)
   name = cmassMariana.name + "_mockgrf"+str(iMock)+"_pactf150night20190311"
   ts = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=save, nProc=nProc)
   # output the various stacked profiles (tSZ, kSZ, etc.) from this mock GRF CMB map
   return ts.stackedProfile


print "Stacking on each mock map"
tStart = time()
result = np.array(map(doStacking, range(iMock0, iMock0+nMocks)))
#result = np.array(map(doStacking, [199, 298, 299, 398, 399]))
tStop = time()
print "Finished all stacking: took", (tStop-tStart)/60., "min"



###################################################################################
# Estimate the covariances from the mocks

# for debugging only:
#nMocks = 10

# Load one thumbstack object to access the plotting routines
pathMap = pathOut+"mock_"+str(iMock0)+"_grf_f150_daynight.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_mockgrf"+str(iMock0)+"_pactf150night20190311"
ts = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=False, nProc=nProc)


Est = ['tsz_uniformweight', 'tsz_hitweight', 'tsz_varweight', 'ksz_uniformweight', 'ksz_hitweight', 'ksz_varweight', 'ksz_massvarweight']
covStackedProfile = {}
meanStackedProfile = {}
for iEst in range(len(Est)):
#for iEst in range(1):
   est = Est[iEst]
   nRAp = len(result[0][est])
   # shape (nRAp, nMocks)
   profiles = np.array([[result[j][est][i] for j in range(nMocks)] for i in range(nRAp)])
   
   # estimate and save the mean
   meanStackedProfile[est] = np.mean(profiles, axis=-1)
   np.savetxt(pathOut+'mean_'+est+'_mocks'+str(iMock0)+"-"+str(iMock0+nMocks)+'.txt', meanStackedProfile[est])

   # estimate and save the cov
   covStackedProfile[est] = np.cov(profiles, rowvar=True)
   np.savetxt(pathOut+"cov_"+est+"_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".txt", covStackedProfile[est])
   # plot it
   ts.plotCov(covStackedProfile[est], name=est+"_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".pdf")


###################################################################################
# Load and plot the mean and covariances from the mocks


Est = ['tsz_uniformweight', 'tsz_hitweight', 'tsz_varweight', 'ksz_uniformweight', 'ksz_hitweight', 'ksz_varweight', 'ksz_massvarweight']
covStackedProfile = {}
for iEst in range(len(Est)):
#for iEst in [2]:
   est = Est[iEst]
   meanStackedProfile[est] = np.genfromtxt(pathOut+"mean_"+est+"_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".txt")
   covStackedProfile[est] = np.genfromtxt(pathOut+"cov_"+est+"_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".txt")


   # Compare diagonal covariance elements
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   # convert from sr to arcmin^2
   factor = (180.*60./np.pi)**2
   #
   sMocks = np.sqrt(np.diag(covStackedProfile[est]))
   #ax.errorbar(ts.RApArcmin, factor * sMocks, yerr= factor * sMocks / np.sqrt(nMocks),  label=r'mocks')
   ax.fill_between(ts.RApArcmin, factor * sMocks * (1.-1./np.sqrt(nMocks)), factor * sMocks * (1.+1./np.sqrt(nMocks)), facecolor='m', edgecolor='', alpha=0.5, label=r'mocks')
   ax.plot(ts.RApArcmin, factor * ts.sStackedProfile[est], label=r'semi-analytical')
   if est in ts.covBootstrap.keys():
      sBootstrap = np.sqrt(np.diag(ts.covBootstrap[est]))
      ax.plot(ts.RApArcmin, factor * sBootstrap, label=r'bootstrap')
   #
   ax.legend(loc=2)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(est.replace('_', ' '))
   #
   fig.savefig(pathFig+'/compare_std_'+est+'_mocks'+str(iMock0)+"-"+str(iMock0+nMocks)+'.pdf', bbox_inches='tight')
   #plt.show()
   fig.clf()



# Null test: plot the mean
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# convert from sr to arcmin^2
factor = (180.*60./np.pi)**2
#
for iEst in range(len(Est)):
   est = Est[iEst]
   sMocks = np.sqrt(np.diag(covStackedProfile[est]))
   ax.errorbar(ts.RApArcmin*(1.+0.02*iEst/len(Est)), factor * meanStackedProfile[est], yerr= factor * sMocks / np.sqrt(nMocks),  label=est.replace('_', ' '))
#
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T$ [$\mu K\cdot\text{arcmin}^2$]')
#
fig.savefig(pathFig+'/nulltests_mocks'+str(iMock0)+"-"+str(iMock0+nMocks)+'.pdf', bbox_inches='tight')
#plt.show()
fig.clf()

###################################################################################
