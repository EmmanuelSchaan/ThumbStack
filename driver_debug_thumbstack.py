import universe
reload(universe)
from universe import *

import mass_conversion
reload(mass_conversion)
from mass_conversion import *

import catalog
reload(catalog)
from catalog import *

import thumbstack
reload(thumbstack)
from thumbstack import *

import cmb
reload(cmb)
from cmb import *

# running on cori
# 68 cores per knl node, 32 cores per haswell node
#salloc -N 1 --qos=interactive -C haswell -t 04:00:00 -L SCRATCH
#plt.switch_backend('Agg')

##################################################################################

nProc = 32 # 1 haswell node on cori

##################################################################################
##################################################################################
# cosmological parameters

u = UnivMariana()

###################################################################################
# M*-Mh relation

massConversion = MassConversionKravtsov14()
#massConversion.plot()

###################################################################################
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
#cmassMariana = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False)


# Shuffle velocities to kill the 2-halo term
#cmassMarianaVShuffle = cmassMariana.copy(name="cmass_mariana_vshuffle", nameLong="CMASS M Vshuffle")
#np.random.shuffle(cmassMarianaVShuffle.vR)
#cmassMarianaVShuffle.writeCatalog()
#cmassMarianaVShuffle = Catalog(u, massConversion, name="cmass_mariana_vshuffle", nameLong="CMASS M Vshuffle", save=False)

# extract a small square patch with RA in [0., 5.] and DEC in [-5., 0.],
# to speed up the tests: reduces nObj from 777202 to 2172.
#I = np.where((0.<=cmassMarianaVShuffle.RA) * (cmassMarianaVShuffle.RA<=5.) * (-5.<=cmassMarianaVShuffle.DEC) * (cmassMarianaVShuffle.DEC<=0.))
#cmassMarianaVShuffle.extractCatalog(I, name="cmass_mariana_vshuffle_small", nameLong="CMASS M Vshuffle small")
cmassMarianaVShuffleSmall = Catalog(u, massConversion, name="cmass_mariana_vshuffle_small", nameLong="CMASS M Vshuffle small")


###################################################################################
# Generate mock, starting either from a square box map
# or from the full ACT map geometry


## Make mock map using the full pactMap geometry
## path to true hit count map and mask: Planck + ACT
#pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/"
#pathMask = pathIn + "f150_mask_foot_planck_ps_car.fits"
## read maps in common for all mocks
#boxMask = enmap.read_map(pathMask)[0]




## Cannot upgrade map: one cori node runs out of memory when generating the catalog
## Increase the map resolution for testing purposes
#print("Upgrade map resolution")
#print("Initial dimensions "+str(boxMask.shape))
#boxMask = boxMask.upgrade(2)
#print("Final dimensions "+str(boxMask.shape))

# Downgrade the map resolution for testing purposes
#print("Downgrade map resolution")
#print("Initial dimensions "+str(boxMask.shape))
#boxMask = boxMask.downgrade(2)
#print("Final dimensions "+str(boxMask.shape))


# Generate empty square map, then make a mock map
# Generate an empty square map with RA in [0., 5.] and DEC in [-5., 0.] 
box = np.array([[-2., 2.], [0., 0.]]) * utils.degree
#box = np.array([[-5., 0.], [0., 5.]]) * utils.degree
resArcmin = 0.5 #1. #0.5  # 0.1   # map pixel size [arcmin]
shape,wcs = enmap.geometry(pos=box, res=resArcmin * utils.arcmin, proj='car')
# create a mask that keeps the whole area
#boxMask = enmap.zeros((1,) + shape, wcs=wcs)
boxMask = enmap.ones(shape, wcs=wcs)


# create mock map with point sources and Gaussian profiles with sigma=1.5'
#cmassMarianaVShuffleSmall.generateMockMaps(boxMask, sigma=1.5, test=False)

# check that the mock map has non-zero pixels
#pathMap = cmassMarianaVShuffleSmall.pathOut + "mock_count_gauss_car.fits"
#boxMap = enmap.read_map(pathMap)
#print np.sum(np.abs(boxMap))


###################################################################################
###################################################################################
# Read CMB map and mask

# theory power spectrum
cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)


###################################################################################
###################################################################################

import thumbstack
reload(thumbstack)
from thumbstack import *

save = False

# Do it on mocks with shuffled velocities
pathMap = cmassMarianaVShuffleSmall.pathOut + "mock_vel_dirac_car.fits"
boxMap = enmap.read_map(pathMap)
name = cmassMarianaVShuffleSmall.name + "_test_endtoend_vel_dirac_carmanu"
#tsVelDiracVShuffleSmall = ThumbStack(u, cmassMarianaVShuffleSmall, boxMap, boxMask, cmbHit=None, name=name, nameLong=None, save=save, nProc=nProc, filterTypes='all')
tsVelDiracVShuffleSmall = ThumbStack(u, cmassMarianaVShuffleSmall, boxMap, boxMask, cmbHit=None, name=name, nameLong=None, save=save, nProc=nProc, filterTypes='diskring', doMBins=True)

ts = tsVelDiracVShuffleSmall

#ts.plotTszKszContaminationMMax()
#mMax = 1.e15
#prof, sprof  = ts.computeStackedProfile('diskring', 'ksz_varweight', mVir=[1.e6, mMax]) # [map unit * sr]




'''
def checkFilterHistograms(filterType, est, mVir=None, z=[0., 100.], ts=None):
   if ts is None:
      ts = ts
   if mVir is None:
      mVir = [ts.mMin, ts.mMax]

   # select objects that overlap, and reject point sources
   mask = ts.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=mVir, z=z)

   # aperture photometry filters
   t = ts.filtMap[filterType].copy() # [muK * sr]
   t = t[mask, :]
   # convert from sr to arcmin^2
   factor = (180.*60./np.pi)**2
   t *= factor

   #true filter variance for each object and aperture,
   # valid whether or not a hit count map is available
   s2Full = ts.filtVarTrue[filterType][mask, :]
   # Variance from hit count (if available)
   s2Hit = ts.filtHitNoiseStdDev[filterType][mask, :]**2


   for iRAp in range(ts.nRAp):

      path = ts.pathFig + "/histogram_t_"+filterType+"_uniformweight_"+str(iRAp)+".pdf"
      myHistogram(t[:,iRAp], nBins=101, lim=None, S2Theory=[], path=path, plot=False, nameLatex=r'$T_i$ [$\mu$K$\cdot$arcmin$^2$]', semilogx=False, semilogy=True, doGauss=True)

      path = ts.pathFig + "/histogram_t_"+filterType+"_varweight_"+str(iRAp)+".pdf"
      myHistogram(t[:,iRAp] / s2Full[:,iRAp], nBins=101, lim=None, S2Theory=[], path=path, plot=False, nameLatex=r'$T_i/\sigma_T^2_i$', semilogx=False, semilogy=True, doGauss=True)


checkFilterHistograms('diskring', 'ksz_varweight', mVir=None, z=[0., 100.], ts=ts)
'''
