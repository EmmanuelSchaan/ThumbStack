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



#mask = ts.catalogMask()
#
#filterType= 'diskring'
#
## ts has shape (nObj, nRAp)
## sigmas has shape  nRAp
#sigmas = np.std(ts.filtMap[filterType][mask,:], axis=0)
#
## find the cut for the typical number of objects
#nObj = np.sum(mask)
#f = lambda nSigmas: nObj * special.erfc(nSigmas / np.sqrt(2.)) - special.erfc(5. / np.sqrt(2.))
#nSigmasCut = optimize.brentq(f , 0., 1.e2)
#
## shape is (nObj, nRAp)
#newMask = (np.abs(ts.filtMap[filterType][mask,:]) <= nSigmasCut * sigmas[np.newaxis,:])
## take the intersection of the masks
#newMask = np.prod(newMask, axis=1)

###################################################################################
###################################################################################
# test removing the mean temperature in redshift bins

#ts.measureAllMeanTZBins(plot=True, test=False)



