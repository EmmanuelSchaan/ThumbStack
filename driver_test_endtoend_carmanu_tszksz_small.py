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


# Make mock map using the full pactMap geometry
# path to true hit count map and mask: Planck + ACT
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/"
pathMask = pathIn + "f150_mask_foot_planck_ps_car.fits"
# read maps in common for all mocks
boxMask = enmap.read_map(pathMask)[0]




## Cannot upgrade map: one cori node runs out of memory when generating the catalog
## Increase the map resolution for testing purposes
#print("Upgrade map resolution")
#print("Initial dimensions "+str(boxMask.shape))
#boxMask = boxMask.upgrade(2)
#print("Final dimensions "+str(boxMask.shape))

# Downgrade the map resolution for testing purposes
print("Downgrade map resolution")
print("Initial dimensions "+str(boxMask.shape))
boxMask = boxMask.downgrade(2)
print("Final dimensions "+str(boxMask.shape))


## Generate empty square map, then make a mock map
## Generate an empty square map with RA in [0., 5.] and DEC in [-5., 0.] 
#box = np.array([[-5., 5.], [0., 0.]]) * utils.degree
##box = np.array([[-5., 0.], [0., 5.]]) * utils.degree
#resArcmin = 0.25 #1. #0.5  # 0.1   # map pixel size [arcmin]
#shape,wcs = enmap.geometry(pos=box, res=resArcmin * utils.arcmin, proj='car')
## create a mask that keeps the whole area
##boxMask = enmap.zeros((1,) + shape, wcs=wcs)
#boxMask = enmap.ones(shape, wcs=wcs)


# create mock map with point sources and Gaussian profiles with sigma=1.5'
cmassMarianaVShuffleSmall.generateMockMaps(boxMask, sigma=1.5, test=False)

# check that the mock map has non-zero pixels
#pathMap = cmassMarianaVShuffleSmall.pathOut + "mock_count_dirac_car.fits"
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



# Do it on mocks with shuffled velocities
pathMap = cmassMarianaVShuffleSmall.pathOut + "mock_vel_dirac_car.fits"
boxMap = enmap.read_map(pathMap)
name = cmassMarianaVShuffleSmall.name + "_test_endtoend_vel_dirac_carmanu"
tsVelDiracVShuffleSmall = ThumbStack(u, cmassMarianaVShuffleSmall, boxMap, boxMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)

pathMap = cmassMarianaVShuffleSmall.pathOut + "mock_vel_gauss_car.fits"
boxMap = enmap.read_map(pathMap)
name = cmassMarianaVShuffleSmall.name + "_test_endtoend_vel_gauss_carmanu"
tsVelGaussVShuffleSmall = ThumbStack(u, cmassMarianaVShuffleSmall, boxMap, boxMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)




###################################################################################



# Gaussian with sigma = 1.5'
s1 = 0.76   # 0.6
s2 = 1.5 # 1.5
s3 = 1.68   # 1.61
profile1 = tsVelDiracVShuffleSmall.ftheoryGaussianProfile(sigma_cluster=s1) # 0.6
profile2 = tsVelDiracVShuffleSmall.ftheoryGaussianProfile(sigma_cluster=s2) # 1.5
profile3 = tsVelDiracVShuffleSmall.ftheoryGaussianProfile(sigma_cluster=s3) # 1.61
#profilePix = tsVelDiracVShuffleSmall.ftheoryGaussianProfilePixelated(sigma_cluster=1.5, resArcmin=resArcmin) # 1.61
#profilePixPixwin = tsVelDiracVShuffleSmall.ftheoryGaussianProfilePixelated(sigma_cluster=1.5, pixwin=1, resArcmin=resArcmin) # 1.5

fig=plt.figure(0)
ax=fig.add_subplot(111)
#
factor =  (180.*60./np.pi)**2
ax.errorbar(tsVelDiracVShuffleSmall.RApArcmin, factor*tsVelDiracVShuffleSmall.stackedProfile['ksz_uniformweight'], factor*tsVelDiracVShuffleSmall.sStackedProfile['ksz_uniformweight'], fmt='--', c='g', label=r'vel Dirac v-shuffle')
ax.errorbar(tsVelGaussVShuffleSmall.RApArcmin, factor*tsVelGaussVShuffleSmall.stackedProfile['ksz_uniformweight'], factor*tsVelGaussVShuffleSmall.sStackedProfile['ksz_uniformweight'], fmt='-', c='g', label=r'vel Gauss v-shuffle')
#
ax.plot(tsVelDiracVShuffleSmall.RApArcmin, profile1, 'k-', label=r'expected '+str(s1))
ax.plot(tsVelDiracVShuffleSmall.RApArcmin, profile2, 'k-', label=r'expected '+str(s2))
ax.plot(tsVelDiracVShuffleSmall.RApArcmin, profile3, 'k-', label=r'expected '+str(s3))
#ax.plot(tsVelDiracVShuffleSmall.RApArcmin, profilePix, 'k--', label=r'expected, pixelated')
#ax.plot(tsVelDiracVShuffleSmall.RApArcmin, profilePixPixwin, 'c-', label=r'expected, pixelated, pixwin 1')
#
ax.legend(loc=4, fontsize='x-small', labelspacing=0.1)
#
fig.savefig(tsVelDiracVShuffleSmall.pathFig+"/test_mean_stacked_temperature.pdf")
#fig.clf()

plt.show()

# asymptote is 1.92858407e-08


