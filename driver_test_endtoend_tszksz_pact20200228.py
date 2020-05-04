# Last full run for VelDiraacVShuffle and VelGaussVShuffle: 2020/01/01

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
cmassMariana = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False)


# Shuffle velocities to kill the 2-halo term
#cmassMarianaVShuffle = cmassMariana.copy(name="cmass_mariana_vshuffle", nameLong="CMASS M Vshuffle")
#np.random.shuffle(cmassMarianaVShuffle.vR)
#cmassMarianaVShuffle.writeCatalog()
cmassMarianaVShuffle = Catalog(u, massConversion, name="cmass_mariana_vshuffle", nameLong="CMASS M Vshuffle", save=False)



###################################################################################

# keep only the first few objects, to speed things up
#I = range(10000)
#cmassMariana = cmassMariana.extractCatalog(I, name="mini_cmass_mariana", nameLong="mini CMASS M")
#cmassMariana = Catalog(u, massConversion, name="mini_cmass_mariana", nameLong="mini CMASS M", save=False) 

# Shuffle velocities to kill the 2-halo term
#cmassMarianaVShuffle = cmassMariana.copy(name="mini_cmass_mariana_vshuffle", nameLong="mini CMASS M Vshuffle")
#np.random.shuffle(cmassMarianaVShuffle.vR)
#cmassMarianaVShuffle.writeCatalog()
#cmassMarianaVShuffle = Catalog(u, massConversion, name="mini_cmass_mariana_vshuffle", nameLong="mini CMASS M Vshuffle", save=False)


###################################################################################
###################################################################################
# Read CMB mask and hit count

# Hit count
pathHit = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits"
pactHit = enmap.read_map(pathHit)[0]


###################################################################################
###################################################################################
# Generate mock maps

# Point sources and Gaussian profiles with sigma=1.5'
#cmassMariana.generateMockMaps(pactHit, sigma=1.5)

# Same for the catalog with shuffled velocities
#cmassMarianaVShuffle.generateMockMaps(pactHit, sigma=1.5)


###################################################################################
###################################################################################

# CMB Mask
pathMask = "./output/cmb_map/pact20200228/" + "mask_full_foot_gal60_ps.fits"
pactMask = enmap.read_map(pathMask)


###################################################################################
###################################################################################

import thumbstack
reload(thumbstack)
from thumbstack import *


###################################################################################

save = False

pathMap = cmassMariana.pathOut + "mock_count_dirac_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_count_dirac_carmanu"
tsCountDirac = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=save, filterTypes='diskring', nProc=nProc)

pathMap = cmassMariana.pathOut + "mock_count_gauss_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_count_gauss_carmanu"
tsCountGauss = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=save, filterTypes='diskring', nProc=nProc)

pathMap = cmassMariana.pathOut + "mock_vel_dirac_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_vel_dirac_carmanu"
tsVelDirac = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=save, filterTypes='diskring', nProc=nProc)

pathMap = cmassMariana.pathOut + "mock_vel_gauss_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_vel_gauss_carmanu"
tsVelGauss = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=save, filterTypes='diskring', nProc=nProc)


###################################################################################
# Same on mocks with shuffled velocities

pathMap = cmassMarianaVShuffle.pathOut + "mock_vel_dirac_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMarianaVShuffle.name + "_pactf150night20190311_test_endtoend_vel_dirac_carmanu"
tsVelDiracVShuffle = ThumbStack(u, cmassMarianaVShuffle, pactMap, pactMask, pactHit, name=name, nameLong=None, save=save, filterTypes='diskring', nProc=nProc)

pathMap = cmassMarianaVShuffle.pathOut + "mock_vel_gauss_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMarianaVShuffle.name + "_pactf150night20190311_test_endtoend_vel_gauss_carmanu"
tsVelGaussVShuffle = ThumbStack(u, cmassMarianaVShuffle, pactMap, pactMask, pactHit, name=name, nameLong=None, save=save, filterTypes='diskring', nProc=nProc)



###################################################################################
###################################################################################

filterType = 'diskring'

# Theory Gaussian with sigma = 1.5'
s1 = 0.25   # 0.6
s2 = 1.5 # 1.5
s3 = 1.61   # 1.61
#profile1 = tsVelDiracVShuffle.ftheoryGaussianProfile(sigma_cluster=s1, filterType=filterType) # 0.6
profile2 = tsVelDiracVShuffle.ftheoryGaussianProfile(sigma_cluster=s2, filterType=filterType) # 1.5
#profile3 = tsVelDiracVShuffle.ftheoryGaussianProfile(sigma_cluster=s3, filterType=filterType) # 1.61
# Gaussian with sigma = 1.5'
#profile = tsVelDiracVShuffle.ftheoryGaussianProfile(sigma_cluster=1.5) # 1.61
#profilePix = tsVelDiracVShuffle.ftheoryGaussianProfilePixelated(sigma_cluster=1.5, resArcmin=0.5) # 1.61
#profilePixPixwin1 = tsVelDiracVShuffle.ftheoryGaussianProfilePixelated(sigma_cluster=1.5, resArcmin=0.5, pixwin=1) # 1.5
#profilePixPixwin3 = tsVelDiracVShuffle.ftheoryGaussianProfilePixelated(sigma_cluster=1.5, resArcmin=0.5, pixwin=3) # 1.5



###################################################################################
# Get statistical uncertainty on kSZ and tSZ,
# from the main analysis

path = "./output/thumbstack/cmass_mariana_pactf150daynight20200228maskgal60"+"/cov_diskring_ksz_varweight_bootstrap.txt"
sKsz = np.sqrt(np.diag(np.genfromtxt(path)))
path = "./output/thumbstack/cmass_mariana_pactf150daynight20200228maskgal60"+"/cov_diskring_tsz_varweight_bootstrap.txt"
sTsz = np.sqrt(np.diag(np.genfromtxt(path)))



# get the expected integrated tSZ and kSZ in [muK*arcmin^2]

# get the expected integrated y in [muK*arcmin^2]
# expected integrated y in sr
tSZ = np.mean(cmassMariana.integratedY[np.where(cmassMariana.hasM)[0]])
# convert from y profile to dT profile
Tcmb = 2.726   # K
h = 6.63e-34   # SI
kB = 1.38e-23  # SI
def f(nu):
   """frequency dependence for tSZ temperature
   """
   x = h*nu/(kB*Tcmb)
   return x*(np.exp(x)+1.)/(np.exp(x)-1.) -4.
tSZ *= f(150.e9) * Tcmb * 1.e6  # [muK * sr]
# convert to [muK*arcmin^2]
tSZ *= (180.*60./np.pi)**2

# get the expected integrated kSZ in [muK*arcmin^2]
# expected integrated kSZ in muK*sr
kSZ = np.std(cmassMariana.integratedKSZ[np.where(cmassMariana.hasM)[0]])
# convert to [muK*arcmin^2]
kSZ *= (180.*60./np.pi)**2

print "EXpected kSZ and tSZ in [muK*arcmin^2] at 150GHz:"
print kSZ, tSZ


# Rescale the error bars to be relevant for the mocks
# where the kSZ and tSZ integrate to 1 muK*arcmin^2
sKsz /= kSZ
sTsz /= tSZ

###################################################################################



fig=plt.figure(0)
ax=fig.add_subplot(111)
#
factor =  (180.*60./np.pi)**2
#
# Expected signal for an isolated Gaussian
#ax.plot(tsVelDiracVShuffle.RApArcmin, profile1, 'k-', label=r'theory '+str(s1))
ax.plot(tsVelDiracVShuffle.RApArcmin, profile2, 'k-', lw=3, label=r'Theory 1h, Gaussian')# '+str(s2))
#ax.plot(tsVelDiracVShuffle.RApArcmin, profile3, 'k-', label=r'theory '+str(s3))
#ax.plot(tsVelDiracVShuffle.RApArcmin, profile, 'k-', label=r'expected')
#ax.plot(tsVelDiracVShuffle.RApArcmin, profilePix, 'k--', label=r'expected, pixelated')
#ax.plot(tsVelDiracVShuffle.RApArcmin, profilePixPixwin1, 'c-', label=r'expected, pixelated, pixwin 1')
#ax.plot(tsVelDiracVShuffle.RApArcmin, profilePixPixwin3, 'm-', label=r'expected, pixelated, pixwin 3')
#
# Statistical uncertainty on tSZ and kSZ, from main analysis
ax.fill_between(tsVelDirac.RApArcmin, factor * (tsVelGauss.stackedProfile['diskring_ksz_uniformweight'] - sKsz), factor * (tsVelGauss.stackedProfile['diskring_ksz_uniformweight'] + sKsz), facecolor='b', edgecolor='', alpha=0.05)
ax.fill_between(tsVelDirac.RApArcmin, factor * (tsCountGauss.stackedProfile['diskring_tsz_uniformweight'] - sTsz), factor * (tsCountGauss.stackedProfile['diskring_tsz_uniformweight'] + sTsz), facecolor='r', edgecolor='', alpha=0.05)
#
# Measured 1h and 2h
ax.plot(tsVelDiracVShuffle.RApArcmin, factor*tsVelDiracVShuffle.stackedProfile['diskring_ksz_uniformweight'], 'g--', label=r'1h only, Pointlike')
ax.plot(tsVelGaussVShuffle.RApArcmin, factor*tsVelGaussVShuffle.stackedProfile['diskring_ksz_uniformweight'], 'g', label=r'1h only Gaussian')
ax.plot(tsCountDirac.RApArcmin, factor*tsCountDirac.stackedProfile['diskring_tsz_uniformweight'], 'r--', label=r'tSZ Pointlike')
ax.plot(tsCountGauss.RApArcmin, factor*tsCountGauss.stackedProfile['diskring_tsz_uniformweight'], 'r', label=r'tSZ Gaussian')
ax.plot(tsVelDirac.RApArcmin, factor*tsVelDirac.stackedProfile['diskring_ksz_uniformweight'], 'b--', label=r'kSZ Pointlike')
ax.plot(tsVelGauss.RApArcmin, factor*tsVelGauss.stackedProfile['diskring_ksz_uniformweight'], 'b', label=r'kSZ Gaussian')
#
#
ax.legend(loc=4, fontsize='x-small', labelspacing=0.1, handlelength=1.)
ax.set_ylim((0., 1.3))
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T$ [normalized to unity]')
#
print "saving figure to:"
print tsCountDirac.pathFig+"/test_mean_stacked_temperature_diskring_full.pdf"
fig.savefig(tsCountDirac.pathFig+"/test_mean_stacked_temperature_diskring_full.pdf", bbox_inches='tight')
fig.clf()

#plt.show()



