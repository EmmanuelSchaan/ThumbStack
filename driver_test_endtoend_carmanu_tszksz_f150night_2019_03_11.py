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
cmassSMariana = Catalog(u, massConversion, name="cmass_s_mariana", nameLong="CMASS S M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt", save=False)
#cmassSMariana.plotHistograms()
#cmassSMariana.plotFootprint()
#cmassMariana.printProperties()
#
cmassNMariana = Catalog(u, massConversion, name="cmass_n_mariana", nameLong="CMASS N M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt", save=False)
#cmassNMariana.plotHistograms()
#cmassNMariana.plotFootprint()
#cmassMariana.printProperties()
#
# combined catalog
cmassMariana = cmassSMariana.copy(name="cmass_mariana", nameLong="CMASS M")
cmassMariana.addCatalog(cmassNMariana, save=False)
#cmassMariana.plotHistograms()
#cmassMariana.plotFootprint()
#cmassMariana.printProperties()


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
# Generate mock maps for CMASS Mariana

# Gaussian profiles with sigma=1.5'
cmassMariana.generateMockMaps(pactHit, sigma=1.5)


###################################################################################
###################################################################################

import thumbstack
reload(thumbstack)
from thumbstack import *



pathMap = cmassMariana.pathOut + "mock_count_dirac_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_count_dirac_carmanu"
tsCountDirac = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=False, nProc=nProc)

pathMap = cmassMariana.pathOut + "mock_count_gauss_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_count_gauss_carmanu"
tsCountGauss = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=False, nProc=nProc)

pathMap = cmassMariana.pathOut + "mock_vel_dirac_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_vel_dirac_carmanu"
tsVelDirac = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=True, nProc=nProc)

pathMap = cmassMariana.pathOut + "mock_vel_gauss_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_vel_gauss_carmanu"
tsVelGauss = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=True, nProc=nProc)




###################################################################################


'''
mask = tsCmassM.catalogMask(overlap=True, psMask=True, mVir=[1.e6, 1.e17], extraSelection=1.)
print np.where(mask==True)

iObj = 0
filtMap, filtMask, filtNoiseStdDev, diskArea = tsCmassM.analyzeObject(iObj, test=True)
print filtMap / hp.pixelfunc.nside2pixarea(8192)
'''










#tszMeas = tsCmassM.stackedProfile['tsz_uniformweight']
#tszTh = tsCmassM.stackedProfile['tsz_uniformweight_theory_tsz']


#print  tszTh / tszMeas

# Gaussian with sigma = 1.5'
profile = tsCountDirac.ftheoryGaussianProfile(1.5)

fig=plt.figure(0)
ax=fig.add_subplot(111)
#
factor =  (180.*60./np.pi)**2
ax.errorbar(tsCountDirac.RApArcmin, factor*tsCountDirac.stackedProfile['tsz_uniformweight'], factor*tsCountDirac.sStackedProfile['tsz_uniformweight'], fmt='--', c='r', label=r'count Dirac')
ax.errorbar(tsCountGauss.RApArcmin, factor*tsCountGauss.stackedProfile['tsz_uniformweight'], factor*tsCountGauss.sStackedProfile['tsz_uniformweight'], fmt='-', c='r', label=r'count Gauss')
ax.errorbar(tsVelDirac.RApArcmin, factor*tsVelDirac.stackedProfile['ksz_uniformweight'], factor*tsVelDirac.sStackedProfile['ksz_uniformweight'], fmt='--', c='b', label=r'vel Dirac')
ax.errorbar(tsVelGauss.RApArcmin, factor*tsVelGauss.stackedProfile['ksz_uniformweight'], factor*tsVelGauss.sStackedProfile['ksz_uniformweight'], fmt='-', c='b', label=r'vel Gauss')
#
ax.plot(tsCountDirac.RApArcmin, profile, 'k-', label=r'expected')
#
ax.legend(loc=2)
#
fig.savefig(tsCountDirac.pathFig+"/test_mean_stacked_temperature.pdf")
#fig.clf()

plt.show()

# asymptote is 1.92858407e-08
