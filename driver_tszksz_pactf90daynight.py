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
# Kendrick

# CMASS
cmassSKendrick = Catalog(u, massConversion, name="cmass_s_kendrick", nameLong="CMASS S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_S_kendrick.txt", save=False)
#cmassSKendrick.plotHistograms()
#cmassSKendrick.plotFootprint()
#
cmassNKendrick = Catalog(u, massConversion, name="cmass_n_kendrick", nameLong="CMASS N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_N_kendrick.txt", save=False)
#cmassNKendrick.plotHistograms()
#cmassNKendrick.plotFootprint()
#
# combined catalog
cmassKendrick = cmassSKendrick.copy(name="cmass_kendrick", nameLong="CMASS K")
cmassKendrick.addCatalog(cmassNKendrick, save=False)
#cmassKendrick.plotHistograms()
#cmassKendrick.plotFootprint()

# LOWZ
lowzSKendrick = Catalog(u, massConversion, name="lowz_s_kendrick", nameLong="LOWZ S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_S_kendrick.txt", save=False)
#lowzSKendrick.plotHistograms()
#lowzSKendrick.plotFootprint()
#
lowzNKendrick = Catalog(u, massConversion, name="lowz_n_kendrick", nameLong="LOWZ N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_N_kendrick.txt", save=False)
#lowzNKendrick.plotHistograms()
#lowzNKendrick.plotFootprint()
#
# combined catalog
lowzKendrick = lowzSKendrick.copy(name="lowz_kendrick", nameLong="LOWZ K")
lowzKendrick.addCatalog(lowzNKendrick, save=False)
#lowzKendrick.plotHistograms()
#lowzKendrick.plotFootprint()

# BOSS = CMASS + LOWZ
bossKendrick = cmassSKendrick.copy(name="boss_kendrick", nameLong="BOSS K")
bossKendrick.addCatalog(cmassNKendrick, save=False)
bossKendrick.addCatalog(lowzSKendrick, save=False)
bossKendrick.addCatalog(lowzNKendrick, save=False)
#bossKendrick.plotHistograms()
#bossKendrick.plotFootprint()



###################################################################################
###################################################################################
# Read CMB maps

# Planck + ACT 90GHz day and night
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/"
pathMap = pathIn + "f90_daynight_all_map_mono.fits"
pathHit = pathIn + "f90_daynight_all_div_mono.fits"
pathMask = pathIn + "f90_mask_foot_planck_ps_car.fits"
pathPower = pathIn + "f90_power_T_masked.txt"

tStart = time()
print "- Read CMB map, mask and hit count"
pact90Map = enmap.read_map(pathMap)[0]   # keep only temperature
pact90Mask = enmap.read_map(pathMask)
pact90Hit = enmap.read_map(pathHit)

print "Set up interpolations"
#!!!! NO! This prefiltering seems to mess up the map values, in a complicated way...
# This pre-filtering step introduces some ringing in the maps

#utils.interpol_prefilter(pact90Map, inplace=True)
#utils.interpol_prefilter(pact90Mask, inplace=True)   # order=2 seems to reduce ringing at sharp transitions
#utils.interpol_prefilter(pact90Hit, inplace=True)

tStop = time()
print "took", tStop-tStart, "sec"

# measured power spectrum
data = np.genfromtxt(pathPower)  # l, Cl, sCl
data = np.nan_to_num(data)
fCl = interp1d(data[:,0], data[:,1], kind='linear', bounds_error=False, fill_value=0.)

# theory power spectrum
cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)


###################################################################################
###################################################################################
# Stacking


import thumbstack
reload(thumbstack)
from thumbstack import *


name = cmassMariana.name + "_pactf90daynight"
tsCmassM = ThumbStack(u, cmassMariana, pact90Map, pact90Mask, pact90Hit, name=name, nameLong=None, save=True, nProc=nProc)


#mask = tsCmassM.catalogMask(overlap=True, psMask=True)
#tsCmassM.analyzeObject(0, test=True)
#tsCmassM.examineCmbMaps()

# Expected std dev of AP filter, function of disk radius in rad
#fsApActual = lambda r0: cmb1_4.fsigmaDiskRing(r0, thetaIn=None, thetaOut=None, fCl=fCl, lMin=1., lMax=1.e5)
#fsApNoiseless = lambda r0: cmb1_4.fsigmaDiskRing(r0, thetaIn=None, thetaOut=None, fCl=cmb1_4.flensedTT, lMin=1., lMax=1.e5)
#tsCmassM.examineHistograms(fsAp=[fsApActual, fsApNoiseless])

#tsCmassM.measureVarFromHitCount(plot=True)

#cov = tsCmassM.kszCovBootstrap(nSamples=1000, nProc=nProc)
#print np.sqrt(np.diag(cov))


tsCmassM.compareKszEstimators()
tsCmassM.computeSnrTsz()
tsCmassM.computeSnrKsz()
tsCmassM.kszNullTests()
tsCmassM.plotCovTszKsz()



###################################################################################

name = cmassSMariana.name + "_pactf90daynight"
tsCmassSM = ThumbStack(u, cmassSMariana, pact90Map, pact90Mask, pact90Hit, name=name, nameLong=None, save=True, nProc=nProc)
#tsCmassSM.compareKszEstimators()
#tsCmassSM.computeSnrTsz()
#tsCmassSM.computeSnrKsz()
#tsCmassSM.kszNullTests()
#tsCmassSM.plotCovKsz()

name = cmassNMariana.name + "_pactf90daynight"
tsCmassNM = ThumbStack(u, cmassNMariana, pact90Map, pact90Mask, pact90Hit, name=name, nameLong=None, save=True, nProc=nProc)
#tsCmassNM.compareKszEstimators()
#tsCmassNM.computeSnrTsz()
#tsCmassNM.computeSnrKsz()
#tsCmassNM.kszNullTests()
#tsCmassNM.plotCovKsz()

name = cmassSKendrick.name + "_pactf90daynight"
tsCmassSK = ThumbStack(u, cmassSKendrick, pact90Map, pact90Mask, pact90Hit, name=name, nameLong=None, save=True, nProc=nProc)
#tsCmassSK.compareKszEstimators()
#tsCmassSK.computeSnrTsz()
#tsCmassSK.computeSnrKsz()
#tsCmassSK.kszNullTests()
#tsCmassSK.plotCovKsz()

name = cmassNKendrick.name + "_pactf90daynight"
tsCmassNK = ThumbStack(u, cmassNKendrick, pact90Map, pact90Mask, pact90Hit, name=name, nameLong=None, save=True, nProc=nProc)
#tsCmassNK.compareKszEstimators()
#tsCmassNK.computeSnrTsz()
#tsCmassNK.computeSnrKsz()
#tsCmassNK.kszNullTests()
#tsCmassNK.plotCovKsz()

name = cmassKendrick.name + "_pactf90daynight"
tsCmassK = ThumbStack(u, cmassKendrick, pact90Map, pact90Mask, pact90Hit, name=name, nameLong=None, save=True, nProc=nProc)
#tsCmassK.compareKszEstimators()
#tsCmassK.computeSnrTsz()
#tsCmassK.computeSnrKsz()
#tsCmassK.kszNullTests()
#tsCmassK.plotCovKsz()

name = lowzSKendrick.name + "_pactf90daynight"
tsLowzSK = ThumbStack(u, lowzSKendrick, pact90Map, pact90Mask, pact90Hit, name=name, nameLong=None, save=True, nProc=nProc)
#tsLowzSK.compareKszEstimators()
#tsLowzSK.computeSnrTsz()
#tsLowzSK.computeSnrKsz()
#tsLowzSK.kszNullTests()
#tsLowzSK.plotCovKsz()

name = lowzNKendrick.name + "_pactf90daynight"
tsLowzNK = ThumbStack(u, lowzNKendrick, pact90Map, pact90Mask, pact90Hit, name=name, nameLong=None, save=True, nProc=nProc)
#tsLowzNK.compareKszEstimators()
#tsLowzNK.computeSnrTsz()
#tsLowzNK.computeSnrKsz()
#tsLowzNK.kszNullTests()
#tsLowzNK.plotCovKsz()

name = lowzKendrick.name + "_pactf90daynight"
tsLowzK = ThumbStack(u, lowzKendrick, pact90Map, pact90Mask, pact90Hit, name=name, nameLong=None, save=True, nProc=nProc)
#tsLowzK.compareKszEstimators()
#tsLowzK.computeSnrTsz()
#tsLowzK.computeSnrKsz()
#tsLowzK.kszNullTests()
#tsLowzK.plotCovKsz()

name = bossKendrick.name + "_pactf90daynight"
tsBossK = ThumbStack(u, bossKendrick, pact90Map, pact90Mask, pact90Hit, name=name, nameLong=None, save=True, nProc=nProc)
#tsBossK.compareKszEstimators()
#tsBossK.computeSnrTsz()
#tsBossK.computeSnrKsz()
#tsBossK.kszNullTests()
#tsBossK.plotCovKsz()


###################################################################################
# Plot kSZ for the various samples


# plot CMASS Mariana
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# convert from sr to arcmin^2
factor = (180.*60./np.pi)**2
#
# CMASS M
ax.errorbar(tsCmassM.RApArcmin+0.02, factor * tsCmassM.kSZ, factor * np.sqrt(np.diag(tsCmassM.covKsz)), c='r', label=r'CMASS M')
ax.errorbar(tsCmassNM.RApArcmin, factor * tsCmassNM.kSZ, factor * np.sqrt(np.diag(tsCmassNM.covKsz)), fmt='--', c='r', alpha=0.2, label=r'CMASS N M')
ax.errorbar(tsCmassSM.RApArcmin+0.01, factor * tsCmassSM.kSZ, factor * np.sqrt(np.diag(tsCmassSM.covKsz)), fmt='-.', c='r', alpha=0.2, label=r'CMASS S M')
#
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$\T_\text{kSZ}$ [$\mu K\cdot$arcmin$^2$]')
ax.set_ylim((0., 2.))
#
path = tsCmassM.pathFig+"/ksz_cmass_mariana.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()


# plot CMASS Kendrick
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# CMASS K
ax.errorbar(tsCmassK.RApArcmin+0.02, factor * tsCmassK.kSZ, factor * np.sqrt(np.diag(tsCmassK.covKsz)), c='r', label=r'CMASS K')
ax.errorbar(tsCmassNK.RApArcmin, factor * tsCmassNK.kSZ, factor * np.sqrt(np.diag(tsCmassNK.covKsz)), fmt='--', c='r', alpha=0.2, label=r'CMASS N K')
ax.errorbar(tsCmassSK.RApArcmin+0.01, factor * tsCmassSK.kSZ, factor * np.sqrt(np.diag(tsCmassSK.covKsz)), fmt='-.', c='r', alpha=0.2, label=r'CMASS S K')
#
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$\T_\text{kSZ}$ [$\mu K\cdot$arcmin$^2$]')
ax.set_ylim((0., 2.))
#
path = tsCmassK.pathFig+"/ksz_cmass_kendrick.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()


# plot LOWZ Kendrick
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# CMASS K
ax.errorbar(tsLowzK.RApArcmin+0.02, factor * tsLowzK.kSZ, factor * np.sqrt(np.diag(tsLowzK.covKsz)), c='r', label=r'LOWZ K')
ax.errorbar(tsLowzNK.RApArcmin, factor * tsLowzNK.kSZ, factor * np.sqrt(np.diag(tsLowzNK.covKsz)), fmt='--', c='r', alpha=0.2, label=r'LOWZ N K')
ax.errorbar(tsLowzSK.RApArcmin+0.01, factor * tsLowzSK.kSZ, factor * np.sqrt(np.diag(tsLowzSK.covKsz)), fmt='-.', c='r', alpha=0.2, label=r'LOWZ S K')
#
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$\T_\text{kSZ}$ [$\mu K\cdot$arcmin$^2$]')
ax.set_ylim((0., 2.))
#
path = tsLowzK.pathFig+"/ksz_lowz_kendrick.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()

