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
lowzSKendrick = Catalog(u, massConversion, name="lowz_s_kendrick", nameLong="LOWZ S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_S_kendrick.txt", save=True)
#lowzSKendrick.plotHistograms()
#lowzSKendrick.plotFootprint()
#
lowzNKendrick = Catalog(u, massConversion, name="lowz_n_kendrick", nameLong="LOWZ N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_N_kendrick.txt", save=True)
#lowzNKendrick.plotHistograms()
#lowzNKendrick.plotFootprint()
#
# combined catalog
lowzKendrick = lowzSKendrick.copy(name="lowz_kendrick", nameLong="LOWZ K")
lowzKendrick.addCatalog(lowzNKendrick, save=True)
#lowzKendrick.plotHistograms()
#lowzKendrick.plotFootprint()

# BOSS = CMASS + LOWZ
bossKendrick = cmassSKendrick.copy(name="boss_kendrick", nameLong="BOSS K")
bossKendrick.addCatalog(cmassNKendrick, save=True)
bossKendrick.addCatalog(lowzSKendrick, save=True)
bossKendrick.addCatalog(lowzNKendrick, save=True)
#bossKendrick.plotHistograms()
#bossKendrick.plotFootprint()



###################################################################################
###################################################################################
# Read CMB maps

# Planck + ACT 150GHz day and night
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/"
pathMap = pathIn + "f150_daynight_all_map_mono.fits"
pathHit = pathIn + "f150_daynight_all_div_mono.fits"
pathMask = pathIn + "f150_mask_foot_planck_ps_car.fits"
pathPower = pathIn + "f150_power_T_masked.txt"

tStart = time()
print "- Read CMB map, mask and hit count"
pact150Map = enmap.read_map(pathMap)[0]   # keep only temperature
pact150Mask = enmap.read_map(pathMask)
pact150Hit = enmap.read_map(pathHit)

print "Set up interpolations"
#!!!! NO! This prefiltering seems to mess up the map values, in a complicated way...
# This pre-filtering step introduces some ringing in the maps

#utils.interpol_prefilter(pact150Map, inplace=True)
#utils.interpol_prefilter(pact150Mask, inplace=True)   # order=2 seems to reduce ringing at sharp transitions
#utils.interpol_prefilter(pact150Hit, inplace=True)

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

'''
name = cmassMariana.name + "_pactf150daynight"
ts = ThumbStack(u, cmassMariana, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)
'''

















#mask = ts.catalogMask(overlap=True, psMask=True)
#ts.analyzeObject(0, test=True)
#ts.examineCmbMaps()

# Expected std dev of AP filter, function of disk radius in rad
#fsApActual = lambda r0: cmb1_4.fsigmaDiskRing(r0, thetaIn=None, thetaOut=None, fCl=fCl, lMin=1., lMax=1.e5)
#fsApNoiseless = lambda r0: cmb1_4.fsigmaDiskRing(r0, thetaIn=None, thetaOut=None, fCl=cmb1_4.flensedTT, lMin=1., lMax=1.e5)
#ts.examineHistograms(fsAp=[fsApActual, fsApNoiseless])

#ts.measureVarFromHitCount(plot=True)
#ts.compareKszEstimators()

#cov = ts.kszCovBootstrap(nSamples=1000, nProc=nProc)
#print np.sqrt(np.diag(cov))

#ts.computeSnrKsz()

#ts.kszNullTests()







###################################################################################

name = cmassSMariana.name + "_pactf150daynight"
ts = ThumbStack(u, cmassSMariana, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)

name = cmassNMariana.name + "_pactf150daynight"
ts = ThumbStack(u, cmassNMariana, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)

name = cmassSKendrick.name + "_pactf150daynight"
ts = ThumbStack(u, cmassSKendrick, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)

name = cmassNKendrick.name + "_pactf150daynight"
ts = ThumbStack(u, cmassNKendrick, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)

name = cmassKendrick.name + "_pactf150daynight"
ts = ThumbStack(u, cmassKendrick, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)

name = lowzSKendrick.name + "_pactf150daynight"
ts = ThumbStack(u, lowzSKendrick, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)

name = lowzNKendrick.name + "_pactf150daynight"
ts = ThumbStack(u, lowzNKendrick, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)

name = lowzKendrick.name + "_pactf150daynight"
ts = ThumbStack(u, lowzKendrick, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)

name = bossKendrick.name + "_pactf150daynight"
ts = ThumbStack(u, bossKendrick, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)
