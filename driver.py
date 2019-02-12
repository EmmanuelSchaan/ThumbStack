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

#import cmb_map
#reload(cmb_map)
#from cmb_map import *
#
#import diskring_filter
#reload(diskring_filter)
#from diskring_filter import *
#
#import pixel_noise
#reload(pixel_noise)
#from pixel_noise import *
#
#import matched_filter
#reload(matched_filter)
#from matched_filter import *
#
#import posterior_alpha
#reload(posterior_alpha)
#from posterior_alpha import *



# running on cori
# 68 cores per knl node, 32 cores per haswell node
#salloc -N 3 --qos=interactive -C haswell -t 04:00:00 -L SCRATCH


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
# Loading catalogs

###################################################################################
# Mariana

# CMASS
cmassSMariana = Catalog(u, massConversion, name="cmass_s_mariana", nameLong="CMASS S M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt", save=False)
#cmassSMariana.plotHistograms()
#cmassSMariana.plotFootprint()
#
'''
cmassNMariana = Catalog(u, massConversion, name="cmass_n_mariana", nameLong="CMASS N M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt", save=False)
#cmassNMariana.plotHistograms()
#cmassNMariana.plotFootprint()
#
# combined catalog
cmassMariana = cmassSMariana.copy(name="cmass_mariana", nameLong="CMASS M")
cmassMariana.addCatalog(cmassNMariana, save=True)
#cmassMariana.plotHistograms()
#cmassMariana.plotFootprint()
'''

###################################################################################
# Kendrick
'''
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
cmassKendrick.addCatalog(cmassNKendrick, save=True)
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
'''
###################################################################################

#
## Mariana
#cmassSMariana = CMASS_S_Mariana(u, massConversion, save=False)
##cmassSMariana.plotHistograms()
##cmassSMariana.plotFootprint()
##
#cmassNMariana = CMASS_N_Mariana(u, massConversion, save=False)
##cmassNMariana.plotHistograms()
##cmassNMariana.plotFootprint()
##
## combined catalog
#cmassMariana = CMASS_S_Mariana(u, massConversion, save=False)
#cmassMariana.addCatalog(cmassNMariana)
##cmassMariana.plotHistograms()
##cmassMariana.plotFootprint()
#
#
## Kendrick
#
## CMASS
#cmassSKendrick = CMASS_S_Kendrick(u, massConversion, save=False)
##cmassSKendrick.plotHistograms()
##cmassSKendrick.plotFootprint()
##
#cmassNKendrick = CMASS_N_Kendrick(u, massConversion, save=False)
##cmassNKendrick.plotHistograms()
##cmassNKendrick.plotFootprint()
##
## combined catalog
#cmassKendrick = CMASS_S_Kendrick(u, massConversion, save=False)
#cmassKendrick.addCatalog(cmassNKendrick)
##cmassKendrick.plotHistograms()
##cmassKendrick.plotFootprint()
#
## LOWZ
#lowzSKendrick = LOWZ_S_Kendrick(u, massConversion, save=False)
##lowzSKendrick.plotHistograms()
##lowzSKendrick.plotFootprint()
##
#lowzNKendrick = LOWZ_N_Kendrick(u, massConversion, save=False)
##lowzNKendrick.plotHistograms()
##lowzNKendrick.plotFootprint()
##
## combined catalog
#lowzKendrick = LOWZ_S_Kendrick(u, massConversion, save=False)
#lowzKendrick.addCatalog(lowzNKendrick)
##lowzKendrick.plotHistograms()
##lowzKendrick.plotFootprint()
#
## BOSS = CMASS + LOWZ
#bossKendrick = CMASS_S_Kendrick(u, massConversion, save=False)
#bossKendrick.addCatalog(cmassNKendrick)
#bossKendrick.addCatalog(lowzSKendrick)
#bossKendrick.addCatalog(lowzNKendrick)
##lowzKendrick.plotHistograms()
##lowzKendrick.plotFootprint()
#
#
## Other functions coded up
##catalog.printProperties()
##catalog.compareV1dRms()
#


###################################################################################
###################################################################################
# Read CMB maps

# Planck + ACT 150GHz day and night
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/"
pathMap = pathIn + "f150_daynight_all_map_mono.fits"
pathHit = pathIn + "f150_daynight_all_div_mono.fits"
pathMask = pathIn + "f150_mask_foot_planck_ps_car.fits"

tStart = time()
print "- Read CMB map, mask and hit count"
pact150Map = enmap.read_map(pathMap)
pact150Mask = enmap.read_map(pathMask)
pact150Hit = enmap.read_map(pathHit)
print "Set up interpolations"
#pact150Map = utils.interpol_prefilter(pact150Map, inplace=True)
#pact150Mask = utils.interpol_prefilter(pact150Mask, inplace=True)
#pact150Hit = utils.interpol_prefilter(pact150Hit, inplace=True)

utils.interpol_prefilter(pact150Map, inplace=True)
utils.interpol_prefilter(pact150Mask, inplace=True)
utils.interpol_prefilter(pact150Hit, inplace=True)

tStop = time()
print "took", (tStop-tStart)/60., "min"


###################################################################################
###################################################################################
# Stacking


import thumbstack
reload(thumbstack)
from thumbstack import *



name = cmassSMariana.name + "_pactf150daynight"
ts = ThumbStack(u, cmassSMariana, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=True, nProc=nProc)


#ts = ThumbStack(u, cmassSMariana, pathMap=pathMap, pathMask=pathMask, pathHit=pathHit, name=name, nameLong=None, save=False, nProc=nProc)




#ts.examineHistograms()





'''
mask = ts.catalogMask(overlap=True, mVir=None, extraSelection=1.)



ts.filtMap[np.where(np.isfinite(ts.filtMap)==False)] = 0.


      # mask highest masses?
      Y = Y[~mask]
      TSZ = TSZ[~mask]
      KSZ = KSZ[~mask]
      Nobj = len(KSZ)




ts.medianFiltMap = np.median(ts.filtMap, axis=0)






fig=plt.figure(0)
ax=fig.subplot(111)
#
ax.plot(ts.RApMpch, ts.medianFiltMap)
'''







'''
##################################################################################
# theoretical pixel noise power spectrum

pixelNoise = ACT148PixelNoise()
#pixelNoise = PlanckSMICAPixelNoise()

###################################################################################
###################################################################################
# Mariana/Kendrick velocities

#cmbMap = ACTDeep5656()
#cmbMap = ACTDeep56c7v5Jon()
#cmbMap = ACTDeep56c7v5Sigurd()
#cmbMap = ACTEquatorial(part='left')
cmbMap = ACTCoadd()

# CMASS South
#catalog = CMASSSouthSmoo(U, massConversion, meanMass=False, smoo=10, rand=False, save=False)
catalog = CMASSSouthMarianaKendrick(U, massConversion, meanMass=False, save=False)


#catalog = CMASSSouthDR10(U, massConversion, save=False)
#catalog = CMASSSouth(U, massConversion, save=False)
#catalog = CMASSSouthSmooMeanMassFull(U, massConversion, save=False, smoo=10)
#catalog = CMASSSouthDR10Kendrick(U, massConversion, meanMass=True, save=False)
#
# nulls
#catalog = CMASSSouthMarianaKendrick(U, massConversion, shift=True, raShift=-10., decShift=-10., parity=0, save=False)

# LOWZ South
#catalog = LOWZSouthDR10Kendrick(U, massConversion)
#catalog = LOWZSouthMariana(U, massConversion, save=False)
'''

###################################################################################
# BOSS North
'''
cmbMap = ACTBOSSNorth()

#catalog = CMASSNorthDR10Kendrick(U, massConversion, save=False)
catalog = CMASSNorthSmoo(U, massConversion, save=False, smoo=10)
'''

##################################################################################
##################################################################################
# disk minus ring filter
'''
diskRingFilter = DiskRingFilter(U, catalog, cmbMap, save=False)

# disk-ring
#diskRingFilter.fcheckMeanFilterAll()
#diskRingFilter.fbinMassAll()
#diskRingFilter.fmeasureGasFractionAll(mMin=1.e6, mMax=2.e14)
#diskRingFilter.fmeasureSNRmMaxAll()

#diskRingFilter.fplotGasFraction(mMax=1.e14, imin=1)
'''


