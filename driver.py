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







import cmb
reload(cmb)
from cmb import *


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


#name = cmassSMariana.name + "_pactf150daynight"
#tsS = ThumbStack(u, cmassSMariana, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)
#
#name = cmassNMariana.name + "_pactf150daynight"
#tsN = ThumbStack(u, cmassNMariana, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)

name = cmassMariana.name + "_pactf150daynight"
ts = ThumbStack(u, cmassMariana, pact150Map, pact150Mask, pact150Hit, name=name, nameLong=None, save=False, nProc=nProc)



#mask = ts.catalogMask(overlap=True, psMask=True)

#ts.analyzeObject(0, test=True)

#ts.examineCmbMaps()

# Expected std dev of AP filter, function of disk radius in rad
#fsApActual = lambda r0: cmb1_4.fsigmaDiskRing(r0, thetaIn=None, thetaOut=None, fCl=fCl, lMin=1., lMax=1.e5)
#fsApNoiseless = lambda r0: cmb1_4.fsigmaDiskRing(r0, thetaIn=None, thetaOut=None, fCl=cmb1_4.flensedTT, lMin=1., lMax=1.e5)
#ts.examineHistograms(fsAp=[fsApActual, fsApNoiseless])

#ts.measureVarFromHitCount(plot=True)
#ts.measureKSZ()

cov = ts.kszCovBootstrap(nSamples=1000, nProc=nProc)
print np.sqrt(np.diag(cov))



'''
iObj = 0
ra = ts.Catalog.RA[iObj]
dec = ts.Catalog.DEC[iObj]

opos, stampMap, stampMask, stampHit = ts.extractStamp(ra, dec, dxDeg=0.25, dyDeg=0.25, resArcmin=0.25, proj='cea', test=False)

plots = enplot.plot(stampMap, grid=True)


path = ts.pathTestFig+"/test.pdf"
#plots.write(path, plots)
enplot.write(path,plots)


'''









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


