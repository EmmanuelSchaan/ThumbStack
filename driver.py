import universe
reload(universe)
from universe import *

import mass_conversion
reload(mass_conversion)
from mass_conversion import *

import catalog
reload(catalog)
from catalog import *

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

##################################################################################
# cosmological parameters

u = UnivMariana()

###################################################################################
# M*-Mh relation

massConversion = MassConversionKravtsov14()
#massConversion.plot()

###################################################################################
# Loading catalogs

# Mariana
cmassSMariana = CMASS_S_Mariana(u, massConversion, save=False)
#cmassSMariana.plotHistograms()
#cmassSMariana.plotFootprint()
#
cmassNMariana = CMASS_N_Mariana(u, massConversion, save=False)
#cmassNMariana.plotHistograms()
#cmassNMariana.plotFootprint()
#
# combined catalog
cmassMariana = CMASS_S_Mariana(u, massConversion, save=False)
cmassMariana.addCatalog(cmassNMariana)
#cmassMariana.plotHistograms()
#cmassMariana.plotFootprint()


# Kendrick

# CMASS
cmassSKendrick = CMASS_S_Kendrick(u, massConversion, save=False)
#cmassSKendrick.plotHistograms()
#cmassSKendrick.plotFootprint()
#
cmassNKendrick = CMASS_N_Kendrick(u, massConversion, save=False)
#cmassNKendrick.plotHistograms()
#cmassNKendrick.plotFootprint()
#
# combined catalog
cmassKendrick = CMASS_S_Kendrick(u, massConversion, save=False)
cmassKendrick.addCatalog(cmassNKendrick)
#cmassKendrick.plotHistograms()
#cmassKendrick.plotFootprint()

# LOWZ
lowzSKendrick = LOWZ_S_Kendrick(u, massConversion, save=False)
#lowzSKendrick.plotHistograms()
#lowzSKendrick.plotFootprint()
#
lowzNKendrick = LOWZ_N_Kendrick(u, massConversion, save=False)
#lowzNKendrick.plotHistograms()
#lowzNKendrick.plotFootprint()
#
# combined catalog
lowzKendrick = LOWZ_S_Kendrick(u, massConversion, save=False)
lowzKendrick.addCatalog(lowzNKendrick)
#lowzKendrick.plotHistograms()
#lowzKendrick.plotFootprint()

# BOSS = CMASS + LOWZ
bossKendrick = CMASS_S_Kendrick(u, massConversion, save=False)
bossKendrick.addCatalog(cmassNKendrick)
bossKendrick.addCatalog(lowzSKendrick)
bossKendrick.addCatalog(lowzNKendrick)
#lowzKendrick.plotHistograms()
#lowzKendrick.plotFootprint()


# Other functions coded up
#catalog.printProperties()
#catalog.compareV1dRms()






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


