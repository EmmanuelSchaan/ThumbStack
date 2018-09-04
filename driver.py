import universe
reload(universe)
from universe import *

import mass_conversion
reload(mass_conversion)
from mass_conversion import *

#import catalog
#reload(catalog)
#from catalog import *
#
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
# Loading catalog




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


