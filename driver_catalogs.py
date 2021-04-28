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

# for cori
#plt.switch_backend('Agg')


##################################################################################

nProc = 32 # 1 haswell node on cori

##################################################################################
##################################################################################
# cosmological parameters
u = UnivMariana()

# M*-Mh relation
massConversion = MassConversionKravtsov14()
#massConversion.plot()

###################################################################################
###################################################################################
# Prepare galaxy catalogs
'''
# LOWZ full sample
lowzFull = Catalog(u, massConversion, name="lowz20200908full", nameLong="LOWZ Full", pathInCatalog="../../data/lowz/output/lowz_z0p2_0p35.txt", save=True)
lowzFull.plotHistograms()
lowzFull.plotFootprint()
lowzFull.printProperties()

# stellar mass bin 0
lowzMBin0 = Catalog(u, massConversion, name="lowz20200908mbin0", nameLong="LOWZ Bin 0", pathInCatalog="../../data/lowz/output/lowz_z0p2_0p35_mbin0.txt", save=True)
lowzMBin0.plotHistograms()
lowzMBin0.plotFootprint()
lowzMBin0.printProperties()

# stellar mass bin 1
lowzMBin1 = Catalog(u, massConversion, name="lowz20200908mbin1", nameLong="LOWZ Bin 1", pathInCatalog="../../data/lowz/output/lowz_z0p2_0p35_mbin1.txt", save=True)
lowzMBin1.plotHistograms()
lowzMBin1.plotFootprint()
lowzMBin1.printProperties()

# stellar mass bin 2
lowzMBin2 = Catalog(u, massConversion, name="lowz20200908mbin2", nameLong="LOWZ Bin 2", pathInCatalog="../../data/lowz/output/lowz_z0p2_0p35_mbin2.txt", save=True)
lowzMBin2.plotHistograms()
lowzMBin2.plotFootprint()
lowzMBin2.printProperties()
'''


###################################################################################
###################################################################################
# Read catalogs and plot

# LOWZ, full and mass bins
lowzFull = Catalog(u, massConversion, name="lowz20200908full", nameLong="LOWZ Full", save=False)
lowzMBin0 = Catalog(u, massConversion, name="lowz20200908mbin0", nameLong="LOWZ Bin 0", save=False)
lowzMBin1 = Catalog(u, massConversion, name="lowz20200908mbin1", nameLong="LOWZ Bin 1", save=False)
lowzMBin2 = Catalog(u, massConversion, name="lowz20200908mbin2", nameLong="LOWZ Bin 2", save=False)


# mini LOWZ, for debugging
nObj = 5000
lowzShort = Catalog(u, massConversion, name="lowzfull20200908", nameLong="LOWZ Full", save=False, nObj=nObj)
lowzShort.name = "lowzshort20200908"
lowzShort.nameLong = "LOWZ Short"


###################################################################################







