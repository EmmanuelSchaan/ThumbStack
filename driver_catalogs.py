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
# Galaxy catalogs

###################################################################################
# Mariana

# CMASS
cmassSMariana = Catalog(u, massConversion, name="cmass_s_mariana", nameLong="CMASS S M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt", save=True)
#cmassSMariana.plotHistograms()
#cmassSMariana.plotFootprint()
#cmassMariana.printProperties()
#
cmassNMariana = Catalog(u, massConversion, name="cmass_n_mariana", nameLong="CMASS N M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt", save=True)
#cmassNMariana.plotHistograms()
#cmassNMariana.plotFootprint()
#cmassMariana.printProperties()
#
# combined catalog
cmassMariana = cmassSMariana.copy(name="cmass_mariana", nameLong="CMASS M")
cmassMariana.addCatalog(cmassNMariana, save=True)
cmassMariana = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False)
#cmassMariana.plotHistograms()
#cmassMariana.plotFootprint()
#cmassMariana.printProperties()


###################################################################################
# Kendrick

# CMASS
cmassSKendrick = Catalog(u, massConversion, name="cmass_s_kendrick", nameLong="CMASS S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_S_kendrick.txt", save=True)
#cmassSKendrick.plotHistograms()
#cmassSKendrick.plotFootprint()
#
cmassNKendrick = Catalog(u, massConversion, name="cmass_n_kendrick", nameLong="CMASS N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_N_kendrick.txt", save=True)
#cmassNKendrick.plotHistograms()
#cmassNKendrick.plotFootprint()
#
# combined catalog
cmassKendrick = cmassSKendrick.copy(name="cmass_kendrick", nameLong="CMASS K")
cmassKendrick.addCatalog(cmassNKendrick, save=True)
cmassKendrick = Catalog(u, massConversion, name="cmass_kendrick", nameLong="CMASS K", save=False)
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
lowzKendrick = Catalog(u, massConversion, name="lowz_kendrick", nameLong="LOWZ K", save=False)
#lowzKendrick.plotHistograms()
#lowzKendrick.plotFootprint()

# BOSS = CMASS + LOWZ
bossKendrick = cmassSKendrick.copy(name="boss_kendrick", nameLong="BOSS K")
bossKendrick.addCatalog(cmassNKendrick, save=True)
bossKendrick.addCatalog(lowzSKendrick, save=True)
bossKendrick.addCatalog(lowzNKendrick, save=True)
bossKendrick = Catalog(u, massConversion, name="boss_kendrick", nameLong="BOSS K", save=True)
#bossKendrick.plotHistograms()
#bossKendrick.plotFootprint()



