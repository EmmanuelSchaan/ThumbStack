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
plt.switch_backend('Agg')


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

# Mariana
cmassSMariana = Catalog(u, massConversion, name="cmass_s_mariana", nameLong="CMASS S M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt", save=False)
cmassNMariana = Catalog(u, massConversion, name="cmass_n_mariana", nameLong="CMASS N M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt", save=False)
cmassMariana = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False)

# Kendrick
cmassSKendrick = Catalog(u, massConversion, name="cmass_s_kendrick", nameLong="CMASS S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_S_kendrick.txt", save=False)
cmassNKendrick = Catalog(u, massConversion, name="cmass_n_kendrick", nameLong="CMASS N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_N_kendrick.txt", save=False)
cmassKendrick = Catalog(u, massConversion, name="cmass_kendrick", nameLong="CMASS K", save=False)
lowzSKendrick = Catalog(u, massConversion, name="lowz_s_kendrick", nameLong="LOWZ S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_S_kendrick.txt", save=False)
lowzNKendrick = Catalog(u, massConversion, name="lowz_n_kendrick", nameLong="LOWZ N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_N_kendrick.txt", save=False)
lowzKendrick = Catalog(u, massConversion, name="lowz_kendrick", nameLong="LOWZ K", save=False)
bossKendrick = Catalog(u, massConversion, name="boss_kendrick", nameLong="BOSS K", save=False)


###################################################################################
###################################################################################
# Read CMB maps


pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/"
#
pathMap = pathIn + "car_planck545.fits"
pathMask = pathIn + "car_planck_gal_ps_mask.fits"
#
pactMap = enmap.read_map(pathMap)
pactMask = enmap.read_map(pathMask)


###################################################################################
###################################################################################
# Stacking


import thumbstack
reload(thumbstack)
from thumbstack import *


name = cmassMariana.name + "_pactf150night20190311"
tsCmassM = ThumbStack(u, cmassMariana, pactMap, pactMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)


###################################################################################

name = cmassSMariana.name + "_pactf150night20190311"
tsCmassSM = ThumbStack(u, cmassSMariana, pactMap, pactMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)

name = cmassNMariana.name + "_pactf150night20190311"
tsCmassNM = ThumbStack(u, cmassNMariana, pactMap, pactMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)



name = cmassSKendrick.name + "_pactf150night20190311"
tsCmassSK = ThumbStack(u, cmassSKendrick, pactMap, pactMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)

name = cmassNKendrick.name + "_pactf150night20190311"
tsCmassNK = ThumbStack(u, cmassNKendrick, pactMap, pactMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)

name = cmassKendrick.name + "_pactf150night20190311"
tsCmassK = ThumbStack(u, cmassKendrick, pactMap, pactMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)

name = lowzSKendrick.name + "_pactf150night20190311"
tsLowzSK = ThumbStack(u, lowzSKendrick, pactMap, pactMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)

name = lowzNKendrick.name + "_pactf150night20190311"
tsLowzNK = ThumbStack(u, lowzNKendrick, pactMap, pactMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)

name = lowzKendrick.name + "_pactf150night20190311"
tsLowzK = ThumbStack(u, lowzKendrick, pactMap, pactMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)

name = bossKendrick.name + "_pactf150night20190311"
tsBossK = ThumbStack(u, bossKendrick, pactMap, pactMask, cmbHit=None, name=name, nameLong=None, save=True, nProc=nProc)


