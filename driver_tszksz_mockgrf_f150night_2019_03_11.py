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
'''
# CMASS
cmassSMariana = Catalog(u, massConversion, name="cmass_s_mariana", nameLong="CMASS S M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt", save=False)
cmassNMariana = Catalog(u, massConversion, name="cmass_n_mariana", nameLong="CMASS N M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt", save=False)
# combined catalog
#cmassMariana = cmassSMariana.copy(name="cmass_mariana", nameLong="CMASS M")
#cmassMariana.addCatalog(cmassNMariana, save=True)
cmassMariana = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False)
'''

###################################################################################
# Kendrick

# CMASS
cmassSKendrick = Catalog(u, massConversion, name="cmass_s_kendrick", nameLong="CMASS S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_S_kendrick.txt", save=False)
cmassNKendrick = Catalog(u, massConversion, name="cmass_n_kendrick", nameLong="CMASS N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_N_kendrick.txt", save=False)
# combined catalog
#cmassKendrick = cmassSKendrick.copy(name="cmass_kendrick", nameLong="CMASS K")
#cmassKendrick.addCatalog(cmassNKendrick, save=True)
cmassKendrick = Catalog(u, massConversion, name="cmass_kendrick", nameLong="CMASS K", save=False)
'''
# LOWZ
lowzSKendrick = Catalog(u, massConversion, name="lowz_s_kendrick", nameLong="LOWZ S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_S_kendrick.txt", save=False)
lowzNKendrick = Catalog(u, massConversion, name="lowz_n_kendrick", nameLong="LOWZ N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_N_kendrick.txt", save=False)
# combined catalog
#lowzKendrick = lowzSKendrick.copy(name="lowz_kendrick", nameLong="LOWZ K")
#lowzKendrick.addCatalog(lowzNKendrick, save=True)
lowzKendrick = Catalog(u, massConversion, name="lowz_kendrick", nameLong="LOWZ K", save=False)

# BOSS = CMASS + LOWZ
#bossKendrick = cmassKendrick.copy(name="boss_kendrick", nameLong="BOSS K")
#bossKendrick.addCatalog(lowzKendrick, save=True)
bossKendrick = Catalog(u, massConversion, name="boss_kendrick", nameLong="BOSS K", save=False)
'''


###################################################################################
###################################################################################
# Read CMB maps

# directory of GRF mocks
pathGRF = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/output/cmb_map/mocks_grf_planck_act_coadd_2019_03_11/"
# path to mean power spectrum of GRF mocks
pathPower = pathGRF + "mean_cl.txt"


# path to true hit count map and mask: Planck + ACT 150GHz day and night
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/"
pathHit = pathIn + "act_planck_f150_prelim_div_mono.fits"
pathMask = pathIn + "f150_mask_foot_planck_ps_car.fits"

# read maps in common for all mocks
pactMask = enmap.read_map(pathMask)
pactHit = enmap.read_map(pathHit)

## measured power spectrum
#data = np.genfromtxt(pathPower)  # l, Cl, sCl
#data = np.nan_to_num(data)
#fCl = interp1d(data[:,0], data[:,1], kind='linear', bounds_error=False, fill_value=0.)

# theory power spectrum
cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)


###################################################################################


import thumbstack
reload(thumbstack)
from thumbstack import *


def analyzeMock(iMock):
   print "- analyzing mock", iMock

   # path to each GRF mock map
   pathMap = pathGRF + "mock_"+str(iMock)+"_grf_f150_daynight.fits"

   # Read mock map
   pactMap = enmap.read_map(pathMap)#[0]   # keep only temperature

   # Stacking
   name = cmassSKendrick.name + "_pactf150night20190311_mock"+str(iMock)
   tsCmassSK = ThumbStack(u, cmassSKendrick, pactMap, pactMask, pactHit, name=name, nameLong=None, save=True, nProc=nProc)

   return


###################################################################################
# Analyze all the mocks

# set of mocks
nMocks = 1000
IMocks = range(100)
#IMocks = [0]

# run analysis serially
tStart = time()
map(analyzeMock, IMocks)
tStop = time()
print "took", (tStop-tStart)/60., "min"


###################################################################################
