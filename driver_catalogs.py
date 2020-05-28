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
'''
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
bossKendrick = Catalog(u, massConversion, name="boss_kendrick", nameLong="BOSS K", save=False)
#bossKendrick.plotHistograms()
#bossKendrick.plotFootprint()
'''

###################################################################################
# Null test: vMariana - vKendrick
# get the intersection of the catalogs
'''
cmassMariana = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False)
cmassKendrick = Catalog(u, massConversion, name="cmass_kendrick", nameLong="CMASS K", save=False)

cmassMKDiff = cmassMariana.copy(name="cmass_mk_diff", nameLong="CMASS M-K")
cmassMKDiff.intersectCatalog(cmassKendrick, vDiff=True, save=True, nProc=32)
cmassMKDiff = Catalog(u, massConversion, name="cmass_mk_diff", nameLong="CMASS M-K", save=False)
'''

###################################################################################
# mini CMASS Mariana, for debugging
'''
nObj = 50000
cmassMarianaShort = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False, nObj=nObj)
cmassMarianaShort.name = "cmass_mariana_short"
'''


###################################################################################
###################################################################################
# Read catalogs and plot

cmassMariana = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False)
cmassKendrick = Catalog(u, massConversion, name="cmass_kendrick", nameLong="CMASS K", save=False)
lowzKendrick = Catalog(u, massConversion, name="lowz_kendrick", nameLong="LOWZ K", save=False)


###################################################################################
# plot dn/dz
'''
# read overlap flags
# for PACT maps
path = "./output/thumbstack/cmass_mariana_pactf150daynight20200228/overlap_flag.txt"
overlapCmassM = np.genfromtxt(path).astype(bool)
path = "./output/thumbstack/cmass_kendrick_pactf150daynight20200228/overlap_flag.txt"
overlapCmassK = np.genfromtxt(path).astype(bool)
path = "./output/thumbstack/lowz_kendrick_pactf150daynight20200228/overlap_flag.txt"
overlapLowzK = np.genfromtxt(path).astype(bool)

print "CMASS M:", np.sum(overlapCmassM), "galaxies overlap out of", len(overlapCmassM), "ie", 1.*np.sum(overlapCmassM)/len(overlapCmassM),"%"
print "CMASS K:", np.sum(overlapCmassK), "galaxies overlap out of", len(overlapCmassK), "ie", 1.*np.sum(overlapCmassK)/len(overlapCmassK),"%"
print "LOWZ K:", np.sum(overlapLowzK), "galaxies overlap out of", len(overlapLowzK), "ie", 1.*np.sum(overlapLowzK)/len(overlapLowzK),"%"


# Bins for histograms
nBins = 51
Bins = np.linspace(0., 0.7, nBins)
binwidth = Bins[1:] - Bins[:-1]


fig = plt.figure(0)
ax = fig.add_subplot(111)
#
# CMASS Mariana
histX = np.histogram(cmassMariana.Z[overlapCmassM], Bins)[0].astype(float)
ax.bar(Bins[:-1], histX, binwidth, color='b', alpha=0.5, label=r'CMASS M ('+np.sum(overlapCmassM).astype(str)+' total)')
#
# CMASS Kendrick
histX = np.histogram(cmassKendrick.Z[overlapCmassK], Bins)[0].astype(float)
ax.bar(Bins[:-1], histX, binwidth, color='c', alpha=0.5, label=r'CMASS K ('+np.sum(overlapCmassK).astype(str)+' total)')
#
# LOWZ Kendrick
histX = np.histogram(lowzKendrick.Z[overlapLowzK], Bins)[0].astype(float)
ax.bar(Bins[:-1], histX, binwidth, color='y', alpha=0.5, label=r'LOWZ K ('+np.sum(overlapLowzK).astype(str)+' total)')
#
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlim((Bins[0], Bins[-1]))
ax.set_xlabel(r'$z$')
ax.set_ylabel(r'$N_\text{galaxies}$')
#
path = cmassMariana.pathFig + "/summary_dndz.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()
'''


###################################################################################
# Plot overlap CMASS M and ACT

'''
# AdvACT hit map to superimpose
# PACT day+night, 20200228, Planck Galactic masks 60%
pathMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits"
pathMask = "./output/cmb_map/pact20200228/" + "mask_full_foot_gal60_ps.fits"
pathHit = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits"
#
#hMap = enmap.read_map(pathMap)[0] * enmap.read_map(pathMask)
hMap = enmap.read_map(pathHit)[0] 
hMap *= enmap.read_map(pathMask)[0]
hMap = np.log(np.abs(hMap)+1.e-5)
hMap = enmap.to_healpix(hMap)

cmassMariana.plotFootprint(hMap=hMap)
'''

###################################################################################
# Generate vtk files for visualization

'''
cmassMariana.writeVtk()
cmassKendrick.writeVtk()
lowzKendrick.writeVtk()
cmassMKDiff.writeVtk()
'''

