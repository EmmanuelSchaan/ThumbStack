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
cmassSMariana = Catalog(u, massConversion, name="cmass_s_mariana", nameLong="CMASS S M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt", save=False)
cmassNMariana = Catalog(u, massConversion, name="cmass_n_mariana", nameLong="CMASS N M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt", save=False)
# combined catalog
cmassMariana = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False)
#cmassMariana.plotHistograms()
#cmassMariana.plotFootprint()
#cmassMariana.printProperties()


###################################################################################
# Kendrick

# CMASS
cmassSKendrick = Catalog(u, massConversion, name="cmass_s_kendrick", nameLong="CMASS S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_S_kendrick.txt", save=False)
cmassNKendrick = Catalog(u, massConversion, name="cmass_n_kendrick", nameLong="CMASS N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_N_kendrick.txt", save=False)
# combined catalog
cmassKendrick = Catalog(u, massConversion, name="cmass_kendrick", nameLong="CMASS K", save=False)

# LOWZ
lowzSKendrick = Catalog(u, massConversion, name="lowz_s_kendrick", nameLong="LOWZ S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_S_kendrick.txt", save=False)
lowzNKendrick = Catalog(u, massConversion, name="lowz_n_kendrick", nameLong="LOWZ N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_N_kendrick.txt", save=False)
# combined catalog
lowzKendrick = Catalog(u, massConversion, name="lowz_kendrick", nameLong="LOWZ K", save=False)

# BOSS = CMASS + LOWZ
bossKendrick = Catalog(u, massConversion, name="boss_kendrick", nameLong="BOSS K", save=False)


###################################################################################
###################################################################################
# plot dn/dz

# read overlap flags
# for PACT maps
path = "./output/thumbstack/cmass_mariana_pactf150daynight20200228/overlap_flag.txt"
overlapCmassM = np.genfromtxt(path).astype(bool)
path = "./output/thumbstack/cmass_kendrick_pactf150daynight20200228/overlap_flag.txt"
overlapCmassK = np.genfromtxt(path).astype(bool)
path = "./output/thumbstack/lowz_kendrick_pactf150daynight20200228/overlap_flag.txt"
overlapLowzK = np.genfromtxt(path).astype(bool)

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
