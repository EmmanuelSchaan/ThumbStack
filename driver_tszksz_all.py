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

# M*-Mh relation
massConversion = MassConversionKravtsov14()
#massConversion.plot()


###################################################################################
# Galaxy catalogs

print("Read galaxy catalogs")
tStart = time()

catalogs = {
      "cmass_s_mariana": Catalog(u, massConversion, name="cmass_s_mariana", nameLong="CMASS S M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt", save=False),
      "cmass_s_mariana": Catalog(u, massConversion, name="cmass_n_mariana", nameLong="CMASS N M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt", save=False),
      "cmass_mariana": Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False),
      #
      "cmass_s_kendrick": Catalog(u, massConversion, name="cmass_s_kendrick", nameLong="CMASS S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_S_kendrick.txt", save=False),
      "cmass_n_kendrick": Catalog(u, massConversion, name="cmass_n_kendrick", nameLong="CMASS N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_N_kendrick.txt", save=False),
      "cmass_kendrick": Catalog(u, massConversion, name="cmass_kendrick", nameLong="CMASS K", save=False),
      "lowz_s_kendrick": Catalog(u, massConversion, name="lowz_s_kendrick", nameLong="LOWZ S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_S_kendrick.txt", save=False),
      "lowz_n_kendrick": Catalog(u, massConversion, name="lowz_n_kendrick", nameLong="LOWZ N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_N_kendrick.txt", save=False),
      "lowz_kendrick": Catalog(u, massConversion, name="lowz_kendrick", nameLong="LOWZ K", save=False),
      "boss_kendrick": Catalog(u, massConversion, name="boss_kendrick", nameLong="BOSS K", save=False)
      }

tStop = time()
print("took "+str(round((tStop-tStart)/60., 2))+" min")


###################################################################################
# Read CMB maps

class cmbMap(object):

   def __init__(self, pathMap, pathMask, pathHit=None,  name="test"):
      self.name = name
      self.pathMap = pathMap
      self.pathMask = pathMask
      self.pathHit = pathHit

   def map(self):
      result = enmap.read_map(self.pathMap)
      # if the map contains polarization, keep only temperature
      if len(result.shape)>2:
         result = result[0]
      return result

   def mask(self):
      return enmap.read_map(self.pathMask)
   
   def hit(self):
      if self.pathHit is None:
         return None
      else:
         return enmap.read_map(self.pathHit)


print("Read CMB maps")
tStart = time()

cmbMaps = {
      "pactf150night20190311": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "act_planck_f150_prelim_map_mono.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "act_planck_f150_prelim_div_mono.fits",
         name="pactf150night20190311"),
      "pactf090night20190311": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "act_planck_f090_prelim_map_mono.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f090_mask_foot_planck_ps_car.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "act_planck_f090_prelim_div_mono.fits",
         name="pactf090night20190311"),
      #
      "plancksmica18": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck_smicacmb.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck_nonact_gal_ps_mask.fits",
         name="plancksmica18"),
      "plancksmicanosz18": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck_smicacmbnosz.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck_nonact_gal_ps_mask.fits",
         name="plancksmicanosz18"),
      "planck54518": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck545.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck_gal_ps_mask.fits",
         name="planck54518"),
      #
      "tilecpactcmbkszd56": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_deep56/" + "tilec_single_tile_deep56_cmb_map_v1.0.0_rc_joint.fits",
         "./output/cmb_map/tilec_pact_cmbksz_d56/" + "mask_foot_planck_ps_car.fits",
         name="tilecpactcmbkszd56"),
      # !!! missing y and y-no-CIB for D56
      "tilecpactcmbkszbossn": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_boss/" + "tilec_single_tile_boss_cmb_map_v1.0.0_rc_joint.fits",
         "./output/cmb_map/tilec_pact_cmbksz_boss/" + "mask_foot_planck_ps_car.fits",
         name="tilecpactcmbkszbossn"),
      "tilecpactybossn": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_boss/" + "tilec_single_tile_boss_comptony_map_v1.0.0_rc_joint.fits",
         "./output/cmb_map/tilec_pact_cmbksz_boss/" + "mask_foot_planck_ps_car.fits",
         name="tilecpactybossn"),
      "tilecpactynocibbossn": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_boss/" + "tilec_single_tile_boss_comptony_deprojects_cib_map_v1.0.0_rc_joint.fits",
         "./output/cmb_map/tilec_pact_cmbksz_boss/" + "mask_foot_planck_ps_car.fits",
         name="tilecpactynocibbossn")
      }

tStop = time()
print("took "+str(round((tStop-tStart)/60., 2))+" min")


###################################################################################
###################################################################################
# Stacking


import thumbstack
reload(thumbstack)
from thumbstack import *


save = True


#ts = {}
for cmbMapKey in cmbMaps.keys():
   cmbMap = cmbMaps[cmbMapKey].map()
   cmbMask = cmbMaps[cmbMapKey].mask()
   cmbHit = cmbMaps[cmbMapKey].hit()

   for catalogKey in catalogs.keys():
      catalog = catalogs[catalogKey]
      name = catalogs[catalogKey].name + "_" + cmbMaps[cmbMapKey].name
      #ts[cmbMapKey, catalogKey] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc)
      ts = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc)


###################################################################################












# Plot kSZ for the various samples

'''
# plot CMASS Mariana
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# convert from sr to arcmin^2
factor = (180.*60./np.pi)**2
#
# CMASS M
ax.errorbar(tsCmassM.RApArcmin+0.02, factor * tsCmassM.kSZ, factor * np.sqrt(np.diag(tsCmassM.covKsz)), c='r', label=r'CMASS M')
ax.errorbar(tsCmassNM.RApArcmin, factor * tsCmassNM.kSZ, factor * np.sqrt(np.diag(tsCmassNM.covKsz)), fmt='--', c='r', alpha=0.2, label=r'CMASS N M')
ax.errorbar(tsCmassSM.RApArcmin+0.01, factor * tsCmassSM.kSZ, factor * np.sqrt(np.diag(tsCmassSM.covKsz)), fmt='-.', c='r', alpha=0.2, label=r'CMASS S M')
#
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
#ax.set_ylim((0., 2.))
#
path = tsCmassM.pathFig+"/ksz_cmass_mariana.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()


# plot CMASS Kendrick
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# CMASS K
ax.errorbar(tsCmassK.RApArcmin+0.02, factor * tsCmassK.kSZ, factor * np.sqrt(np.diag(tsCmassK.covKsz)), c='r', label=r'CMASS K')
ax.errorbar(tsCmassNK.RApArcmin, factor * tsCmassNK.kSZ, factor * np.sqrt(np.diag(tsCmassNK.covKsz)), fmt='--', c='r', alpha=0.2, label=r'CMASS N K')
ax.errorbar(tsCmassSK.RApArcmin+0.01, factor * tsCmassSK.kSZ, factor * np.sqrt(np.diag(tsCmassSK.covKsz)), fmt='-.', c='r', alpha=0.2, label=r'CMASS S K')
#
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
ax.set_ylim((0., 10.))
#
path = tsCmassK.pathFig+"/ksz_cmass_kendrick.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()


# plot LOWZ Kendrick
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.errorbar(tsLowzK.RApArcmin+0.02, factor * tsLowzK.kSZ, factor * np.sqrt(np.diag(tsLowzK.covKsz)), c='r', label=r'LOWZ K')
ax.errorbar(tsLowzNK.RApArcmin, factor * tsLowzNK.kSZ, factor * np.sqrt(np.diag(tsLowzNK.covKsz)), fmt='--', c='r', alpha=0.2, label=r'LOWZ N K')
ax.errorbar(tsLowzSK.RApArcmin+0.01, factor * tsLowzSK.kSZ, factor * np.sqrt(np.diag(tsLowzSK.covKsz)), fmt='-.', c='r', alpha=0.2, label=r'LOWZ S K')
#
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
#ax.set_ylim((0., 2.))
#
path = tsLowzK.pathFig+"/ksz_lowz_kendrick.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()



# plot BOSS Kendrick
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# CMASS K
ax.errorbar(tsBossK.RApArcmin, factor * tsBossK.kSZ, factor * np.sqrt(np.diag(tsBossK.covKsz)), fmt='--', c='r', label=r'BOSS K')
#
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
#ax.set_xlim((0., 2.))
#
path = tsBossK.pathFig+"/ksz_kendrick.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()



###################################################################################
# Plot tSZ for the various samples


# plot CMASS Mariana
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# convert from sr to arcmin^2
factor = (180.*60./np.pi)**2
#
# CMASS M
ax.errorbar(tsCmassM.RApArcmin+0.02, factor * tsCmassM.tSZ, factor * np.sqrt(np.diag(tsCmassM.covTsz)), c='r', label=r'CMASS M')
ax.errorbar(tsCmassNM.RApArcmin, factor * tsCmassNM.tSZ, factor * np.sqrt(np.diag(tsCmassNM.covTsz)), fmt='--', c='r', alpha=0.2, label=r'CMASS N M')
ax.errorbar(tsCmassSM.RApArcmin+0.01, factor * tsCmassSM.tSZ, factor * np.sqrt(np.diag(tsCmassSM.covTsz)), fmt='-.', c='r', alpha=0.2, label=r'CMASS S M')
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{tSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
#ax.set_ylim((0., 2.))
#
path = tsCmassM.pathFig+"/tsz_cmass_mariana.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()


# plot CMASS Kendrick
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# CMASS K
ax.errorbar(tsCmassK.RApArcmin+0.02, factor * tsCmassK.tSZ, factor * np.sqrt(np.diag(tsCmassK.covTsz)), c='r', label=r'CMASS K')
ax.errorbar(tsCmassNK.RApArcmin, factor * tsCmassNK.tSZ, factor * np.sqrt(np.diag(tsCmassNK.covTsz)), fmt='--', c='r', alpha=0.2, label=r'CMASS N K')
ax.errorbar(tsCmassSK.RApArcmin+0.01, factor * tsCmassSK.tSZ, factor * np.sqrt(np.diag(tsCmassSK.covTsz)), fmt='-.', c='r', alpha=0.2, label=r'CMASS S K')
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{tSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
#ax.set_ylim((0., 10.))
#
path = tsCmassK.pathFig+"/tsz_cmass_kendrick.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()


# plot LOWZ Kendrick
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# CMASS K
ax.errorbar(tsLowzK.RApArcmin+0.02, factor * tsLowzK.tSZ, factor * np.sqrt(np.diag(tsLowzK.covTsz)), c='r', label=r'LOWZ K')
ax.errorbar(tsLowzNK.RApArcmin, factor * tsLowzNK.tSZ, factor * np.sqrt(np.diag(tsLowzNK.covTsz)), fmt='--', c='r', alpha=0.2, label=r'LOWZ N K')
ax.errorbar(tsLowzSK.RApArcmin+0.01, factor * tsLowzSK.tSZ, factor * np.sqrt(np.diag(tsLowzSK.covTsz)), fmt='-.', c='r', alpha=0.2, label=r'LOWZ S K')
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{tSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
#ax.set_ylim((0., 2.))
#
path = tsLowzK.pathFig+"/tsz_lowz_kendrick.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()


# plot BOSS Kendrick
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# CMASS K
ax.errorbar(tsBossK.RApArcmin, factor * tsBossK.tSZ, factor * np.sqrt(np.diag(tsBossK.covTsz)), fmt='--', c='r', label=r'BOSS K')
#
ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{tSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
#ax.set_xlim((0., 2.))
#
path = tsBossK.pathFig+"/tsz_kendrick.pdf"
fig.savefig(path, bbox_inches='tight')
fig.clf()
'''
