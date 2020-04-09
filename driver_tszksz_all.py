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


pathFig = '/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/summary_plots'


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
      #"cmass_s_mariana": Catalog(u, massConversion, name="cmass_s_mariana", nameLong="CMASS S M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt", save=False),
      #"cmass_n_mariana": Catalog(u, massConversion, name="cmass_n_mariana", nameLong="CMASS N M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt", save=False),
      "cmass_mariana": Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False),
      #
      #"cmass_s_kendrick": Catalog(u, massConversion, name="cmass_s_kendrick", nameLong="CMASS S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_S_kendrick.txt", save=False),
      #"cmass_n_kendrick": Catalog(u, massConversion, name="cmass_n_kendrick", nameLong="CMASS N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_N_kendrick.txt", save=False),
      "cmass_kendrick": Catalog(u, massConversion, name="cmass_kendrick", nameLong="CMASS K", save=False),
      #"lowz_s_kendrick": Catalog(u, massConversion, name="lowz_s_kendrick", nameLong="LOWZ S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_S_kendrick.txt", save=False),
      #"lowz_n_kendrick": Catalog(u, massConversion, name="lowz_n_kendrick", nameLong="LOWZ N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_N_kendrick.txt", save=False),
      "lowz_kendrick": Catalog(u, massConversion, name="lowz_kendrick", nameLong="LOWZ K", save=False),
      #"boss_kendrick": Catalog(u, massConversion, name="boss_kendrick", nameLong="BOSS K", save=False)
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
      result = enmap.read_map(self.pathMask)
      # if the map contains polarization, keep only temperature
      if len(result.shape)>2:
         result = result[0]
      return result
   
   def hit(self):
      if self.pathHit is None:
         return None
      else:
         result = enmap.read_map(self.pathHit)
         # if the map contains polarization, keep only temperature
         if len(result.shape)>2:
            result = result[0]
         return result


print("Read CMB maps")
tStart = time()

cmbMaps = {
      # Coadded PACT maps, daynight, 20200228
#      "coaddcmb20200228": cmbMap("./output/cmb_map/coadd/" + "coadd_hitmap_weighted.fits",
#         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
#         "./output/cmb_map/coadd/" + "coadd_hitmap_weighted_ivar.fits",
#         name="coaddcmb20200228"),
#      "coaddcmb20200228lmax9e3": cmbMap("./output/cmb_map/coadd/" + "coadd_hitmap_weighted_lmax9000.fits",
#         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
#         "./output/cmb_map/coadd/" + "coadd_hitmap_weighted_lmax9000_ivar.fits",
#         name="coaddcmb20200228lmax9e3"),
#      "pactf150daynight20200228": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits",
#         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
#         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits",
#         name="pactf150daynight20200228"),
      "pactf90daynight20200228": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_map.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_ivar.fits",
         name="pactf90daynight20200228"),
      #
      # TileC v1.2, reconvolved to 1.4' beam, combining BOSS N and D56
#      "tilecpactcmbksz": cmbMap("./output/cmb_map/tilec_pact_cmbksz_v1.2.0/" + "tilec_reconv14_map.fits",
#      "./output/cmb_map/tilec_pact_cmbksz_v1.2.0/" + "mask_foot_planck_ps_car.fits",
#      name="tilecpactcmbksz"),
#      "tilecpactcmbksznoy": cmbMap("./output/cmb_map/tilec_pact_cmbksznoy_v1.2.0/" + "tilec_reconv14_map.fits",
#      "./output/cmb_map/tilec_pact_cmbksznoy_v1.2.0/" + "mask_foot_planck_ps_car.fits",
#      name="tilecpactcmbksznoy"),
      #
      ## TileC v1.1, reconvolved to 1.4' beam, combining BOSS N and D56
      #"tilecpactcmbksz_14": cmbMap("./output/cmb_map/tilec_pact_cmbksz/" + "tilec_reconv14_map.fits",
      #"./output/cmb_map/tilec_pact_cmbksz/" + "mask_foot_planck_ps_car.fits",
      #name="tilecpactcmbksz_14"),
      #
      #
      ## PACT maps, night only, 20200228
      #"pactf220night20200228": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f220_night_map.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f220_night_ivar.fits",
      #   name="pactf220night20200228"),
      #"pactf150night20200228": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_night_map.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_night_ivar.fits",
      #   name="pactf150night20200228"),
      #"pactf090night20200228": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_night_map.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_night_ivar.fits",
      #   name="pactf090night20200228"),
      #
      #
      ## Old PACT maps
      #"pactf150night20190311": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "act_planck_f150_prelim_map_mono.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "act_planck_f150_prelim_div_mono.fits",
      #   name="pactf150night20190311"),
      #"pactf090night20190311": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "act_planck_f090_prelim_map_mono.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f090_mask_foot_planck_ps_car.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "act_planck_f090_prelim_div_mono.fits",
      #   name="pactf090night20190311"),
      #
      #
      ## Planck maps
      #"plancksmica18": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck_smicacmb.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck_nonact_gal_ps_mask.fits",
      #   name="plancksmica18"),
      #"plancksmicanosz18": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck_smicacmbnosz.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck_nonact_gal_ps_mask.fits",
      #   name="plancksmicanosz18"),
      #"planck54518": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck545.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_maps_2018_PR3/" + "car_planck_gal_ps_mask.fits",
      #   name="planck54518"),
      #
      #
      # TileC maps
      #"tilecpactcmbkszd56": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_deep56/" + "tilec_single_tile_deep56_cmb_map_v1.0.0_rc_joint.fits",
      #   "./output/cmb_map/tilec_pact_cmbksz_d56/" + "mask_foot_planck_ps_car.fits",
      #   name="tilecpactcmbkszd56"),
      #"tilecpactyd56": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_deep56/" + "tilec_single_tile_deep56_comptony_map_v1.0.0_rc_joint.fits",
      #   "./output/cmb_map/tilec_pact_cmbksz_d56/" + "mask_foot_planck_ps_car.fits",
      #   name="tilecpactyd56"),
      #"tilecpactynocibd56": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_deep56/" + "tilec_single_tile_deep56_comptony_deprojects_cib_map_v1.0.0_rc_joint.fits",
      #   "./output/cmb_map/tilec_pact_cmbksz_d56/" + "mask_foot_planck_ps_car.fits",
      #   name="tilecpactynocibd56"),
      #"tilecpactcmbksznoyd56": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_deep56/" + "tilec_single_tile_deep56_cmb_deprojects_comptony_map_v1.0.0_rc_joint.fits",
      #   "./output/cmb_map/tilec_pact_cmbksz_d56/" + "mask_foot_planck_ps_car.fits",
      #   name="tilecpactcmbksznoyd56"),
      #"tilecpactcmbksznocibd56": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_deep56/" + "tilec_single_tile_deep56_cmb_deprojects_cib_map_v1.0.0_rc_joint.fits",
      #   "./output/cmb_map/tilec_pact_cmbksz_d56/" + "mask_foot_planck_ps_car.fits",
      #   name="tilecpactcmbksznocibd56"),
      # 
      #"tilecpactcmbkszbossn": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_boss/" + "tilec_single_tile_boss_cmb_map_v1.0.0_rc_joint.fits",
      #   "./output/cmb_map/tilec_pact_cmbksz_boss/" + "mask_foot_planck_ps_car.fits",
      #   name="tilecpactcmbkszbossn"),
      #"tilecpactybossn": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_boss/" + "tilec_single_tile_boss_comptony_map_v1.0.0_rc_joint.fits",
      #   "./output/cmb_map/tilec_pact_cmbksz_boss/" + "mask_foot_planck_ps_car.fits",
      #   name="tilecpactybossn"),
      #"tilecpactynocibbossn": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_boss/" + "tilec_single_tile_boss_comptony_deprojects_cib_map_v1.0.0_rc_joint.fits",
      #   "./output/cmb_map/tilec_pact_cmbksz_boss/" + "mask_foot_planck_ps_car.fits",
      #   name="tilecpactynocibbossn"),
      #"tilecpactcmbksznoybossn": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_boss/" + "tilec_single_tile_boss_cmb_deprojects_comptony_map_v1.0.0_rc_joint.fits",
      #   "./output/cmb_map/tilec_pact_cmbksz_boss/" + "mask_foot_planck_ps_car.fits",
      #   name="tilecpactcmbksznoybossn"),
      #"tilecpactcmbksznocibbossn": cmbMap("/global/cscratch1/sd/msyriac/data/depot/tilec/v1.0.0_rc_20190919/map_v1.0.0_rc_joint_boss/" + "tilec_single_tile_boss_cmb_deprojects_cib_map_v1.0.0_rc_joint.fits",
      #   "./output/cmb_map/tilec_pact_cmbksz_boss/" + "mask_foot_planck_ps_car.fits",
      #   name="tilecpactcmbksznocibbossn")
      }

tStop = time()
print("took "+str(round((tStop-tStart)/60., 2))+" min")


###################################################################################
###################################################################################
# Do the stacking


import thumbstack
reload(thumbstack)
from thumbstack import *


save = True


for cmbMapKey in cmbMaps.keys():
   cmbMap = cmbMaps[cmbMapKey].map()
   cmbMask = cmbMaps[cmbMapKey].mask()
   cmbHit = cmbMaps[cmbMapKey].hit()
   cmbName = cmbMaps[cmbMapKey].name
   print("Analyzing map "+cmbName)

   for catalogKey in catalogs.keys():
      catalog = catalogs[catalogKey]
      print("Analyzing catalog "+catalog.name)
      name = catalog.name + "_" + cmbName
      try:
         ts = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc, doMBins=True)
      except:
         ts = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=True, nProc=nProc, doMBins=True)



###################################################################################
# Create summary plots

"""
# loop over maps
for cmbMapKey in cmbMaps.keys():
   cmbName = cmbMaps[cmbMapKey].name
   print("Loading map "+cmbMaps[cmbMapKey].name)
   cmbMap = cmbMaps[cmbMapKey].map()
   cmbMask = cmbMaps[cmbMapKey].mask()
   cmbHit = cmbMaps[cmbMapKey].hit()

   # load stacking on all catalogs
   print("Loading all catalogs stacked on this map:")
   ts = {}
   for catalogKey in catalogs.keys():
      catalog = catalogs[catalogKey]
      name = catalogs[catalogKey].name + "_" + cmbName
      try:
         ts[catalogKey] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=False, nProc=nProc)
      except:
         print("Could not load stacking for catalog "+catalog.name)

   print("Plotting all stacked profiles on map "+cmbMaps[cmbMapKey].name)
   # loop over estimators 
   # the estimators available depend only on the presence of a hit map
   ts0 = ts[ts.keys()[0]]
   for est in ts0.Est:
      print("Estimator "+est)
      
      #plot CMASS Mariana: N, S, N+S
      tsArr = [ts['cmass_n_mariana'], ts['cmass_s_mariana'], ts['cmass_mariana']]
      ts0.plotStackedProfile('diskring', [est], name='diskring_'+est+'_'+cmbName+'_cmassm', pathDir=pathFig, theory=False, tsArr=tsArr, plot=False)

      #plot CMASS Kendrick: N, S, N+S
      tsArr = [ts['cmass_n_kendrick'], ts['cmass_s_kendrick'], ts['cmass_kendrick']]
      ts0.plotStackedProfile('diskring', [est], name='diskring_'+est+'_'+cmbName+'_cmassk', pathDir=pathFig, theory=False, tsArr=tsArr, plot=False)

      #plot LOWZ Kendrick: N, S, N+S
      tsArr = [ts['lowz_n_kendrick'], ts['lowz_s_kendrick'], ts['lowz_kendrick']]
      ts0.plotStackedProfile('diskring', [est], name='diskring_'+est+'_'+cmbName+'_lowzk', pathDir=pathFig, theory=False, tsArr=tsArr, plot=False)

      #plot BOSS Kendrick: CMASS, LOWZ, CMASS+LOWZ
      tsArr = [ts['cmass_kendrick'], ts['lowz_kendrick'], ts['boss_kendrick']]
      ts0.plotStackedProfile('diskring', [est], name='diskring_'+est+'_'+cmbName+'_bossk', pathDir=pathFig, theory=False, tsArr=tsArr, plot=False)
"""




