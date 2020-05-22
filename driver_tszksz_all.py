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


pathFig = '/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/figures/thumbstack/summary_plots/'


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


# Mini CMASS Mariana, for debugging
nObj = 10000#50000
cmassMarianaShort = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False, nObj=nObj)
cmassMarianaShort.name = "cmass_mariana_short"



catalogs = {
      #"cmass_mariana_short": cmassMarianaShort,
      #
      #"cmass_s_mariana": Catalog(u, massConversion, name="cmass_s_mariana", nameLong="CMASS S M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt", save=False),
      #"cmass_n_mariana": Catalog(u, massConversion, name="cmass_n_mariana", nameLong="CMASS N M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt", save=False),
      "cmass_mariana": Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", rV=0.5, save=False),
      #
      #"cmass_s_kendrick": Catalog(u, massConversion, name="cmass_s_kendrick", nameLong="CMASS S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_S_kendrick.txt", save=False),
      #"cmass_n_kendrick": Catalog(u, massConversion, name="cmass_n_kendrick", nameLong="CMASS N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_N_kendrick.txt", save=False),
      "cmass_kendrick": Catalog(u, massConversion, name="cmass_kendrick", nameLong="CMASS K", rV=0.7, save=False),
      #"lowz_s_kendrick": Catalog(u, massConversion, name="lowz_s_kendrick", nameLong="LOWZ S K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_S_kendrick.txt", save=False),
      #"lowz_n_kendrick": Catalog(u, massConversion, name="lowz_n_kendrick", nameLong="LOWZ N K", pathInCatalog="../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_N_kendrick.txt", save=False),
      "lowz_kendrick": Catalog(u, massConversion, name="lowz_kendrick", nameLong="LOWZ K", rV=0.7, save=False),
      #"boss_kendrick": Catalog(u, massConversion, name="boss_kendrick", nameLong="BOSS K", save=False),
      #"cmass_mk_diff": Catalog(u, massConversion, name="cmass_mk_diff", nameLong="CMASS M-K", rV=1., save=False),
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
      # PACT day+night, 20200228, Planck Galactic masks 60%
      "pactf150daynight20200228maskgal60": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits",
         "./output/cmb_map/pact20200228/" + "mask_full_foot_gal60_ps.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits",
         name="pactf150daynight20200228maskgal60"),
      "pactf90daynight20200228maskgal60": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_map.fits",
         "./output/cmb_map/pact20200228/" + "mask_full_foot_gal60_ps.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_ivar.fits",
         name="pactf90daynight20200228maskgal60"),
      "pactf150reconvto90minus90daynight20200228maskgal60": cmbMap("./output/cmb_map/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150reconvto90_minus_f090_daynight_map.fits",
         "./output/cmb_map/pact20200228/" + "mask_full_foot_gal60_ps.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits",
         name="pactf150reconvto90minus90daynight20200228maskgal60"),
      #
      #
      #"pactf150daynight20200228maskgal60reconvto90": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvto90.fits",
      #   "./output/cmb_map/pact20200228/" + "mask_full_foot_gal60_ps.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits",
      #   name="pactf150daynight20200228maskgal60reconvto90"),
      #
      #
      ## PACT day+night, 20200228, Planck Galactic masks 70%
      #"pactf150daynight20200228": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits",
      #   name="pactf150daynight20200228"),
      #"pactf90daynight20200228": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_map.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_ivar.fits",
      #   name="pactf90daynight20200228"),
      ## Coadded PACT maps, daynight, 20200228, Planck Galactic mask 70%
      #"coaddcmb20200228": cmbMap("./output/cmb_map/coadd/" + "coadd_hitmap_weighted.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
      #   "./output/cmb_map/coadd/" + "coadd_hitmap_weighted_ivar.fits",
      #   name="coaddcmb20200228"),
      #"coaddcmb20200228lmax9e3": cmbMap("./output/cmb_map/coadd/" + "coadd_hitmap_weighted_lmax9000.fits",
      #   "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
      #   "./output/cmb_map/coadd/" + "coadd_hitmap_weighted_lmax9000_ivar.fits",
      #   name="coaddcmb20200228lmax9e3"),
      #
      # TileC v1.2, reconvolved to 1.4' beam, combining BOSS N and D56
      "tilecpactynocmb": cmbMap("./output/cmb_map/tilec_pact_ynocmb_v1.2.0/" + "tilec_reconv14_map.fits",
      "./output/cmb_map/tilec_pact_ynocmb_v1.2.0/" + "mask_full_foot_gal_ps.fits",
      name="tilecpactynocmb"),
      "tilecpactyminusynocib": cmbMap("./output/cmb_map/tilec_pact_yminusynocib_v1.2.0/" + "tilec_reconv14_map.fits",
      "./output/cmb_map/tilec_pact_yminusynocib_v1.2.0/" + "mask_full_foot_gal_ps.fits",
      name="tilecpactyminusynocib"),
      #
      "tilecpacty": cmbMap("./output/cmb_map/tilec_pact_y_v1.2.0/" + "tilec_reconv14_map.fits",
      "./output/cmb_map/tilec_pact_y_v1.2.0/" + "mask_full_foot_gal_ps.fits",
      name="tilecpacty"),
      "tilecpactynocib": cmbMap("./output/cmb_map/tilec_pact_ynocib_v1.2.0/" + "tilec_reconv14_map.fits",
      "./output/cmb_map/tilec_pact_ynocib_v1.2.0/" + "mask_full_foot_gal_ps.fits",
      name="tilecpactynocib"),
      "tilecpactcmbksz": cmbMap("./output/cmb_map/tilec_pact_cmbksz_v1.2.0/" + "tilec_reconv14_map.fits",
      "./output/cmb_map/tilec_pact_cmbksz_v1.2.0/" + "mask_full_foot_gal_ps.fits",
      name="tilecpactcmbksz"),
      "tilecpactcmbksznoy": cmbMap("./output/cmb_map/tilec_pact_cmbksznoy_v1.2.0/" + "tilec_reconv14_map.fits",
      "./output/cmb_map/tilec_pact_cmbksznoy_v1.2.0/" + "mask_full_foot_gal_ps.fits",
      name="tilecpactcmbksznoy"),
      #
      # kSZ pipeline check
      "pactf150daynight20200228maskgal60_minus_tilecpactcmbksz": cmbMap("./output/cmb_map/pactf150daynight20200228maskgal60_minus_tilec_pact_cmbksz/" + "tilec_reconv14_map.fits", "./output/cmb_map/pactf150daynight20200228maskgal60_minus_tilec_pact_cmbksz/" + "mask_full_foot_gal_ps.fits", name="pactf150daynight20200228maskgal60_minus_tilecpactcmbksz"),
      #
      # kSZ dust contamination test
      "pactf150daynight20200228maskgal60_minus_tilecpactcmbksznocib": cmbMap("./output/cmb_map/pactf150daynight20200228maskgal60_minus_tilec_pact_cmbksznocib/" + "tilec_reconv14_map.fits", "./output/cmb_map/pactf150daynight20200228maskgal60_minus_tilec_pact_cmbksznocib/" + "mask_full_foot_gal_ps.fits", name="pactf150daynight20200228maskgal60_minus_tilecpactcmbksznocib"),
      #
      # tSZ pipeline (map) check
      "pactf150daynight20200228maskgal60_minus_tilecpactymuk": cmbMap("./output/cmb_map/pactf150daynight20200228maskgal60_minus_tilec_pact_ymuk/" + "tilec_reconv14_map.fits", "./output/cmb_map/pactf150daynight20200228maskgal60_minus_tilec_pact_ymuk/" + "mask_full_foot_gal_ps.fits", name="pactf150daynight20200228maskgal60_minus_tilecpactymuk"),
      #
      # tSZ internal dust test
      "pactf150reconvto90minus90noydaynight20200228maskgal60": cmbMap("./output/cmb_map/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150reconvto90_minus_f090_noy_daynight_map.fits",
         "./output/cmb_map/pact20200228/" + "mask_full_foot_gal60_ps.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits",
         name="pactf150reconvto90minus90noydaynight20200228maskgal60"),
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


catalogCombi = {
      "pactf150daynight20200228maskgal60": ['cmass_kendrick', 'lowz_kendrick', 'cmass_mariana'],
      "pactf90daynight20200228maskgal60": ['cmass_kendrick', 'lowz_kendrick', 'cmass_mariana'],
      "pactf150reconvto90minus90daynight20200228maskgal60": ['cmass_kendrick'],
      #
      "tilecpactynocmb": ['cmass_kendrick'],
      "tilecpactyminusynocib": ['cmass_kendrick'],
      #
      "tilecpacty": ['cmass_kendrick'],
      "tilecpactynocib": ['cmass_kendrick', 'lowz_kendrick'],
      "tilecpactcmbksz": ['cmass_kendrick'],
      "tilecpactcmbksznoy": ['cmass_kendrick'],
      #
      "pactf150daynight20200228maskgal60_minus_tilecpactcmbksz": ['cmass_kendrick'],
      "pactf150daynight20200228maskgal60_minus_tilecpactymuk": ['cmass_kendrick'],
      "pactf150reconvto90minus90noydaynight20200228maskgal60": ['cmass_kendrick'],
      "pactf150daynight20200228maskgal60_minus_tilecpactcmbksznocib": ['cmass_kendrick'],
      }

#catalogCombi = {
#      "pactf150daynight20200228maskgal60": ['cmass_mariana', 'cmass_kendrick', 'lowz_kendrick', 'cmass_mk_diff'],
#      "pactf90daynight20200228maskgal60": [],#['cmass_mariana', 'cmass_kendrick', 'lowz_kendrick'],
#      "pactf150reconvto90minus90daynight20200228maskgal60": ['cmass_mariana', 'cmass_kendrick'],
#      #
#      "tilecpactynocmb": [],#['cmass_kendrick'],
#      "tilecpactyminusynocib": [],#['cmass_kendrick'],
#      #
#      "tilecpacty": [],#['cmass_kendrick'],
#      "tilecpactynocib": [],#['cmass_kendrick'],
#      "tilecpactcmbksz": [],#['cmass_kendrick'],
#      "tilecpactcmbksznoy": [],#['cmass_kendrick'],
#      #
#      "pactf150daynight20200228maskgal60_minus_tilecpactcmbksz": [],#['cmass_kendrick'],
#      "pactf150daynight20200228maskgal60_minus_tilecpactymuk": [],#['cmass_kendrick'],
#      "pactf150reconvto90minus90noydaynight20200228maskgal60": [],
#      "pactf150daynight20200228maskgal60_minus_tilecpactcmbksznocib": ['cmass_mariana'],#['cmass_kendrick'],
#      }



###################################################################################
###################################################################################
# Compute all the stacked profiles
'''
import thumbstack
reload(thumbstack)
from thumbstack import *


save = True


for cmbMapKey in cmbMaps.keys():
#for cmbMapKey in ['pactf150daynight20200228maskgal60', 'pactf90daynight20200228maskgal60', 'tilecpactynocib']:
#for cmbMapKey in ['tilecpactynocib', 'pactf90daynight20200228maskgal60','pactf150daynight20200228maskgal60']:
   cmbMap = cmbMaps[cmbMapKey].map()
   cmbMask = cmbMaps[cmbMapKey].mask()
   cmbHit = cmbMaps[cmbMapKey].hit()
   cmbName = cmbMaps[cmbMapKey].name
   print("Analyzing map "+cmbName)

   for catalogKey in catalogCombi[cmbMapKey]:
      catalog = catalogs[catalogKey]
      print("Analyzing catalog "+catalog.name)
      name = catalog.name + "_" + cmbName

#      ts = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc, doMBins=False, doBootstrap=False, doVShuffle=False)

      try:
         ts = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=False, nProc=nProc, doMBins=True, doBootstrap=True, doVShuffle=False)
      except:
         ts = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=True, nProc=nProc, doMBins=True, doBootstrap=True, doVShuffle=False)

      #ts.plotAllStackedProfiles()
'''
###################################################################################
###################################################################################
#  PACT 90 and 150: stacks and joint cov

import thumbstack
reload(thumbstack)
from thumbstack import *

#ts[freq] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc, doMBins=True)

save = False


for catalogKey in ['cmass_kendrick', 'lowz_kendrick']:
   catalog = catalogs[catalogKey]
   print("Analyzing catalog "+catalog.name)
   
   ts = {}
   for freq in ['90', '150']:
      cmbMapKey = "pactf"+freq+"daynight20200228maskgal60"
      cmbMap = cmbMaps[cmbMapKey].map()
      cmbMask = cmbMaps[cmbMapKey].mask()
      cmbHit = cmbMaps[cmbMapKey].hit()
      cmbName = cmbMaps[cmbMapKey].name
      print("Analyzing map "+cmbName)
      name = catalog.name + "_" + cmbName

      ts[freq] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc, doMBins=False, doBootstrap=False, doVShuffle=False)
#      try:
#         ts[freq] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc, doMBins=False, doBootstrap=False, doVShuffle=False)
#      except:
#         ts[freq] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=True, nProc=nProc, doMBins=False, doBootstrap=False, doVShuffle=False)


   ###################################################################################
   # Joint covariance between 150 and 90

   # compute the joint cov
   if True:
      ts['150'].saveAllCovBootstrapTwoStackedProfiles(ts['90'])
   ts['150'].plotAllCovTwoStackedProfiles(ts['90'])

   ###################################################################################
   # Summary kSZ and tSZ at 150 and 90


   # kSZ plot at 150 and 90
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   # convert from sr to arcmin^2
   factor = (180.*60./np.pi)**2
   #
   ax.axhline(0., c='k', lw=1)
   #
   ax.errorbar(ts['150'].RApArcmin, factor * ts['150'].stackedProfile["diskring_ksz_varweight"], factor * ts['150'].sStackedProfile["diskring_ksz_varweight"], fmt='-', c='royalblue', label='150GHz')
   ax.errorbar(ts['150'].RApArcmin + 0.05, factor * ts['90'].stackedProfile["diskring_ksz_varweight"], factor * ts['90'].sStackedProfile["diskring_ksz_varweight"], fmt='-', c='darkviolet', label='90GHz')
   #
   ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(r'kSZ profile')
   ax.set_ylim((0., 10.))
   #
   #path = ts['150'].pathFig+"summary_ksz_150_90_"+catalogKey+".pdf"
   path = pathFig+"summary_ksz_150_90_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()


   # tSZ + dust plot at 150 and 90
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   # convert from sr to arcmin^2
   factor = (180.*60./np.pi)**2
   #
   ax.axhline(0., c='k', lw=1)
   #
   ax.errorbar(ts['150'].RApArcmin, factor * ts['150'].stackedProfile["diskring_tsz_varweight"], factor * ts['150'].sStackedProfile["diskring_tsz_varweight"], fmt='-', c='royalblue', label='150GHz')
   ax.errorbar(ts['90'].RApArcmin + 0.05, factor * ts['90'].stackedProfile["diskring_tsz_varweight"], factor * ts['90'].sStackedProfile["diskring_tsz_varweight"], fmt='-', c='darkviolet', label='90GHz')
   #
   ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_{\text{tSZ} + \text{dust}}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(r'tSZ + dust profile')
   #ax.set_ylim((0., 2.))
   #
   #path = ts['150'].pathFig+"summary_tsz_150_90_"+catalogKey+".pdf"
   path = pathFig+"summary_tsz_150_90_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()




###################################################################################
###################################################################################
# Read all the stacked profiles
'''

# convert from sr to arcmin^2
factor = (180.*60./np.pi)**2

# conversion from y to muK at 150 GHz
Tcmb = 2.726   # K
h = 6.63e-34   # SI
kB = 1.38e-23  # SI
def f(nu):
   """frequency dependence for tSZ temperature
   """
   x = h*nu/(kB*Tcmb)
   return x*(np.exp(x)+1.)/(np.exp(x)-1.) -4.
yTomuK150 = f(150.e9) * Tcmb * 1.e6  # [muK * sr]




pathThumb = './output/thumbstack/'

rKsz = {}
ksz = {}
sKsz = {}
#
rTsz = {}
tsz = {}
sTsz = {}

for cmbMapKey in cmbMaps.keys():
   cmbName = cmbMaps[cmbMapKey].name

   for catalogKey in catalogCombi[cmbMapKey]:
      catalog = catalogs[catalogKey]
      name = catalog.name + "_" + cmbName

      # read the stacked kSZ profile
      try:
         data = np.genfromtxt(pathThumb + name + "/diskring_ksz_varweight_measured.txt")
      except:
         data = np.genfromtxt(pathThumb + name + "/diskring_ksz_uniformweight_measured.txt")
      rKsz[name] = data[:,0]
      ksz[name] = data[:,1] * factor
      sKsz[name] = data[:,2] * factor

      # read the stacked tSZ profile
      try:
         data = np.genfromtxt(pathThumb + name + "/diskring_tsz_varweight_measured.txt")
      except:
         data = np.genfromtxt(pathThumb + name + "/diskring_tsz_uniformweight_measured.txt")
      rTsz[name] = data[:,0]
      tsz[name] = data[:,1] * factor
      sTsz[name] = data[:,2] * factor




# read the stacks on mock GRFs, to compare
pathMockGRF = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/output/cmb_map/mocks_grf_planck_act_coadd_2019_03_11/"
iMock0 = 0
nMocks = 800
#
est = 'ksz_varweight'
meanStackedKszGRF = np.genfromtxt(pathMockGRF+"mean_diskring_"+est+"_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".txt") * factor
covStackedKszGRF = np.genfromtxt(pathMockGRF+"cov_diskring_"+est+"_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".txt") * factor**2
sStackedKszGRF = np.sqrt(np.diag(covStackedKszGRF)) / np.sqrt(nMocks) 
#
est = 'tsz_varweight'
meanStackedTszGRF = np.genfromtxt(pathMockGRF+"mean_diskring_"+est+"_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".txt") * factor
covStackedTszGRF = np.genfromtxt(pathMockGRF+"cov_diskring_"+est+"_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".txt") * factor**2
sStackedTszGRF = np.sqrt(np.diag(covStackedTszGRF)) / np.sqrt(nMocks)




# kSZ: v-shuffle mean
data = np.genfromtxt(pathThumb + "cmass_kendrick_pactf150daynight20200228maskgal60/" + "diskring_ksz_varweight_vshufflemean.txt")
rKsz150VShuffleMean = data[:,0]
ksz150VShuffleMean = data[:,1] * factor
sKsz150VShuffleMean = data[:,2] * factor



###################################################################################
###################################################################################
# kSZ null tests

rAp = rKsz['cmass_kendrick_pactf150daynight20200228maskgal60']
#
# fiducial uncertainty
ksz150 = ksz['cmass_kendrick_pactf150daynight20200228maskgal60']
sKsz150 = sKsz['cmass_kendrick_pactf150daynight20200228maskgal60']
# Mariana's velocities
ksz150Mariana = ksz['cmass_mariana_pactf150daynight20200228maskgal60']
sKsz150Mariana = sKsz['cmass_mariana_pactf150daynight20200228maskgal60']
#
# 150 - tilec cmb, as a consistency check
ksz150MinusTilecCmb = ksz['cmass_kendrick_pactf150daynight20200228maskgal60_minus_tilecpactcmbksz']
sKsz150MinusTilecCmb = sKsz['cmass_kendrick_pactf150daynight20200228maskgal60_minus_tilecpactcmbksz']
#
# 150 reconv to 90 minus 90, to check consistency
ksz150Reconv90Minus90 = ksz['cmass_kendrick_pactf150reconvto90minus90daynight20200228maskgal60']
sKsz150Reconv90Minus90 = sKsz['cmass_kendrick_pactf150reconvto90minus90daynight20200228maskgal60']
#
# tilec y no cmb, to check for tSZ contamination
kszYNoCmb = ksz['cmass_kendrick_tilecpactynocmb'] * yTomuK150
sKszYNoCmb = sKsz['cmass_kendrick_tilecpactynocmb'] * yTomuK150
#
# 150 - tilec cmb no cib, to check for dust contamination
ksz150MinusCmbNoCib = ksz['cmass_kendrick_pactf150daynight20200228maskgal60_minus_tilecpactcmbksznocib']
sKsz150MinusCmbNoCib = sKsz['cmass_kendrick_pactf150daynight20200228maskgal60_minus_tilecpactcmbksznocib']




# kSZ pipeline null tests
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
#
ax.axhline(0., c='k', lw=1)
#
# Uncertainty band
ax.fill_between(rAp, - sKsz150, sKsz150, edgecolor='', facecolor='gray', alpha=0.5, label=r'statistical error')
#
# V-shuffle mean
ax.errorbar(rAp, ksz150VShuffleMean, yerr=sKsz150VShuffleMean, fmt='--', label='mean of 100 v-shuffles')
#
# Average of many mocks
ax.errorbar(rAp + 0.025, meanStackedKszGRF, yerr=sStackedKszGRF, fmt='--', label=r'mean of '+str(nMocks)+' mocks')
#
# Mariana - Kendrick
#ax.errorbar(rAp, ksz150MKDiff, yerr=sKsz150MKDiff, fmt='-', c='r', label=r'$v_\text{Mariana} - v_\text{Kendrick}$')
#ax.errorbar(rAp, ksz150, yerr=sKsz150, fmt='-', label='K')
#ax.errorbar(rAp, ksz150Mariana, yerr=sKsz150Mariana, fmt='-', label='M')
ax.errorbar(rAp + 0.05, (ksz150-ksz150Mariana), yerr=sKsz150, fmt='-', label=r'$v_\text{Kendrick} - v_\text{Mariana}$')
#
# 150 - tilec cmb
ax.errorbar(rAp + 0.075, ksz150MinusTilecCmb, yerr=sKsz150MinusTilecCmb, fmt='-', label='150 - TileC CMB/kSZ')
#
# 150 reconv to 90 minus 90
ax.errorbar(rAp + 0.1, ksz150Reconv90Minus90, yerr=sKsz150Reconv90Minus90, fmt='-', label='150\' - 90')
#
ax.set_ylim((-10., 15.))
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
ax.set_title(r'kSZ pipeline null tests')
ax.set_ylim((-6., 6.))
#
path = pathFig + "pipenulltests_ksz_150_cmass.pdf"
fig.savefig(path, bbox_inches='tight')
#plt.show()
fig.clf()




# kSZ foreground null tests
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# convert from sr to arcmin^2
factor = (180.*60./np.pi)**2
#
ax.axhline(0., c='k', lw=1)
#
# Uncertainty band
ax.fill_between(rAp, - sKsz150, sKsz150, edgecolor='', facecolor='gray', alpha=0.5, label=r'statistical error')
#
# kSZ on TileC y no CMB map
ax.errorbar(rAp, kszYNoCmb, yerr=sKszYNoCmb, label=r'TileC y no CMB')
#
# cmbksz no cib, to check for dust
ax.errorbar(rAp + 0.025, ksz150MinusCmbNoCib, yerr=sKsz150MinusCmbNoCib, fmt='-', label=r'150 - TileC CMB/kSZ no CIB')
#
ax.set_ylim((-10., 15.))
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
ax.set_title(r'kSZ foreground null tests')
#ax.set_ylim((-2., 2.))
#
path = pathFig + "fgnulltests_ksz_150_cmass.pdf"
fig.savefig(path, bbox_inches='tight')
#plt.show()
fig.clf()



###################################################################################
###################################################################################
# tSZ null tests

rAp = rTsz['cmass_kendrick_pactf150daynight20200228maskgal60']
#
# fiducial uncertainty
sTsz150 = sTsz['cmass_kendrick_pactf150daynight20200228maskgal60'] 
#
# 150 - y, to check for map consistency
tsz150MinusY = tsz['cmass_kendrick_pactf150daynight20200228maskgal60_minus_tilecpactymuk']
sTsz150MinusY = sTsz['cmass_kendrick_pactf150daynight20200228maskgal60_minus_tilecpactymuk']
#
# y - y no CIB
tszYMinusYNoCib = tsz['cmass_kendrick_tilecpactyminusynocib'] * yTomuK150
sTszYMinusYNoCib = sTsz['cmass_kendrick_tilecpactyminusynocib'] * yTomuK150
#
# 150' - 90, after rescaling 90 to null tSZ
# in order to check for the dust contamination
tsz150Reconv90Minus90NoY = tsz['cmass_kendrick_pactf150reconvto90minus90noydaynight20200228maskgal60']
sTsz150Reconv90Minus90NoY = sTsz['cmass_kendrick_pactf150reconvto90minus90noydaynight20200228maskgal60']



# tSZ pipeline test
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.axhline(0., c='k', lw=1)
#
# Uncertainty band
ax.fill_between(rAp, - sTsz150, sTsz150, edgecolor='', facecolor='gray', alpha=0.5, label=r'statistical error')
#
# mean of GRF mocks
ax.errorbar(rAp, meanStackedTszGRF, yerr=sStackedTszGRF, fmt='--', label=r'mean of '+str(nMocks)+' mocks')
#
# 150 - tilec y
ax.errorbar(rAp, tsz150MinusY, yerr=sTsz150MinusY, fmt='--', label=r'150 - TileC y ')
#
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{tSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
ax.set_title(r'tSZ pipeline null tests')
ax.set_ylim((-6., 6.))
#
path = pathFig + "pipenulltests_tsz_150_cmass.pdf"
fig.savefig(path, bbox_inches='tight')
#plt.show()
fig.clf()



# dust contamination to tSZ
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.axhline(0., c='k', lw=1)
#
# Uncertainty band
ax.fill_between(rAp, - sTsz150, sTsz150, edgecolor='', facecolor='gray', alpha=0.5, label=r'statistical error')
#
# 150 - tilec y
ax.errorbar(rAp, tsz150MinusY, yerr=sTsz150MinusY, fmt='--', label=r'150 - TileC y')
#
# y - y no CIB
ax.errorbar(rAp, tszYMinusYNoCib, yerr=sTszYMinusYNoCib, fmt='--', label=r'TileC y - y no CIB')
#
# 150' - 90 rescaled to null y
ax.errorbar(rAp, tsz150Reconv90Minus90NoY, yerr=sTsz150Reconv90Minus90NoY, fmt='--', label=r"150\' - 90 no y")
#
ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_\text{dust}$ [$\mu K\cdot\text{arcmin}^2$]')
ax.set_title(r'Dust emission')
ax.set_ylim((-6., 6.))
#
path = pathFig + "fgnulltests_tsz_150_cmass.pdf"
fig.savefig(path, bbox_inches='tight')
#plt.show()
fig.clf()



###################################################################################
###################################################################################
# summary tSZ plot


rAp = rTsz['cmass_kendrick_pactf150daynight20200228maskgal60']

# PACT 150
tsz150 = tsz['cmass_kendrick_pactf150daynight20200228maskgal60']
sTsz150 = sTsz['cmass_kendrick_pactf150daynight20200228maskgal60']
#
# PACT 90
tsz90 = tsz['cmass_kendrick_pactf90daynight20200228maskgal60']
sTsz90 = sTsz['cmass_kendrick_pactf90daynight20200228maskgal60']
#
# TileC y no Cib
tszYNoCib = tsz['cmass_kendrick_tilecpactynocib'] * yTomuK150
sTszYNoCib = sTsz['cmass_kendrick_tilecpactynocib'] * yTomuK150


# tSZ-only from the TileC y no CIB map
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# convert from sr to arcmin^2
factor = (180.*60./np.pi)**2
#
ax.axhline(0., c='k', lw=1)
#
# PACT 150
#ax.errorbar(rAp, tsz150, yerr=sTsz150, fmt='-', c='royalblue', label='150GHz')
# PACT 90
#ax.errorbar(rAp, tsz90, yerr=sTsz90, fmt='-', c='darkviolet', label='90GHz')
# Tilec y no CIB
ax.errorbar(rAp, tszYNoCib, yerr=sTszYNoCib, fmt='-', c='r', label='TileC y no CIB')
#
ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_{\text{tSZ}}$ [$\mu K\cdot\text{arcmin}^2$]')
ax.set_title(r'tSZ profile')
#ax.set_ylim((0., 2.))
#
path = pathFig+"summary_tsz_"+catalogKey+".pdf"
fig.savefig(path, bbox_inches='tight')
#plt.show()
fig.clf()


# tSZ + dust plot at 150 and 90
fig=plt.figure(0)
ax=fig.add_subplot(111)
#
# convert from sr to arcmin^2
factor = (180.*60./np.pi)**2
#
ax.axhline(0., c='k', lw=1)
#
# PACT 150
ax.errorbar(rAp, tsz150, yerr=sTsz150, fmt='-', c='royalblue', label='150GHz')
# PACT 90
ax.errorbar(rAp, tsz90, yerr=sTsz90, fmt='-', c='darkviolet', label='90GHz')
# Tilec y no CIB
ax.errorbar(rAp, tszYNoCib, yerr=sTszYNoCib, fmt='-', label='TileC y no CIB')
#
ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
ax.set_xlabel(r'$R$ [arcmin]')
ax.set_ylabel(r'$T_{\text{tSZ} + \text{dust}}$ [$\mu K\cdot\text{arcmin}^2$]')
ax.set_title(r'tSZ + dust profile')
#ax.set_ylim((0., 2.))
#
path = pathFig+"comparison_tsz_150_90_"+catalogKey+".pdf"
fig.savefig(path, bbox_inches='tight')
#plt.show()
fig.clf()


'''



###################################################################################
###################################################################################





