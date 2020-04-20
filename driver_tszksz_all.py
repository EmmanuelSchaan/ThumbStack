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
      # PACT day+night, 20200228
      "pactf150daynight20200228maskgal60": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits",
         "./output/cmb_map/pact20200228/" + "mask_full_foot_gal60_ps.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits",
         name="pactf150daynight20200228maskgal60"),
      "pactf90daynight20200228maskgal60": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_map.fits",
         "./output/cmb_map/pact20200228/" + "mask_full_foot_gal60_ps.fits",
         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_ivar.fits",
         name="pactf90daynight20200228maskgal60"),
      #
#      # PACT day+night, 20200228, old masks
#      "pactf150daynight20200228": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits",
#         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
#         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits",
#         name="pactf150daynight20200228"),
#      "pactf90daynight20200228": cmbMap("/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_map.fits",
#         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
#         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f090_daynight_ivar.fits",
#         name="pactf90daynight20200228"),
#      # Coadded PACT maps, daynight, 20200228
#      "coaddcmb20200228": cmbMap("./output/cmb_map/coadd/" + "coadd_hitmap_weighted.fits",
#         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
#         "./output/cmb_map/coadd/" + "coadd_hitmap_weighted_ivar.fits",
#         name="coaddcmb20200228"),
#      "coaddcmb20200228lmax9e3": cmbMap("./output/cmb_map/coadd/" + "coadd_hitmap_weighted_lmax9000.fits",
#         "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/" + "f150_mask_foot_planck_ps_car.fits",
#         "./output/cmb_map/coadd/" + "coadd_hitmap_weighted_lmax9000_ivar.fits",
#         name="coaddcmb20200228lmax9e3"),
      #
      # TileC v1.2, reconvolved to 1.4' beam, combining BOSS N and D56
#      "tilecpactcmbksz": cmbMap("./output/cmb_map/tilec_pact_cmbksz_v1.2.0/" + "tilec_reconv14_map.fits",
#      "./output/cmb_map/tilec_pact_cmbksz_v1.2.0/" + "mask_full_foot_gal_ps.fits",
#      name="tilecpactcmbksz"),
#      "tilecpactcmbksznoy": cmbMap("./output/cmb_map/tilec_pact_cmbksznoy_v1.2.0/" + "tilec_reconv14_map.fits",
#      "./output/cmb_map/tilec_pact_cmbksznoy_v1.2.0/" + "mask_full_foot_gal_ps.fits",
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

'''
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
'''

###################################################################################
# PACT 90 and 150: stacks and joint cov

import thumbstack
reload(thumbstack)
from thumbstack import *

#ts[freq] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc, doMBins=True)

save = False


for catalogKey in catalogs.keys():#[::-1]:
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

      ts[freq] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc, doMBins=True)
#      try:
#         ts[freq] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc, doMBins=True)
#      except:
#         ts[freq] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=True, nProc=nProc, doMBins=True)


   ###################################################################################
   # Joint covariance between 150 and 90

   # compute the joint cov
   if save:
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
   path = ts['150'].pathFig+"/summary_ksz_150_90_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()


   # tSZ plot at 150 and 90
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
   ax.set_ylabel(r'$T_\text{tSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(r'tSZ profile')
   #ax.set_ylim((0., 2.))
   #
   path = ts['150'].pathFig+"/summary_tsz_150_90_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()


###################################################################################
###################################################################################
# Null tests

# read the stacks on mock GRFs, to compare
pathMockGRF = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/output/cmb_map/mocks_grf_planck_act_coadd_2019_03_11/"
iMock0 = 0
nMocks = 800
est = 'ksz_varweight'
#
meanStackedGRF = np.genfromtxt(pathMockGRF+"mean_diskring_"+est+"_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".txt")
covStackedGRF = np.genfromtxt(pathMockGRF+"cov_diskring_"+est+"_mocks"+str(iMock0)+"-"+str(iMock0+nMocks)+".txt")
sStackedGRF = np.sqrt(np.diag(covStackedGRF)) / np.sqrt(nMocks)


      
for freq in ['90', '150']:
   cmbMapKey = "pactf"+freq+"daynight20200228maskgal60"
   cmbMap = cmbMaps[cmbMapKey].map()
   cmbMask = cmbMaps[cmbMapKey].mask()
   cmbHit = cmbMaps[cmbMapKey].hit()
   cmbName = cmbMaps[cmbMapKey].name
   print("Analyzing map "+cmbName)

   ts = {}
   for catalogKey in ['cmass_mariana', 'cmass_kendrick']:#catalogs.keys():
      catalog = catalogs[catalogKey]
      print("Analyzing catalog "+catalog.name)
      name = catalog.name + "_" + cmbName
      ts[catalogKey] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc, doMBins=True)



   # kSZ plot
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   # convert from sr to arcmin^2
   factor = (180.*60./np.pi)**2
   #
   ax.axhline(0., c='k', lw=1)
   #
   # Uncertainty band
   ax.fill_between(ts['cmass_mariana'].RApArcmin, - factor * ts['cmass_mariana'].sStackedProfile["diskring_ksz_varweight"], factor * ts['cmass_mariana'].sStackedProfile["diskring_ksz_varweight"], edgecolor='', facecolor='gray', alpha=0.5, label=r'statistical error')
   #
   # V-shuffle mean
   ax.errorbar(ts['cmass_mariana'].RApArcmin, factor * ts['cmass_mariana'].stackedProfile["diskring_ksz_varweight_vshufflemean"], factor * ts['cmass_mariana'].sStackedProfile["diskring_ksz_varweight_vshufflemean"], fmt='-', c='b', label='mean of v-shuffles')
   #
   # Mariana - Kendrick
   ax.plot(ts['cmass_mariana'].RApArcmin, factor * (ts['cmass_mariana'].sStackedProfile["diskring_ksz_varweight"] - ts['cmass_kendrick'].sStackedProfile["diskring_ksz_varweight"]), 'b-', label=r'$v_\text{Mariana} - v_\text{Kendrick}$', c='r')
   #
   # Average of many mocks
   ax.errorbar(ts['cmass_mariana'].RApArcmin + 0.05, factor * meanStackedGRF, yerr=factor*sStackedGRF, fmt='-', c='g', label=r'mean of '+str(nMocks)+' mocks')
   #
   ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_ylim((-2., 2.))
   #
   path = ts['cmass_mariana'].pathFig+"/nulltests_ksz_"+freq+"_cmass.pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()











