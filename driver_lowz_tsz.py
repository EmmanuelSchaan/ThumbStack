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


pathFig = "/global/cscratch1/sd/eschaan/project_ucsc/code/ThumbStack/figures/thumstack/summary_plots/"

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


tStop = time()
print("took "+str(round((tStop-tStart)/60., 2))+" min")


catalogs = {
      "lowz20200908full": Catalog(u, massConversion, name="lowz20200908full", nameLong="LOWZ Full", save=False),
      "lowz20200908mbin0": Catalog(u, massConversion, name="lowz20200908mbin0", nameLong="LOWZ Bin 0", save=False),
      "lowz20200908mbin1": Catalog(u, massConversion, name="lowz20200908mbin1", nameLong="LOWZ Bin 1", save=False),
      "lowz20200908mbin2": Catalog(u, massConversion, name="lowz20200908mbin2", nameLong="LOWZ Bin 2", save=False)
      }


# mini LOWZ, for debugging
nObj = 5000
lowzShort = Catalog(u, massConversion, name="lowzfull20200908", nameLong="LOWZ Full", save=False, nObj=nObj)
lowzShort.name = "lowzshort20200908"
lowzShort.nameLong = "LOWZ Short"



###################################################################################
# Read CMB maps

class cmbMap(object):

   def __init__(self, pathMap, pathMask, pathHit=None,  name="test", nu=150.e9, unitLatex=r'$\mu$K'):
      self.name = name
      self.pathMap = pathMap
      self.pathMask = pathMask
      self.pathHit = pathHit

      self.nu = nu
      self.unitLatex = unitLatex

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
      "pactf150daynight20200228maskgal60r2": cmbMap("/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits",
         "./output/cmb_map/pact20200228_r2/" + "mask_full_foot_gal60_ps.fits",
         "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits",
         nu=150.e9, unitLatex=r'$\mu$K',
         name="pactf150daynight20200228maskgal60r2"),
      "pactf90daynight20200228maskgal60r2": cmbMap("/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f090_daynight_map.fits",
         "./output/cmb_map/pact20200228_r2/" + "mask_full_foot_gal60_ps.fits",
         "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f090_daynight_ivar.fits",
         nu=90.e9, unitLatex=r'$\mu$K',
         name="pactf90daynight20200228maskgal60r2"),
#      "pactf150reconvto90minus90daynight20200228maskgal60r2": cmbMap("./output/cmb_map/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150reconvto90_minus_f090_daynight_map.fits",
#         "./output/cmb_map/pact20200228_r2/" + "mask_full_foot_gal60_ps.fits",
#         "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits",
#         nu=150.e9, unitLatex=r'$\mu$K'
#         name="pactf150reconvto90minus90daynight20200228maskgal60r2"),
      #
      # TileC v1.2, reconvolved to 1.4' beam, combining BOSS N and D56
#      "tilecpactynocmb": cmbMap("./output/cmb_map/tilec_pact_ynocmb_v1.2.0/" + "tilec_map.fits",
#      "./output/cmb_map/tilec_pact_ynocmb_v1.2.0/" + "mask_full_foot_gal_ps.fits",
#      nu=0., unitLatex='',
#      name="tilecpactynocmb"),
#      "tilecpactyminusynocib": cmbMap("./output/cmb_map/tilec_pact_yminusynocib_v1.2.0/" + "diff_map.fits",
#      "./output/cmb_map/tilec_pact_yminusynocib_v1.2.0/" + "mask_full_foot_gal_ps.fits",
#      name="tilecpactyminusynocib"),
      #
      "tilecpacty": cmbMap("./output/cmb_map/tilec_pact_y_v1.2.0/" + "tilec_map.fits",
      "./output/cmb_map/tilec_pact_y_v1.2.0/" + "mask_full_foot_gal_ps.fits",
      nu=0., unitLatex='',
      name="tilecpacty"),
      "tilecpactynocib": cmbMap("./output/cmb_map/tilec_pact_ynocib_v1.2.0/" + "tilec_map.fits",
      "./output/cmb_map/tilec_pact_ynocib_v1.2.0/" + "mask_full_foot_gal_ps.fits",
      nu=0., unitLatex='',
      name="tilecpactynocib"),
      #"tilecpactcmbksz": cmbMap("./output/cmb_map/tilec_pact_cmbksz_v1.2.0/" + "tilec_reconv1.4_map.fits",
      #"./output/cmb_map/tilec_pact_cmbksz_v1.2.0/" + "mask_full_foot_gal_ps.fits",
      #nu=150.e9, unitLatex=r'$\mu$K',
      #name="tilecpactcmbksz"),
      #"tilecpactcmbksznoy": cmbMap("./output/cmb_map/tilec_pact_cmbksznoy_v1.2.0/" + "tilec_reconv2.4_map.fits",
      #"./output/cmb_map/tilec_pact_cmbksznoy_v1.2.0/" + "mask_full_foot_gal_ps.fits",
      #nu=150.e9, unitLatex=r'$\mu$K',
      #name="tilecpactcmbksznoy"),
      #
      # kSZ pipeline check
#      "pactf150daynight20200228maskgal60r2_minus_tilecpactcmbksz": cmbMap("./output/cmb_map/pactf150daynight20200228maskgal60r2_minus_tilec_pact_cmbksz/" + "diff_map.fits",
#         "./output/cmb_map/pactf150daynight20200228maskgal60r2_minus_tilec_pact_cmbksz/" + "mask_full_foot_gal_ps.fits",
#         nu=150.e9, unitLatex=r'$\mu$K',
#         name="pactf150daynight20200228maskgal60r2_minus_tilecpactcmbksz"),
      #
      # kSZ dust contamination test
#      "pactf150daynight20200228maskgal60r2_minus_tilecpactcmbksznocib": cmbMap("./output/cmb_map/pactf150daynight20200228maskgal60r2_minus_tilec_pact_cmbksznocib/" + "diff_map.fits",
#         "./output/cmb_map/pactf150daynight20200228maskgal60r2_minus_tilec_pact_cmbksznocib/" + "mask_full_foot_gal_ps.fits",
#         nu=150.e9, unitLatex=r'$\mu$K',
#         name="pactf150daynight20200228maskgal60r2_minus_tilecpactcmbksznocib"),
      #
      # tSZ pipeline (map) check
#      "pactf150daynight20200228maskgal60r2_minus_tilecpactymuk": cmbMap("./output/cmb_map/pactf150daynight20200228maskgal60r2_minus_tilec_pact_ymuk/" + "diff_map.fits",
#         "./output/cmb_map/pactf150daynight20200228maskgal60r2_minus_tilec_pact_ymuk/" + "mask_full_foot_gal_ps.fits",
#         nu=150.e9, unitLatex=r'$\mu$K',
#         name="pactf150daynight20200228maskgal60r2_minus_tilecpactymuk"),
      #
      # To have tau and y at the same TileC deproj beam
#      "pactf150daynight20200228maskgal60r2reconvtotilecdeproj": cmbMap("/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map_reconvtotilecdeproj.fits", 
#         "./output/cmb_map/pact20200228_r2/" + "mask_full_foot_gal60_ps.fits", 
#         nu=150.e9, unitLatex=r'$\mu$K',
#         name="pactf150daynight20200228maskgal60r2reconvtotilecdeproj"),
#      "pactf90daynight20200228maskgal60r2reconvtotilecdeproj": cmbMap("/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f090_daynight_map_reconvtotilecdeproj.fits", 
#         "./output/cmb_map/pact20200228_r2/" + "mask_full_foot_gal60_ps.fits", 
#         nu=90.e9, unitLatex=r'$\mu$K',
#         name="pactf90daynight20200228maskgal60r2reconvtotilecdeproj"),
      #
      # daynight VS night null test
#      "pactf150daynight20200228maskgal60r2_minus_night": cmbMap("./output/cmb_map/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_minus_night_map.fits",
#         "./output/cmb_map/pact20200228_r2/" + "mask_full_foot_gal60_ps.fits",
#         nu=150.e9, unitLatex=r'$\mu$K',
#         name="pactf150daynight20200228maskgal60r2_minus_night"),
      }

tStop = time()
print("took "+str(round((tStop-tStart)/60., 2))+" min")



###################################################################################


catalogCombi = {
      "pactf150daynight20200228maskgal60r2reconvtotilecdeproj": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      "pactf90daynight20200228maskgal60r2reconvtotilecdeproj": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      #
      "pactf150daynight20200228maskgal60r2": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      "pactf90daynight20200228maskgal60r2": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      "pactf150reconvto90minus90daynight20200228maskgal60r2": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      #
      "tilecpactynocmb": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      "tilecpactyminusynocib": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      #
      "tilecpacty": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      "tilecpactynocib": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      "tilecpactcmbksz": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      #
      "pactf150daynight20200228maskgal60r2_minus_tilecpactcmbksz": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      "pactf150daynight20200228maskgal60r2_minus_tilecpactymuk": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      "pactf150daynight20200228maskgal60r2_minus_tilecpactcmbksznocib": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      #
      "pactf150daynight20200228maskgal60r2_minus_night": ['lowz20200908full', 'lowz20200908mbin0', 'lowz20200908mbin1', 'lowz20200908mbin2'],
      }



###################################################################################
###################################################################################
# Compute all the stacked profiles

import thumbstack
reload(thumbstack)
from thumbstack import *


save = True


#for cmbMapKey in cmbMaps.keys():
#for cmbMapKey in cmbMaps.keys()[::-1]:
#for cmbMapKey in cmbMaps.keys()[:len(cmbMaps.keys())//2]:
for cmbMapKey in cmbMaps.keys()[len(cmbMaps.keys())//2:]:
#for cmbMapKey in ['pactf150daynight20200228maskgal60r2']:
#for cmbMapKey in ['pactf90daynight20200228maskgal60r2']:
#for cmbMapKey in ['pactf150daynight20200228maskgal60r2reconvtotilecdeproj', 'pactf90daynight20200228maskgal60r2reconvtotilecdeproj']:
#for cmbMapKey in ['tilecpacty', 'tilecpactynocib']:
   cmbMap = cmbMaps[cmbMapKey].map()
   cmbMask = cmbMaps[cmbMapKey].mask()
   cmbHit = cmbMaps[cmbMapKey].hit()
   cmbName = cmbMaps[cmbMapKey].name
   cmbNu = cmbMaps[cmbMapKey].nu
   cmbUnitLatex = cmbMaps[cmbMapKey].unitLatex
   print("Analyzing map "+cmbName)

   for catalogKey in catalogCombi[cmbMapKey]:
   #for catalogKey in ['lowz20200908mbin2']:
      catalog = catalogs[catalogKey]
      print("Analyzing catalog "+catalog.name)
      name = catalog.name + "_" + cmbName

      ts = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc, doMBins=False, doBootstrap=True, doVShuffle=False, cmbNu=cmbNu, cmbUnitLatex=cmbUnitLatex)


###################################################################################
###################################################################################
# Stacked cutout maps for PR

import thumbstack
reload(thumbstack)
from thumbstack import *


save = False


for cmbMapKey in cmbMaps.keys():
#for cmbMapKey in cmbMaps.keys()[::-1]:
#for cmbMapKey in ['pactf150daynight20200228maskgal60r2', 'pactf90daynight20200228maskgal60r2', 'tilecpactynocib', 'tilecpacty']:
#for cmbMapKey in ['pactf150daynight20200228maskgal60r2', 'pactf90daynight20200228maskgal60r2']:
#for cmbMapKey in ['tilecpactynocib', 'tilecpacty']:
   cmbMap = cmbMaps[cmbMapKey].map()
   cmbMask = cmbMaps[cmbMapKey].mask()
   cmbHit = cmbMaps[cmbMapKey].hit()
   cmbName = cmbMaps[cmbMapKey].name
   cmbNu = cmbMaps[cmbMapKey].nu
   cmbUnitLatex = cmbMaps[cmbMapKey].unitLatex
   print("Analyzing map "+cmbName)

   for catalogKey in catalogCombi[cmbMapKey]:
   #for catalogKey in ['lowz20200908full']:
      catalog = catalogs[catalogKey]
      print("Analyzing catalog "+catalog.name)
      name = catalog.name + "_" + cmbName

      print("C'est parti")

      try:
         ts = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=save, nProc=nProc, doMBins=False, doBootstrap=True, doVShuffle=False, doStackedMap=True, cmbNu=cmbNu, cmbUnitLatex=cmbUnitLatex)
      except: pass



###################################################################################
###################################################################################
# Joint cov 90-150 for PACT and PACT reconv to tilec deproj

'''
import thumbstack
reload(thumbstack)
from thumbstack import *



for catalogKey in ['cmass_kendrick', 'lowz_kendrick']:
#for catalogKey in ['cmass_kendrick']:
   catalog = catalogs[catalogKey]
   print("Analyzing catalog "+catalog.name)


   # PACT 90 - 150 reconv to tilec deproj
   ts = {}
   for freq in ['90', '150']:
      cmbMapKey = "pactf"+freq+"daynight20200228maskgal60r2reconvtotilecdeproj"
      cmbMap = cmbMaps[cmbMapKey].map()
      cmbMask = cmbMaps[cmbMapKey].mask()
      cmbHit = cmbMaps[cmbMapKey].hit()
      cmbName = cmbMaps[cmbMapKey].name
      cmbNu = cmbMaps[cmbMapKey].nu
      cmbUnitLatex = cmbMaps[cmbMapKey].unitLatex
      name = catalog.name + "_" + cmbName
      #
      ts[freq] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=False, nProc=nProc, doMBins=False, doBootstrap=False, doVShuffle=False, cmbNu=cmbNu, cmbUnitLatex=cmbUnitLatex)
   # compute the joint cov
   if True:
      ts['150'].saveAllCovBootstrapTwoStackedProfiles(ts['90'])
   ts['150'].plotAllCovTwoStackedProfiles(ts['90'])

   
   # PACT 90 - 150 
   ts = {}
   for freq in ['90', '150']:
      cmbMapKey = "pactf"+freq+"daynight20200228maskgal60r2"
      cmbMap = cmbMaps[cmbMapKey].map()
      cmbMask = cmbMaps[cmbMapKey].mask()
      cmbHit = cmbMaps[cmbMapKey].hit()
      cmbName = cmbMaps[cmbMapKey].name
      cmbNu = cmbMaps[cmbMapKey].nu
      cmbUnitLatex = cmbMaps[cmbMapKey].unitLatex
      name = catalog.name + "_" + cmbName
      #
      ts[freq] = ThumbStack(u, catalog, cmbMap, cmbMask, cmbHit, name, nameLong=None, save=False, nProc=nProc, doMBins=False, doBootstrap=False, doVShuffle=False, cmbNu=cmbNu, cmbUnitLatex=cmbUnitLatex)
   # compute the joint cov
   if True:
      ts['150'].saveAllCovBootstrapTwoStackedProfiles(ts['90'])
   ts['150'].plotAllCovTwoStackedProfiles(ts['90'])
   
'''



###################################################################################
###################################################################################
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
covKsz = {}
#
rTsz = {}
tsz = {}
sTsz = {}
covTsz = {}

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

      # read the kSZ cov
      path = pathThumb + name + '/cov_diskring_ksz_varweight_bootstrap.txt'
      if os.path.isfile(path):
         covKsz[name] = np.genfromtxt(path) * factor**2
      else:
         path = pathThumb + name + '/cov_diskring_ksz_uniformweight_bootstrap.txt'
         if os.path.isfile(path):
            covKsz[name] = np.genfromtxt(path) * factor**2


      # read the stacked tSZ profile
      try:
         data = np.genfromtxt(pathThumb + name + "/diskring_tsz_varweight_measured.txt")
      except:
         data = np.genfromtxt(pathThumb + name + "/diskring_tsz_uniformweight_measured.txt")
      rTsz[name] = data[:,0]
      tsz[name] = data[:,1] * factor
      sTsz[name] = data[:,2] * factor

      # read the tSZ cov
      path = pathThumb + name + '/cov_diskring_tsz_varweight_bootstrap.txt'
      if os.path.isfile(path):
         covTsz[name] = np.genfromtxt(path) * factor**2
      else:
         path = pathThumb + name + '/cov_diskring_tsz_uniformweight_bootstrap.txt'
         if os.path.isfile(path):
            covTsz[name] = np.genfromtxt(path) * factor**2


# read the joint 150-90 cov for ksz and tsz
for catalogKey in ['cmass_kendrick', 'lowz_kendrick']:
   catalog = catalogs[catalogKey]

   # PACT 150 - 90
   #
   path = pathThumb + catalog.name + "_pactf150daynight20200228maskgal60r2" + "/cov_diskring_ksz_varweight_joint_"+catalog.name+"_pactf150daynight20200228maskgal60r2_"+catalog.name+"_pactf90daynight20200228maskgal60r2_bootstrap.txt"
   covKsz[catalog.name+'_joint15090'] = np.genfromtxt(path) * factor**2
   #
   path = pathThumb + catalog.name + "_pactf150daynight20200228maskgal60r2" + "/cov_diskring_tsz_varweight_joint_"+catalog.name+"_pactf150daynight20200228maskgal60r2_"+catalog.name+"_pactf90daynight20200228maskgal60r2_bootstrap.txt"
   covTsz[catalog.name+'_joint15090'] = np.genfromtxt(path) * factor**2

   # PACT 150 - 90 reconv to tilec
   #
#   path = pathThumb + catalog.name + "_pactf150daynight20200228maskgal60r2reconvtotilecdeproj" + "/cov_diskring_ksz_varweight_joint_"+catalog.name+"_pactf150daynight20200228maskgal60r2reconvtotilecdeproj_"+catalog.name+"_pactf90daynight20200228maskgal60r2reconvtotilecdeproj_bootstrap.txt"
   path = pathThumb + catalog.name + "_pactf150daynight20200228maskgal60r2reconvtotilecdeproj" + "/cov_diskring_ksz_uniformweight_joint_"+catalog.name+"_pactf150daynight20200228maskgal60r2reconvtotilecdeproj_"+catalog.name+"_pactf90daynight20200228maskgal60r2reconvtotilecdeproj_bootstrap.txt"
   covKsz[catalog.name+'_joint15090reconvtotilecdeproj'] = np.genfromtxt(path) * factor**2
   #
#   path = pathThumb + catalog.name + "_pactf150daynight20200228maskgal60r2reconvtotilecdeproj" + "/cov_diskring_tsz_varweight_joint_"+catalog.name+"_pactf150daynight20200228maskgal60r2reconvtotilecdeproj_"+catalog.name+"_pactf90daynight20200228maskgal60r2reconvtotilecdeproj_bootstrap.txt"
   path = pathThumb + catalog.name + "_pactf150daynight20200228maskgal60r2reconvtotilecdeproj" + "/cov_diskring_tsz_uniformweight_joint_"+catalog.name+"_pactf150daynight20200228maskgal60r2reconvtotilecdeproj_"+catalog.name+"_pactf90daynight20200228maskgal60r2reconvtotilecdeproj_bootstrap.txt"
   covTsz[catalog.name+'_joint15090reconvtotilecdeproj'] = np.genfromtxt(path) * factor**2



# read the theory curves from Stefania
for catalogKey in ['cmass_kendrick', 'lowz_kendrick']:
   catalog = catalogs[catalogKey]

   path = './input/stefania/theory_curves/'
   
   # NFW profiles at 90 and 150 GHz
   data = np.genfromtxt(path+'nfw_'+catalogKey+'_mw_trunc1rvir.txt')
   ksz[catalog.name+'_150_nfw'] = data[:,1] * factor
   ksz[catalog.name+'_90_nfw'] = data[:,2] * factor

   # Best fit theory curves
   if catalog.name=='cmass_kendrick':
      #
      data = np.genfromtxt(path+'cmass_kendrick_ksz_best_150_90.txt')
      ksz[catalog.name+'_150_theory'] = data[:,1] * factor
      ksz[catalog.name+'_90_theory'] = data[:,2] * factor
      #
      data = np.genfromtxt(path+'cmass_kendrick_tsz+dust_best_150_90.txt')
      tsz[catalog.name+'_150_theory'] = data[:,1] * factor
      tsz[catalog.name+'_90_theory'] = data[:,2] * factor
      #
      data = np.genfromtxt(path+'cmass_kendrick_y_nocib_best.txt')
      tsz[catalog.name+'_ynocib_theory'] = data[:,1] * factor




# Compare to Schaan+16 measurement
#
# CMASS K
path = "./input/ksz_schaan+16/diagonal_cov_paper/alpha_ksz_kendrick.txt"
data = np.genfromtxt(path)
rKsz['cmass_kendrick_150_schaan+16'] = data[:,0] # [arcmin]
alphaK = data[:,1]   # [dimless]
sAlphaK = data[:,2]  # [dimless]
# convert to the new measurement unit
cmassK = catalogs['cmass_kendrick']
ksz['cmass_kendrick_150_schaan+16'] = alphaK * np.std(cmassK.vR) / 3.e5 / cmassK.rV * 2.726e6 * np.mean(cmassK.integratedTau) * (180.*60./np.pi)**2
sKsz['cmass_kendrick_150_schaan+16'] = sAlphaK * np.std(cmassK.vR) / 3.e5  / cmassK.rV * 2.726e6 * np.mean(cmassK.integratedTau) * (180.*60./np.pi)**2
#
# CMASS M
path = "./input/ksz_schaan+16/diagonal_cov_paper/alpha_ksz_mariana.txt"
data = np.genfromtxt(path)
rKsz['cmass_mariana_150_schaan+16'] = data[:,0] # [arcmin]
alphaM = data[:,1]   # [dimless]
sAlphaM = data[:,2]  # [dimless]
# convert to the new measurement unit
cmassM = catalogs['cmass_mariana']
ksz['cmass_mariana_150_schaan+16'] = alphaM * np.std(cmassM.vR) / 3.e5 / cmassM.rV * 2.726e6 * np.mean(cmassM.integratedTau) * (180.*60./np.pi)**2
sKsz['cmass_mariana_150_schaan+16'] = sAlphaM * np.std(cmassM.vR) / 3.e5 / cmassM.rV * 2.726e6 * np.mean(cmassM.integratedTau) * (180.*60./np.pi)**2




# read the stacks on mock GRFs, to compare
pathMockGRF = "/global/cscratch1/sd/eschaan/project_ucsc/code/thumbstack/output/cmb_map/mocks_grf_planck_act_coadd_2019_03_11/"
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
data = np.genfromtxt(pathThumb + "cmass_kendrick_pactf150daynight20200228maskgal60r2/" + "diskring_ksz_varweight_vshufflemean.txt")
rKsz150VShuffleMean = data[:,0]
ksz150VShuffleMean = data[:,1] * factor
sKsz150VShuffleMean = data[:,2] * factor
'''


###################################################################################
###################################################################################
###################################################################################
# Compute significances

'''
for catalogKey in ['cmass_kendrick', 'lowz_kendrick']:
   
   # kSZ
   d = np.concatenate((ksz[catalogKey+'_pactf150daynight20200228maskgal60r2'],
                       ksz[catalogKey+'_pactf90daynight20200228maskgal60r2']))
   cov = covKsz[catalogKey+'_joint15090']
   tNFW = np.concatenate((ksz[catalogKey+'_150_nfw'],
                       ksz[catalogKey+'_90_nfw']))
   if catalogKey=='cmass_kendrick':
      t = np.concatenate((ksz[catalogKey+'_150_theory'],
                          ksz[catalogKey+'_90_theory']))
   else:
      t = np.zeros_like(d)

   print(catalogKey+", joint ksz 150-90, best fit VS null")
   computeSnr(d, t, cov)
   print(catalogKey+", joint ksz 150-90, best fit VS NFW")
   computeSnr(d - tNFW, t - tNFW, cov)
   print(catalogKey+", joint ksz 150-90, null VS NFW")
   computeSnr(d - tNFW, 0. - tNFW, cov)


   # tSZ
   d = tsz[catalogKey+'_tilecpactynocib']
   cov = covTsz[catalogKey+'_tilecpactynocib']
   if catalogKey=='cmass_kendrick':
      t = tsz[catalogKey+'_ynocib_theory']
   else:
      t = np.zeros_like(d)   

   print(catalogKey+", tsz on tilec y no cib")
   computeSnr(d, t, cov)

   
   # tSZ  + dust
   print(catalogKey+", joint tsz+dust 150-90")
   # joint tsz + dust,
   # discard the first data point
   d = np.concatenate((tsz[catalogKey+'_pactf150daynight20200228maskgal60r2'][1:],
      tsz[catalogKey+'_pactf90daynight20200228maskgal60r2'][1:]))
   I = np.concatenate((range(1,9), range(9+1, 9+9)))
   J = np.ix_(I,I)
   cov = covTsz[catalogKey+'_joint15090'][J]
   if catalogKey=='cmass_kendrick':
      t = np.concatenate((tsz[catalogKey+'_150_theory'][1:],
         tsz[catalogKey+'_90_theory'][1:]))
   else:
      t = np.zeros_like(d)
   computeSnr(d, t, cov)
'''






###################################################################################
###################################################################################
###################################################################################
###################################################################################
# Generate all the plots
'''

for catalogKey in ['cmass_kendrick', 'lowz_kendrick']:
   catalog = catalogs[catalogKey]

   # true velocity RMS for CMASS K
   zMean = catalog.Z.mean()   # Msun
   vRms = u.v1dRms(0., zMean, W3d_sth) # [km/s]
   tauToKsz = 2.726e6 * (vRms/3.e5)
   # virial radius, before convolving with the beams
   mVirMean = catalog.Mvir.mean()   # [Msun]. Will need to be converted to [Msun/h] below 
   rVir = u.frvir(mVirMean * u.bg.h, zMean)   # [Mpc/h]
   thetaVir = rVir / u.bg.comoving_distance(zMean) * (180.*60./np.pi) # [arcmin]

   print catalogKey
   print "mean z =", zMean
   print "mean chi =", u.bg.comoving_distance(zMean), "Mpc/h"
   print "mean mass =", mVirMean, "Msun"
   print "virial radius =", rVir, "Mpc/h"
   print "virial angle =", thetaVir, "arcmin"
   print "RMS 1d velocity =", vRms, "km/s"

   if catalogKey=='cmass_kendrick':
      catalogTitle = 'CMASS'
      fmt = 'o'
   elif catalogKey=='lowz_kendrick':
      catalogTitle = 'LOWZ'
      fmt = 'o:'#'o--'

   rAp = rKsz[catalogKey+'_pactf150daynight20200228maskgal60r2']


   ###################################################################################
   # Summary ksz plots at 150 and 90


   # PACT 150
   ksz150 = ksz[catalogKey+'_pactf150daynight20200228maskgal60r2']
   sKsz150 = sKsz[catalogKey+'_pactf150daynight20200228maskgal60r2']
   #
   # PACT 90
   ksz90 = ksz[catalogKey+'_pactf90daynight20200228maskgal60r2']
   sKsz90 = sKsz[catalogKey+'_pactf90daynight20200228maskgal60r2']

   # best fit theory curves
   if catalogKey=='cmass_kendrick':
      ksz150Th = ksz[catalog.name+'_150_theory']
   ksz150 = ksz[catalogKey+'_pactf150daynight20200228maskgal60r2']
   sKsz150 = sKsz[catalogKey+'_pactf150daynight20200228maskgal60r2']
   #
   # PACT 90
   ksz90 = ksz[catalogKey+'_pactf90daynight20200228maskgal60r2']
   sKsz90 = sKsz[catalogKey+'_pactf90daynight20200228maskgal60r2']

   # best fit theory curves
   if catalogKey=='cmass_kendrick':
      ksz150Th = ksz[catalog.name+'_150_theory']
      ksz90Th = ksz[catalog.name+'_90_theory']

   # NFW curves
   ksz150NFW = ksz[catalog.name+'_150_nfw']
   ksz90NFW = ksz[catalog.name+'_90_nfw']

   # NFW curves
   ksz150NFW = ksz[catalog.name+'_150_nfw']
   ksz90NFW = ksz[catalog.name+'_90_nfw']



   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   ax.minorticks_on()
   #
   ax.axhline(0., c='k', lw=1)
   #
   # virial radius
   ax.axvline(np.sqrt(thetaVir**2 + 1.3**2/(8.*np.log(2.))), color='blue', alpha=0.1)
   ax.axvline(np.sqrt(thetaVir**2 + 2.1**2/(8.*np.log(2.))), color='purple', alpha=0.1)
   #
   # data
   ax.errorbar(rAp, ksz150, sKsz150, fmt=fmt, c='blue', label='150 GHz DR5')
   ax.errorbar(rAp + 0.05, ksz90, sKsz90, fmt=fmt, c='purple', label='98 GHz DR5')
   #
   # best fit theory curves
   if catalogKey=='cmass_kendrick':
      ax.plot(rAp, ksz150Th, 'blue', label=r'Joint best fit profile')
      ax.plot(rAp, ksz90Th, 'purple')
      #
      # NFW profiles 
      ax.plot(rAp, ksz150NFW, ls='--', c='blue', label=r'NFW')
      ax.plot(rAp, ksz90NFW, ls='--', c='purple')
   #
   # Gaussian profiles to guide the eye
   fwhm = np.array([1.3, 2.1, 6.])   # arcmin
   sigma = fwhm / np.sqrt(8. * np.log(2.))
   for s in sigma:
      #y = ftheoryGaussianProfile(s, filterType='diskring')
      y = (1. - np.exp(-0.5*rAp**2/s**2))**2
      # normalize to the last data point at 150 GHz
      y *= ksz150[-1]
      ax.plot(rAp, y, 'gray', alpha=0.3)
   ax.plot([], [], 'gray', alpha=0.3, label=r'Gaussian profiles')
   #
   ax.legend(loc=4, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(catalogTitle + r' kSZ profile', x=0.5, y=1.25)
   ax.set_yscale('log', nonposy='clip')
   ax.set_ylim((0.03, 70.))
   #
   # make extra abscissa with disk comoving size in Mpc/h
   ax2 = ax.twiny()
   ax2.minorticks_on()
   ticks = ax.get_xticks()
   ax2.set_xticks(ticks)
   newticks = np.array(ticks) * np.pi/(180.*60.)*u.bg.comoving_distance(zMean)  # disk radius in Mpc/h
   newticks = np.round(newticks, 2)
   ax2.set_xticklabels(newticks)
   ax2.set_xlim(ax.get_xlim())
   ax2.set_xlabel(r'Comoving radius [Mpc/h] at $z=$'+str(round(zMean,2)), fontsize=20)
   ax2.xaxis.set_label_coords(0.5, 1.15)
   #
   # extra ordinate to convert kSZ to tau
   ax3 = ax.twinx()
   ax3.minorticks_on()
   ax3.set_yscale('log', nonposy='clip')
   ylim = ax.get_ylim()
   ylim = np.array(ylim) / tauToKsz
   ax3.set_ylim(ylim)
   ax3.set_ylabel(r'Integrated $\tau_\text{CAP}$ [arcmin$^2$]', fontsize=20)
   #
   path = pathFig+"summary_ksz_150_90_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()



   ###################################################################################
   # ksz: comparison plot with Schaan+16

   if catalogKey=='cmass_kendrick':

      # Schaan+16 ksz data, converted to the new unit
      rK = rKsz['cmass_kendrick_150_schaan+16']
      kszKS16 = ksz['cmass_kendrick_150_schaan+16']
      sKszKS16 = sKsz['cmass_kendrick_150_schaan+16']
      #
      rM = rKsz['cmass_mariana_150_schaan+16']
      kszMS16 = ksz['cmass_mariana_150_schaan+16']
      sKszMS16 = sKsz['cmass_mariana_150_schaan+16']



      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      ax.minorticks_on()
      #
      ax.axhline(0., c='k', lw=1)
      #
      # virial radius
      ax.axvline(np.sqrt(thetaVir**2 + 1.3**2/(8.*np.log(2.))), color='blue', alpha=0.1)
      #
      # data
      ax.errorbar(rAp, ksz150, sKsz150, fmt=fmt, c='k', label='This work 150 GHz CMASS K')
      #
      # comparison with Schaan+16
      ax.errorbar(rM, kszMS16, yerr=sKszMS16, c='blue', label=r'Schaan+16 150 GHz CMASS M')
      ax.errorbar(rK + 0.05, kszKS16, yerr=sKszKS16, c='blue', alpha=0.3, label=r'Schaan+16 150 GHz CMASS K')
      #
#      # Theory curves if available
#      if catalogKey=='cmass_kendrick':
#         # best fit theory curves
#         ax.plot(rAp, ksz150Th, 'blue', label=r'Joint best fit profile')
#         ax.plot(rAp, ksz90Th, 'purple')
#         # NFW profiles 
#         ax.plot(rAp, ksz150NFW, ls='--', c='blue', label=r'NFW')
#         ax.plot(rAp, ksz90NFW, ls='--', c='purple')
      #
      ax.legend(loc=4, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
      ax.set_title(catalogTitle + r' kSZ profile', x=0.5, y=1.25)
      ax.set_yscale('log', nonposy='clip')
      ax.set_ylim((0.03, 40.))
      #
      # make extra abscissa with disk comoving size in Mpc/h
      ax2 = ax.twiny()
      ax2.minorticks_on()
      ticks = ax.get_xticks()
      ax2.set_xticks(ticks)
      newticks = np.array(ticks) * np.pi/(180.*60.)*u.bg.comoving_distance(zMean)  # disk radius in Mpc/h
      newticks = np.round(newticks, 2)
      ax2.set_xticklabels(newticks)
      ax2.set_xlim(ax.get_xlim())
      ax2.set_xlabel(r'Comoving radius [Mpc/h] at $z=$'+str(round(zMean,2)), fontsize=20)
      ax2.xaxis.set_label_coords(0.5, 1.15)
      #
      # extra ordinate to convert kSZ to tau
      ax3 = ax.twinx()
      ax3.minorticks_on()
      ax3.set_yscale('log', nonposy='clip')
      ylim = ax.get_ylim()
      ylim = np.array(ylim) / tauToKsz
      ax3.set_ylim(ylim)
      ax3.set_ylabel(r'Integrated $\tau_\text{CAP}$ [arcmin$^2$]', fontsize=20)
      #
      path = pathFig+"summary_ksz_150_90_"+catalogKey+"_vs_schaan+16.pdf"
      fig.savefig(path, bbox_inches='tight')
      #plt.show()
      fig.clf()



   ###################################################################################
   # summary tSZ + dust at 150 and 90

   # PACT 150
   tsz150 = tsz[catalogKey+'_pactf150daynight20200228maskgal60r2']
   sTsz150 = sTsz[catalogKey+'_pactf150daynight20200228maskgal60r2']
   #
   # PACT 90
   tsz90 = tsz[catalogKey+'_pactf90daynight20200228maskgal60r2']
   sTsz90 = sTsz[catalogKey+'_pactf90daynight20200228maskgal60r2']

   if catalogKey=='cmass_kendrick':
      # best fit theory curves
      tsz150Th = tsz[catalog.name+'_150_theory']
      tsz90Th = tsz[catalog.name+'_90_theory']


   # tSZ + dust plot at 150 and 90
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   ax.minorticks_on()
   #
   ax.axhline(0., c='k', lw=1)
   #
   # virial radius
   ax.axvline(np.sqrt(thetaVir**2 + 1.3**2/(8.*np.log(2.))), color='blue', alpha=0.1)
   ax.axvline(np.sqrt(thetaVir**2 + 2.1**2/(8.*np.log(2.))), color='purple', alpha=0.1)
   #
   ax.errorbar(rAp, tsz150, sTsz150, fmt=fmt, c='blue', label='150 GHz DR5')
   ax.errorbar(rAp + 0.05, tsz90, sTsz90, fmt=fmt, c='purple', label='98 GHz DR5')
   #
   # Theory curves if available
   if catalogKey=='cmass_kendrick':
      ax.plot(rAp, tsz150Th, 'blue', label=r'Joint best fit profile')
      ax.plot(rAp, tsz90Th, 'purple')
   #
   ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_{\text{tSZ} + \text{dust}}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(catalogTitle + r' tSZ + dust profile', x=0.5, y=1.25)
   ax.set_yscale('symlog')
   #
   # make extra abscissa with disk comoving size in Mpc/h
   ax2 = ax.twiny()
   ax2.minorticks_on()
   ticks = ax.get_xticks()
   ax2.set_xticks(ticks)
   newticks = np.array(ticks) * np.pi/(180.*60.)*u.bg.comoving_distance(zMean)  # disk radius in Mpc/h
   newticks = np.round(newticks, 2)
   ax2.set_xticklabels(newticks)
   ax2.set_xlim(ax.get_xlim())
   ax2.set_xlabel(r'Comoving radius [Mpc/h] at $z=$'+str(round(zMean,2)), fontsize=20)
   ax2.xaxis.set_label_coords(0.5, 1.15)
   #
   path = pathFig+"summary_tsz_150_90_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()



   ###################################################################################
   # Summary tSZ on TileC y no CIB 


   # TileC y no Cib
   tszYNoCib = tsz[catalogKey+'_tilecpactynocib'] * yTomuK150
   sTszYNoCib = sTsz[catalogKey+'_tilecpactynocib'] * yTomuK150

   # best fit theory curves
   if catalogKey=='cmass_kendrick':
      tszYNoCibTh = tsz[catalog.name+'_ynocib_theory'] * yTomuK150

   
   # tSZ-only from the TileC y no CIB map
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   ax.minorticks_on()
   #
   ax.axhline(0., c='k', lw=1)
   #
   # virial radius
   ax.axvline(np.sqrt(thetaVir**2 + 2.4**2/(8.*np.log(2.))), color='r', alpha=0.1)
   #
   # Tilec y no CIB
   ax.errorbar(rAp, tszYNoCib, yerr=sTszYNoCib, fmt=fmt, c='r', label='ILC y no CIB DR4')
   #
   # Theory curves if available
   if catalogKey=='cmass_kendrick':
      ax.plot(rAp, tszYNoCibTh, 'r')
   #
   ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_{\text{tSZ}}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(catalogTitle + r' tSZ profile', x=0.5, y=1.25)
   ax.set_yscale('symlog', linthreshy=1e-1, linscaley=0.1)
   #ax.set_ylim((0., 2.))
   #
   # make extra abscissa with disk comoving size in Mpc/h
   ax2 = ax.twiny()
   ax2.minorticks_on()
   ticks = ax.get_xticks()
   ax2.set_xticks(ticks)
   newticks = np.array(ticks) * np.pi/(180.*60.)*u.bg.comoving_distance(zMean)  # disk radius in Mpc/h
   newticks = np.round(newticks, 2)
   ax2.set_xticklabels(newticks)
   ax2.set_xlim(ax.get_xlim())
   ax2.set_xlabel(r'Comoving radius [Mpc/h] at $z=$'+str(round(zMean,2)), fontsize=20)
   ax2.xaxis.set_label_coords(0.5, 1.15)
   #
   # extra ordinate to convert tSZ to y units
   ax3 = ax.twinx()
   ax3.minorticks_on()
   ax3.set_yscale('symlog', linthreshy=1e-1 / np.abs(yTomuK150), linscaley=0.1)
   ylim = ax.get_ylim()
   ylim = np.array(ylim) / np.abs(yTomuK150)
   ax3.set_ylim(ylim)
   ax3.set_ylabel(r'Integrated $-y_\text{CAP}$ [arcmin$^2$]', fontsize=20)
   #
   path = pathFig+"summary_tsz_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()

   
   
   ###################################################################################
   # Comparison: tSZ + dust plot on 150, 90, TileC


   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   ax.minorticks_on()
   #
   ax.axhline(0., c='k', lw=1)
   #
   # virial radius
   ax.axvline(np.sqrt(thetaVir**2 + 1.3**2/(8.*np.log(2.))), color='blue', alpha=0.1)
   ax.axvline(np.sqrt(thetaVir**2 + 2.1**2/(8.*np.log(2.))), color='purple', alpha=0.1)
   ax.axvline(np.sqrt(thetaVir**2 + 2.4**2/(8.*np.log(2.))), color='r', alpha=0.1)
   #
   # Tilec y no CIB
   ax.errorbar(rAp, tszYNoCib, yerr=sTszYNoCib, fmt=fmt, c='r', label='ILC y no CIB DR4')
   # PACT 150
   ax.errorbar(rAp, tsz150, yerr=sTsz150, fmt='-', c='blue', label='150 GHz DR5')
   # PACT 90
   ax.errorbar(rAp, tsz90, yerr=sTsz90, fmt='-', c='purple', label='98 GHz DR5')
   #
   ax.legend(loc=3, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_{\text{tSZ}}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(catalogTitle + r' tSZ profile', x=0.5, y=1.25)
   ax.set_yscale('symlog')
   #ax.set_ylim((0., 2.))
   #
   # make extra abscissa with disk comoving size in Mpc/h
   ax2 = ax.twiny()
   ax2.minorticks_on()
   ticks = ax.get_xticks()
   ax2.set_xticks(ticks)
   newticks = np.array(ticks) * np.pi/(180.*60.)*u.bg.comoving_distance(zMean)  # disk radius in Mpc/h
   newticks = np.round(newticks, 2)
   ax2.set_xticklabels(newticks)
   ax2.set_xlim(ax.get_xlim())
   ax2.set_xlabel(r'Comoving radius [Mpc/h] at $z=$'+str(round(zMean,2)), fontsize=20)
   ax2.xaxis.set_label_coords(0.5, 1.15)
   #
   #
   path = pathFig+"comparison_tsz_150_90_tilec_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()



   ###################################################################################
   ###################################################################################
   # kSZ null tests

   
   # fiducial uncertainty
   ksz150 = ksz[catalogKey+'_pactf150daynight20200228maskgal60r2']
   sKsz150 = sKsz[catalogKey+'_pactf150daynight20200228maskgal60r2']
   # Mariana's velocities
   if catalogKey=='cmass_kendrick':
      catalogMariana = 'cmass_mariana'
      ksz150Mariana = ksz[catalogMariana+'_pactf150daynight20200228maskgal60r2']
      sKsz150Mariana = sKsz[catalogMariana+'_pactf150daynight20200228maskgal60r2']
   #
   # 150 - tilec cmb, as a consistency check
   ksz150MinusTilecCmb = ksz[catalogKey+'_pactf150daynight20200228maskgal60r2_minus_tilecpactcmbksz']
   sKsz150MinusTilecCmb = sKsz[catalogKey+'_pactf150daynight20200228maskgal60r2_minus_tilecpactcmbksz']
   #
   # 150 reconv to 90 minus 90, to check consistency
   ksz150Reconv90Minus90 = ksz[catalogKey+'_pactf150reconvto90minus90daynight20200228maskgal60r2']
   sKsz150Reconv90Minus90 = sKsz[catalogKey+'_pactf150reconvto90minus90daynight20200228maskgal60r2']
   #
   # tilec y no cmb, to check for tSZ contamination
   kszYNoCmb = ksz[catalogKey+'_tilecpactynocmb'] * yTomuK150
   sKszYNoCmb = sKsz[catalogKey+'_tilecpactynocmb'] * yTomuK150
   #
   # 150 - tilec cmb no cib, to check for dust contamination
   ksz150MinusCmbNoCib = ksz[catalogKey+'_pactf150daynight20200228maskgal60r2_minus_tilecpactcmbksznocib']
   sKsz150MinusCmbNoCib = sKsz[catalogKey+'_pactf150daynight20200228maskgal60r2_minus_tilecpactcmbksznocib']
   #
   # 150 daynight minus night, for null test
   if catalogKey=='cmass_kendrick':
      ksz150DaynightMinusNight = ksz[catalogKey+'_pactf150daynight20200228maskgal60r2_minus_night']
      sKsz150DaynightMinusNight = sKsz[catalogKey+'_pactf150daynight20200228maskgal60r2_minus_night']




   # kSZ pipeline null tests
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   ax.minorticks_on()
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
   ax.errorbar(rAp + 0.05, meanStackedKszGRF, yerr=sStackedKszGRF, fmt='--', label=r'mean of '+str(nMocks)+' mocks')
   #
   if catalogKey=='cmass_kendrick':
      # Mariana - Kendrick
      #ax.errorbar(rAp, ksz150MKDiff, yerr=sKsz150MKDiff, fmt='-', c='r', label=r'$v_\text{Mariana} - v_\text{Kendrick}$')
      #ax.errorbar(rAp, ksz150, yerr=sKsz150, fmt='-', label='K')
      #ax.errorbar(rAp, ksz150Mariana, yerr=sKsz150Mariana, fmt='-', label='M')
      #ax.errorbar(rAp + 0.05, (ksz150-ksz150Mariana), yerr=sKsz150, fmt='-', label=r'$v_\text{Kendrick} - v_\text{Mariana}$')
      ax.plot(rAp + 0.1, (ksz150-ksz150Mariana), '-', label=r'$v_\text{fiducial} - v_\text{M}$')
   #
   # 150 - tilec cmb
   ax.errorbar(rAp + 0.15, ksz150MinusTilecCmb, yerr=sKsz150MinusTilecCmb, fmt='-', label='150 - ILC CMB/kSZ')
   #
   # 150 reconv to 90 minus 90
   ax.errorbar(rAp + 0.2, ksz150Reconv90Minus90, yerr=sKsz150Reconv90Minus90, fmt='-', label='150 - 98')
   #
   # 150 daynight minus night
   if catalogKey=='cmass_kendrick':
      ax.errorbar(rAp + 0.25, ksz150DaynightMinusNight, yerr=sKsz150DaynightMinusNight, fmt='-', label='150 daynight - night')
   #
   ax.set_ylim((-10., 15.))
   ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(catalogTitle + r' kSZ pipeline null tests', x=0.5, y=1.25)
   ax.set_ylim((-6., 6.))
   #
   # make extra abscissa with disk comoving size in Mpc/h
   ax2 = ax.twiny()
   ax2.minorticks_on()
   ticks = ax.get_xticks()
   ax2.set_xticks(ticks)
   newticks = np.array(ticks) * np.pi/(180.*60.)*u.bg.comoving_distance(catalogs[catalogKey].Z.mean())  # disk radius in Mpc/h
   newticks = np.round(newticks, 2)
   ax2.set_xticklabels(newticks)
   ax2.set_xlim(ax.get_xlim())
   ax2.set_xlabel(r'Comoving radius [Mpc/h] at $z=$'+str(round(catalogs[catalogKey].Z.mean(),2)), fontsize=20)
   ax2.xaxis.set_label_coords(0.5, 1.15)
   #
   path = pathFig + "pipenulltests_ksz_150_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()




   # kSZ foreground null tests
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   ax.minorticks_on()
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
   ax.errorbar(rAp, kszYNoCmb, yerr=sKszYNoCmb, label=r'ILC y no CMB')
   #
   # cmbksz no cib, to check for dust
   ax.errorbar(rAp + 0.05, ksz150MinusCmbNoCib, yerr=sKsz150MinusCmbNoCib, fmt='-', label=r'150 - ILC CMB/kSZ no CIB')
   #
   ax.set_ylim((-10., 15.))
   ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(catalogTitle + r' kSZ foreground null tests', x=0.5, y=1.25)
   #ax.set_ylim((-2., 2.))
   #
   # make extra abscissa with disk comoving size in Mpc/h
   ax2 = ax.twiny()
   ax2.minorticks_on()
   ticks = ax.get_xticks()
   ax2.set_xticks(ticks)
   newticks = np.array(ticks) * np.pi/(180.*60.)*u.bg.comoving_distance(catalogs[catalogKey].Z.mean())  # disk radius in Mpc/h
   newticks = np.round(newticks, 2)
   ax2.set_xticklabels(newticks)
   ax2.set_xlim(ax.get_xlim())
   ax2.set_xlabel(r'Comoving radius [Mpc/h] at $z=$'+str(round(catalogs[catalogKey].Z.mean(),2)), fontsize=20)
   ax2.xaxis.set_label_coords(0.5, 1.15)
   #
   path = pathFig + "fgnulltests_ksz_150_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()



   ###################################################################################
   ###################################################################################
   # tSZ null tests

   # fiducial uncertainty
   sTsz150 = sTsz[catalogKey+'_pactf150daynight20200228maskgal60r2'] 
   #
   # 150 - y, to check for map consistency
   tsz150MinusY = tsz[catalogKey+'_pactf150daynight20200228maskgal60r2_minus_tilecpactymuk']
   sTsz150MinusY = sTsz[catalogKey+'_pactf150daynight20200228maskgal60r2_minus_tilecpactymuk']
   #
   # y - y no CIB
   tszYMinusYNoCib = tsz[catalogKey+'_tilecpactyminusynocib'] * yTomuK150
   sTszYMinusYNoCib = sTsz[catalogKey+'_tilecpactyminusynocib'] * yTomuK150
   #
   # 150' - 90, after rescaling 90 to null tSZ
   # in order to check for the dust contamination
#   tsz150Reconv90Minus90NoY = tsz[catalogKey+'_pactf150reconvto90minus90noydaynight20200228maskgal60r2']
#   sTsz150Reconv90Minus90NoY = sTsz[catalogKey+'_pactf150reconvto90minus90noydaynight20200228maskgal60r2']



   # tSZ pipeline test
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   ax.minorticks_on()
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
   ax.errorbar(rAp, tsz150MinusY, yerr=sTsz150MinusY, fmt='--', label=r'150 - ILC y ')
   #
   ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_\text{tSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(catalogTitle + r' tSZ pipeline null tests', x=0.5, y=1.25)
   ax.set_ylim((-6., 6.))
   #
   # make extra abscissa with disk comoving size in Mpc/h
   ax2 = ax.twiny()
   ax2.minorticks_on()
   ticks = ax.get_xticks()
   ax2.set_xticks(ticks)
   newticks = np.array(ticks) * np.pi/(180.*60.)*u.bg.comoving_distance(catalogs[catalogKey].Z.mean())  # disk radius in Mpc/h
   newticks = np.round(newticks, 2)
   ax2.set_xticklabels(newticks)
   ax2.set_xlim(ax.get_xlim())
   ax2.set_xlabel(r'Comoving radius [Mpc/h] at $z=$'+str(round(catalogs[catalogKey].Z.mean(),2)), fontsize=20)
   ax2.xaxis.set_label_coords(0.5, 1.15)
   #
   path = pathFig + "pipenulltests_tsz_150_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()



   # dust contamination to tSZ
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   ax.minorticks_on()
   #
   ax.axhline(0., c='k', lw=1)
   #
   # Uncertainty band
   ax.fill_between(rAp, - sTsz150, sTsz150, edgecolor='', facecolor='gray', alpha=0.5, label=r'statistical error')
   #
   # 150 - tilec y
   ax.errorbar(rAp, tsz150MinusY, yerr=sTsz150MinusY, fmt='--', label=r'150 - ILC y')
   #
   # y - y no CIB
   ax.errorbar(rAp, tszYMinusYNoCib, yerr=sTszYMinusYNoCib, fmt='--', label=r'ILC y - y no CIB')
   #
#   # 150' - 90 rescaled to null y
#   ax.errorbar(rAp, tsz150Reconv90Minus90NoY, yerr=sTsz150Reconv90Minus90NoY, fmt='--', label=r"150\' - 90 no y")
   #
   ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_\text{dust}$ [$\mu K\cdot\text{arcmin}^2$]')
   ax.set_title(catalogTitle + r' Dust emission', x=0.5, y=1.25)
   ax.set_ylim((-6., 6.))
   #
   # make extra abscissa with disk comoving size in Mpc/h
   ax2 = ax.twiny()
   ax2.minorticks_on()
   ticks = ax.get_xticks()
   ax2.set_xticks(ticks)
   newticks = np.array(ticks) * np.pi/(180.*60.)*u.bg.comoving_distance(catalogs[catalogKey].Z.mean())  # disk radius in Mpc/h
   newticks = np.round(newticks, 2)
   ax2.set_xticklabels(newticks)
   ax2.set_xlim(ax.get_xlim())
   ax2.set_xlabel(r'Comoving radius [Mpc/h] at $z=$'+str(round(catalogs[catalogKey].Z.mean(),2)), fontsize=20)
   ax2.xaxis.set_label_coords(0.5, 1.15)
   #
   path = pathFig + "fgnulltests_tsz_150_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()



   ###################################################################################
   ###################################################################################
   # Electron temperature


   # y
   y = tsz[catalogKey+'_tilecpactynocib'] * (180.*60./np.pi)**2   # [arcmin^2]
   sY = sTsz[catalogKey+'_tilecpactynocib'] * (180.*60./np.pi)**2   # [arcmin^2]

   # units for conversion
   me = 9.10938291e-31  # electron mass [kg]
   c = 299792458. # [m/s]
   kB = 1.38e-23  # [SI] = [m2 kg s-2 K-1]

   # compute Virial temperature
   # from van den Bosch's lecture note slides
   tVir = 3.6e5 * (mVirMean / 1.e11)**(2./3.) # [K]

   # From 150 GHz
   # convert ksz to integrated tau
   inTau150 = ksz[catalogKey+'_pactf150daynight20200228maskgal60r2reconvtotilecdeproj'].copy()  # [muK*sr]
   inTau150 *= (180.*60./np.pi)**2   # [muK * arcmin^2
   inTau150 /= tauToKsz # [arcmin^2]
   # uncertainty on integrated tau
   sInTau150 = sKsz[catalogKey+'_pactf150daynight20200228maskgal60r2reconvtotilecdeproj'].copy()  # [muK*sr]
   sInTau150 *= (180.*60./np.pi)**2   # [muK * arcmin^2]
   sInTau150 /= tauToKsz  # [arcmin^2]
   # from tau and y, get temperature
   tE150 = y / inTau150 # = kB*Te / (me*c^2) [dimless]
   # convert electron temperature to K
   tE150 *= me * c**2 / kB
   # uncertainty on temperature: assume y and tau to be uncorrelated
   sTE150 = tE150 * np.sqrt( (sY/y)**2 + (sInTau150/inTau150)**2 )


   # From 90 GHz
   # convert ksz to integrated tau
   inTau90 = ksz[catalogKey+'_pactf90daynight20200228maskgal60r2reconvtotilecdeproj'].copy()  # [muK*sr]
   inTau90 *= (180.*60./np.pi)**2   # [muK * arcmin^2
   inTau90 /= 2.726e6 * (vRms/3.e5)
   # uncertainty on integrated tau
   sInTau90 = sKsz[catalogKey+'_pactf90daynight20200228maskgal60r2reconvtotilecdeproj'].copy()  # [muK*sr]
   sInTau90 *= (180.*60./np.pi)**2   # [muK * arcmin^2]
   sInTau90 /= 2.726e6 * (vRms/3.e5)  # [arcmin^2]
   # from tau and y, get temperature
   tE90 = y / inTau90 # = kB*Te / (me*c^2) [dimless]
   # convert electron temperature to K
   tE90 *= me * c**2 / kB
   # uncertainty on temperature: assume y and tau to be uncorrelated
   sTE90 = tE90 * np.sqrt( (sY/y)**2 + (sInTau90/inTau90)**2 )


   # save as a table
   data = np.zeros((len(rAp), 11))
   data[:,0] = rAp   # [arcmin]
   data[:,1] = y  # [arcmin^2]
   data[:,2] = sY  # [arcmin^2]
   data[:,3] = inTau90  # [arcmin^2]
   data[:,4] = sInTau90  # [arcmin^2]
   data[:,5] = inTau150  # [arcmin^2]
   data[:,6] = sInTau150  # [arcmin^2]
   data[:,7] = tE90  # [K]
   data[:,8] = sTE90  # [K]
   data[:,9] = tE150  # [K]
   data[:,10] = sTE150  # [K]
   np.savetxt(pathFig+'summary_y_tau_te_150_90_'+catalogKey+'.txt', data)


   # Plot the electron temperature
   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   ax.minorticks_on()
   #
   # Expected virial temperature
   ax.axhline(tVir, color='k', ls='--', label=r'$T_\text{vir}$')
   #
   # data
   ax.errorbar(rAp, tE150, sTE150, fmt='o:', c='blue', label=r'from 150 GHz DR5')
   ax.errorbar(rAp + 0.05, tE90, sTE90, fmt='o:', c='purple', label=r'from 98 GHz DR5')
   #
   # virial radius
   ax.axvline(np.sqrt(thetaVir**2 + 2.4**2/(8.*np.log(2.))), color='gray', alpha=0.5)
   #
   ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
   #ax.set_ylim((0., 3.e7))
   ax.set_ylim((1.e6, 1.e8))
   ax.set_yscale('log', nonposy='clip')
   ax.set_xlabel(r'$R$ [arcmin]')
   ax.set_ylabel(r'$T_e \equiv \left(\frac{m_e c^2}{k_B}\right)\left(\frac{y_\text{CAP}}{\tau_\text{CAP}}\right)$ [K]')
   ax.set_title(catalogTitle + r' electron temperature', x=0.5, y=1.25)
   #
   # make extra abscissa with disk comoving size in Mpc/h
   ax2 = ax.twiny()
   ax2.minorticks_on()
   ticks = ax.get_xticks()
   ax2.set_xticks(ticks)
   newticks = np.array(ticks) * np.pi/(180.*60.)*u.bg.comoving_distance(zMean)  # disk radius in Mpc/h
   newticks = np.round(newticks, 2)
   ax2.set_xticklabels(newticks)
   ax2.set_xlim(ax.get_xlim())
   ax2.set_xlabel(r'Comoving radius [Mpc/h] at $z=$'+str(round(zMean,2)), fontsize=20)
   ax2.xaxis.set_label_coords(0.5, 1.15)
   #
   path = pathFig+"electron_temperature_150_90_"+catalogKey+".pdf"
   fig.savefig(path, bbox_inches='tight')
   #plt.show()
   fig.clf()
'''


###################################################################################
###################################################################################
###################################################################################
###################################################################################
# Plot the 2d stacked cutouts


###################################################################################
# map object for plotting

# Reproduce the cutoutGeometry in thumbstack
rApMaxArcmin = 6. # [arcmin]
resCutoutArcmin = 0.25  # [arcmin]
# choose postage stamp size to fit the largest ring
dArcmin = np.ceil(2. * rApMaxArcmin * np.sqrt(2.))
#
nx = np.floor((dArcmin / resCutoutArcmin - 1.) / 2.) + 1.
nX = np.int(2 * nx + 1)
dxDeg = (2. * nx + 1.) * resCutoutArcmin / 60.
# create a map object for plotting
baseMap = FlatMap(nX=nX, nY=nX, sizeX=dxDeg*np.pi/180., sizeY=dxDeg*np.pi/180.)


def plot2dCutout(baseMap, data, title=None, cbTitle=None, cmap='viridis', path=None):
   vMin = np.min(data.flatten())
   vMax = np.max(data.flatten())

   fig=plt.figure(0)
   ax=fig.add_subplot(111)
   #
   # pcolor wants x and y to be edges of cell,
   # ie one more element, and offset by half a cell
   x = baseMap.dX * (np.arange(baseMap.nX+1) - 0.5)
   x -= x[baseMap.nX//2]
   y = baseMap.dY * (np.arange(baseMap.nY+1) - 0.5)
   y -= y[baseMap.nX//2]
   x,y = np.meshgrid(x, y, indexing='ij')
   #
   cp=ax.pcolormesh(x*180.*60./np.pi, y*180.*60./np.pi, data, linewidth=0, rasterized=True)
   #cp=ax.pcolormesh(x*180./np.pi, y*180./np.pi, data, linewidth=0, rasterized=True)
   #
   # get limits for later
   lim = ax.get_xlim()
   #
   # choose color map: jet, summer, winter, Reds, gist_gray, YlOrRd, bwr, seismic
   cp.set_cmap(cmap)
   #cp.set_clim(0.,255.)
   #cp.set_clim(-3.*sigma, 3.*sigma)
   cp.set_clim(vMin, vMax)
   cb=fig.colorbar(cp)
   #
   plt.axis('scaled')
   ax.set_xlim(np.min(x)*180.*60./np.pi, np.max(x)*180.*60./np.pi)
   ax.set_ylim(np.min(y)*180.*60./np.pi, np.max(y)*180.*60./np.pi)
   ax.set_xlabel('$x$ [arcmin]')
   ax.set_ylabel('$y$ [arcmin]')
   if title is not None:
      ax.set_title(title, x=0.5, y=1.25)
   if cbTitle is not None:
      #cb.set_label(cbTitle)
      cb.ax.set_title(cbTitle)
   #
   # make extra abscissa with disk comoving size in Mpc/h
   ax2 = ax.twiny()
   ax2.minorticks_on()
   ticks = ax.get_xticks()
   ax2.set_xticks(ticks)
   newticks = np.array(ticks) * np.pi/(180.*60.)*u.bg.comoving_distance(zMean)  # disk radius in Mpc/h
   newticks = np.round(newticks, 2)
   ax2.set_xticklabels(newticks)
   ax2.set_xlim(ax.get_xlim())
   ax2.set_xlabel(r'Comoving radius [Mpc/h] at $z=$'+str(round(zMean,2)), fontsize=20)
   ax2.xaxis.set_label_coords(0.5, 1.15)
   #
   if path is not None:
      print "saving plot to "+path
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
   else:
      plt.show()


###################################################################################

'''
for catalogKey in ['cmass_kendrick', 'lowz_kendrick']:
   catalog = catalogs[catalogKey]

   # true velocity RMS for CMASS K
   zMean = catalog.Z.mean()   # Msun
   vRms = u.v1dRms(0., zMean, W3d_sth) # [km/s]
   tauToKsz = 2.726e6 * (vRms/3.e5)
   # virial radius, before convolving with the beams
   mVirMean = catalog.Mvir.mean()   # [Msun]. Will need to be converted to [Msun/h] below 
   rVir = u.frvir(mVirMean * u.bg.h, zMean)   # [Mpc/h]
   thetaVir = rVir / u.bg.comoving_distance(zMean) * (180.*60./np.pi) # [arcmin]

   print catalogKey
   print "mean z =", zMean
   print "mean chi =", u.bg.comoving_distance(zMean), "Mpc/h"
   print "mean mass =", mVirMean, "Msun"
   print "virial radius =", rVir, "Mpc/h"
   print "virial angle =", thetaVir, "arcmin"
   print "RMS 1d velocity =", vRms, "km/s"

   if catalogKey=='cmass_kendrick':
      catalogTitle = 'CMASS'
   elif catalogKey=='lowz_kendrick':
      catalogTitle = 'LOWZ'

   rAp = rKsz[catalogKey+'_pactf150daynight20200228maskgal60r2']


   ###################################################################################
   # tSZ + dust, tSZ

   # read the stacked 2d cutouts
   path = "./output/thumbstack/"+catalogKey+"_pactf150daynight20200228maskgal60r2"+"/stackedmap_diskring_tsz_varweight.txt"
   mapTsz150 = np.genfromtxt(path)
   #
   path = "./output/thumbstack/"+catalogKey+"_pactf90daynight20200228maskgal60r2"+"/stackedmap_diskring_tsz_varweight.txt"
   mapTsz90 = np.genfromtxt(path)
   #
   path = "./output/thumbstack/"+catalogKey+"_tilecpacty"+"/stackedmap_diskring_tsz_varweight.txt"
   mapTszTilecy = np.genfromtxt(path) * yTomuK150
   #
   path = "./output/thumbstack/"+catalogKey+"_tilecpactynocib"+"/stackedmap_diskring_tsz_varweight.txt"
   mapTszTilecynocib = np.genfromtxt(path) * yTomuK150


   plot2dCutout(baseMap, mapTsz150, title=catalogTitle+r' tSZ+dust f150', cbTitle=r'$\mu$K', cmap='afmhot', path=pathFig+"stacked_map_"+catalogKey+"_tsz_150.pdf")
   plot2dCutout(baseMap, mapTsz90, title=catalogTitle+r' tSZ+dust f90', cbTitle=r'$\mu$K', cmap='afmhot', path=pathFig+"stacked_map_"+catalogKey+"_tsz_90.pdf")
   plot2dCutout(baseMap, mapTszTilecy, title=catalogTitle+r' tSZ ILC', cbTitle=r'$\mu$K', cmap='afmhot', path=pathFig+"stacked_map_"+catalogKey+"_tsz_tilecy.pdf")
   plot2dCutout(baseMap, mapTszTilecynocib, title=catalogTitle+r' tSZ no dust ILC', cbTitle=r'$\mu$K', cmap='afmhot', path=pathFig+"stacked_map_"+catalogKey+"_tsz_tilecynocib.pdf")


   ###################################################################################
   # kSZ

   # read the stacked 2d cutouts
   path = "./output/thumbstack/"+catalogKey+"_pactf150daynight20200228maskgal60r2"+"/stackedmap_diskring_ksz_varweight.txt"
   mapKsz150 = np.genfromtxt(path)
   #
   path = "./output/thumbstack/"+catalogKey+"_pactf90daynight20200228maskgal60r2"+"/stackedmap_diskring_ksz_varweight.txt"
   mapKsz90 = np.genfromtxt(path)


   plot2dCutout(baseMap, mapKsz150, title=catalogTitle+r' kSZ f150', cbTitle=r'$\mu$K', cmap='viridis', path=pathFig+"stacked_map_"+catalogKey+"_ksz_150.pdf")
   plot2dCutout(baseMap, mapKsz90, title=catalogTitle+r' kSZ f90', cbTitle=r'$\mu$K', cmap='viridis', path=pathFig+"stacked_map_"+catalogKey+"_ksz_90.pdf")
'''






