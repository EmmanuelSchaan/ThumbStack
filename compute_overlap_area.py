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

#cmassMariana = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False)
cmassKendrick = Catalog(u, massConversion, name="cmass_kendrick", nameLong="CMASS K", save=False)
#lowzKendrick = Catalog(u, massConversion, name="lowz_kendrick", nameLong="LOWZ K", save=False)

# read overlap flags
# for PACT maps
#path = "./output/thumbstack/cmass_mariana_pactf150daynight20200228maskgal60r2/overlap_flag.txt"
#overlapCmassM = np.genfromtxt(path).astype(bool)
path = "./output/thumbstack/cmass_kendrick_pactf150daynight20200228maskgal60r2/overlap_flag.txt"
overlapCmassK = np.genfromtxt(path).astype(bool)
#path = "./output/thumbstack/lowz_kendrick_pactf150daynight20200228maskgal60r2/overlap_flag.txt"
#overlapLowzK = np.genfromtxt(path).astype(bool)

#print "CMASS M:", np.sum(overlapCmassM), "galaxies overlap out of", len(overlapCmassM), "ie", 1.*np.sum(overlapCmassM)/len(overlapCmassM),"%"
print "CMASS K:", np.sum(overlapCmassK), "galaxies overlap out of", len(overlapCmassK), "ie", 1.*np.sum(overlapCmassK)/len(overlapCmassK),"%"
#print "LOWZ K:", np.sum(overlapLowzK), "galaxies overlap out of", len(overlapLowzK), "ie", 1.*np.sum(overlapLowzK)/len(overlapLowzK),"%"


###################################################################################
print("Read PACT 150 map")

pathMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits"
# read maps
actMap = enmap.read_map(pathMap)[0]


###################################################################################
# Generate the map

import pointsrcs


def generateCountGaussMap(ra, dec, carMap, sigma):
   '''ra, dec [deg]
   sigma [arcmin]
   '''
   # input for the counts
   srcsCount = np.zeros((len(ra), 3))
   srcsCount[:,0] = dec * np.pi/180.   # [rad]
   srcsCount[:,1] = ra * np.pi/180. # [rad]
   srcsCount[:,2] = 1.
   # create the map
   countGauss = pointsrcs.sim_srcs(carMap.shape, carMap.wcs, srcsCount, sigma*np.pi/(180.*60.))
   return countGauss





#def generateCountDiracMap(ra, dec, carMap):
#   '''ra, dec [deg]
#   '''
#   print "- Generate mock maps"
#   tStart = time()
#   # create empty maps
#   countDirac = carMap.copy()
#   countDirac[:,:] = 0.
#   # get map of ra and dec, just to check
#   posmap = countDirac.posmap()
#
#   for iObj in range(self.nObj):
##      for iObj in range(10):
#      if iObj%100000==0:
#         print "    -", iObj
#      # object coordinates [deg]
#      ra = self.RA[iObj]
#      dec = self.DEC[iObj]
#      # coordinates in [rad]
#      sourcecoord = np.array([dec, ra]) * np.pi/180.
#      # find pixel indices (float) corresponding to ra, dec
#      iY, iX = enmap.sky2pix(countDirac.shape, countDirac.wcs, sourcecoord, safe=True, corner=False)
#      if test:
#         print 'ra, dec =', ra, dec, iY, iX
#         print countDirac.shape
#      # Check that the object is within the map boundaries
#      # before rounding the indices
#      if iX>=0 and iX<=(countDirac.shape[1]-1) and iY>=0 and iY<=(countDirac.shape[0]-1):
#          # nearest pixel
#          # watch out for difference round VS np.round!
#          jY = np.int(round(iY))
#          jX = np.int(round(iX))
#          if test:
#            print("Object "+str(iObj)+" overlaps")
#          # fill the pixel
#          countDirac[jY, jX] = 1.
#          # check that I filled the right pixel
#          if countDirac.at(sourcecoord, prefilter=False, mask_nan=False, order=0)<>1:
#             print "Filled the wrong pixel for  object", iObj
#             print "wanted ra, dec=", ra, dec # [deg]
#             print "chosen closest ra, dec=", posmap[::-1, jY, jX] * 180./np.pi  # [deg]
#             print "difference in arcmin=", (posmap[::-1, jY, jX] * 180./np.pi - np.array([ra, dec]))*60.  # residual in [arcmin]
#             print "ra index=", iX, jX, np.int(np.round(iX)), countDirac.shape[1]
#             print "dec index=", iY, jY, np.int(np.round(iY)), countDirac.shape[0]
#
#   return countDirac














ra = cmassKendrick.RA[overlapCmassK]
dec = cmassKendrick.DEC[overlapCmassK]
sigma = 0.1 # [arcmin]


tStart = time()
countGaussMap = generateCountGaussMap(ra, dec, actMap, sigma)
tStop = time()
print("Took "+str((tStop-tStart)/60.)+" min")


footMap = 1. * (countGaussMap<>0.)
pixArea = footMap.pixsize() * (180./np.pi)**2   # [deg^2
footArea = np.sum(footMap) * pixArea
print('map area is '+str(footArea)+' deg2')






#tStart = time()
#countDiracMap = generateCountDiracMap(ra, dec, actMap)
#tStop = time()
#print("Took "+str((tStop-tStart)/60.)+" min")









