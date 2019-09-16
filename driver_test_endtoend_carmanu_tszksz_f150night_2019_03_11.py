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
#plt.switch_backend('Agg')

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

# CMASS
#cmassSMariana = Catalog(u, massConversion, name="cmass_s_mariana", nameLong="CMASS S M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt", save=False)
#cmassNMariana = Catalog(u, massConversion, name="cmass_n_mariana", nameLong="CMASS N M", pathInCatalog="../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt", save=False)
# combined catalog
#cmassMariana = cmassSMariana.copy(name="cmass_mariana", nameLong="CMASS M")
#cmassMariana.addCatalog(cmassNMariana, save=True)
cmassMariana = Catalog(u, massConversion, name="cmass_mariana", nameLong="CMASS M", save=False)


# Shuffle velocities to kill the 2-halo term
#cmassMarianaVShuffle = cmassMariana.copy(name="cmass_mariana_vshuffle", nameLong="CMASS M Vshuffle")
#np.random.shuffle(cmassMarianaVShuffle.vR)
#cmassMarianaVShuffle.writeCatalog()
cmassMarianaVShuffle = Catalog(u, massConversion, name="cmass_mariana_vshuffle", nameLong="CMASS M Vshuffle", save=False)



###################################################################################

# keep only the first few objects, to speed things up
#I = range(10000)
#cmassMariana = cmassMariana.extractCatalog(I, name="mini_cmass_mariana", nameLong="mini CMASS M")
#cmassMariana = Catalog(u, massConversion, name="mini_cmass_mariana", nameLong="mini CMASS M", save=False) 

# Shuffle velocities to kill the 2-halo term
#cmassMarianaVShuffle = cmassMariana.copy(name="mini_cmass_mariana_vshuffle", nameLong="mini CMASS M Vshuffle")
#np.random.shuffle(cmassMarianaVShuffle.vR)
#cmassMarianaVShuffle.writeCatalog()
#cmassMarianaVShuffle = Catalog(u, massConversion, name="mini_cmass_mariana_vshuffle", nameLong="mini CMASS M Vshuffle", save=False)


###################################################################################
###################################################################################
# Read CMB mask and hit count

# path to true hit count map and mask: Planck + ACT
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/"
pathHit = pathIn + "act_planck_f150_prelim_div_mono.fits"
pathMask = pathIn + "f150_mask_foot_planck_ps_car.fits"
pathPower = pathIn + "f150_power_T_masked.txt"

# read maps in common for all mocks
pactMask = enmap.read_map(pathMask)
pactHit = enmap.read_map(pathHit)

# theory power spectrum
cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)


###################################################################################
###################################################################################
# Generate mock maps

# Point sources and Gaussian profiles with sigma=1.5'
cmassMariana.generateMockMaps(pactHit, sigma=1.5)

# Same for the catalog with shuffled velocities
#cmassMarianaVShuffle.generateMockMaps(pactHit, sigma=1.5)


###################################################################################
###################################################################################



###################################################################################
###################################################################################

import thumbstack
reload(thumbstack)
from thumbstack import *



pathMap = cmassMariana.pathOut + "mock_count_dirac_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_count_dirac_carmanu"
tsCountDirac = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=False, nProc=nProc)

pathMap = cmassMariana.pathOut + "mock_count_gauss_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_count_gauss_carmanu"
tsCountGauss = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=False, nProc=nProc)

pathMap = cmassMariana.pathOut + "mock_vel_dirac_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_vel_dirac_carmanu"
tsVelDirac = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=False, nProc=nProc)

pathMap = cmassMariana.pathOut + "mock_vel_gauss_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMariana.name + "_pactf150night20190311_test_endtoend_vel_gauss_carmanu"
tsVelGauss = ThumbStack(u, cmassMariana, pactMap, pactMask, pactHit, name=name, nameLong=None, save=False, nProc=nProc)

# Same on mocks with shuffled velocities
pathMap = cmassMarianaVShuffle.pathOut + "mock_vel_dirac_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMarianaVShuffle.name + "_pactf150night20190311_test_endtoend_vel_dirac_carmanu"
tsVelDiracVShuffle = ThumbStack(u, cmassMarianaVShuffle, pactMap, pactMask, pactHit, name=name, nameLong=None, save=False, nProc=nProc)

pathMap = cmassMarianaVShuffle.pathOut + "mock_vel_gauss_car.fits"
pactMap = enmap.read_map(pathMap)
name = cmassMarianaVShuffle.name + "_pactf150night20190311_test_endtoend_vel_gauss_carmanu"
tsVelGaussVShuffle = ThumbStack(u, cmassMarianaVShuffle, pactMap, pactMask, pactHit, name=name, nameLong=None, save=False, nProc=nProc)




###################################################################################



# Gaussian with sigma = 1.5'
profile = tsCountDirac.ftheoryGaussianProfile(1.5)

fig=plt.figure(0)
ax=fig.add_subplot(111)
#
factor =  (180.*60./np.pi)**2
ax.errorbar(tsCountDirac.RApArcmin, factor*tsCountDirac.stackedProfile['tsz_uniformweight'], factor*tsCountDirac.sStackedProfile['tsz_uniformweight'], fmt='--', c='r', label=r'count Dirac')
ax.errorbar(tsCountGauss.RApArcmin, factor*tsCountGauss.stackedProfile['tsz_uniformweight'], factor*tsCountGauss.sStackedProfile['tsz_uniformweight'], fmt='-', c='r', label=r'count Gauss')
ax.errorbar(tsVelDirac.RApArcmin, factor*tsVelDirac.stackedProfile['ksz_uniformweight'], factor*tsVelDirac.sStackedProfile['ksz_uniformweight'], fmt='--', c='b', label=r'vel Dirac')
ax.errorbar(tsVelGauss.RApArcmin, factor*tsVelGauss.stackedProfile['ksz_uniformweight'], factor*tsVelGauss.sStackedProfile['ksz_uniformweight'], fmt='-', c='b', label=r'vel Gauss')
ax.errorbar(tsVelDiracVShuffle.RApArcmin, factor*tsVelDiracVShuffle.stackedProfile['ksz_uniformweight'], factor*tsVelDiracVShuffle.sStackedProfile['ksz_uniformweight'], fmt='--', c='g', label=r'vel Dirac v-shuffle')
ax.errorbar(tsVelGaussVShuffle.RApArcmin, factor*tsVelGaussVShuffle.stackedProfile['ksz_uniformweight'], factor*tsVelGaussVShuffle.sStackedProfile['ksz_uniformweight'], fmt='-', c='g', label=r'vel Gauss v-shuffle')
#
ax.plot(tsCountDirac.RApArcmin, profile, 'k-', label=r'expected')
#
ax.legend(loc=4, fontsize='x-small', labelspacing=0.1)
#
fig.savefig(tsCountDirac.pathFig+"/test_mean_stacked_temperature.pdf")
#fig.clf()

plt.show()

# asymptote is 1.92858407e-08

















###################################################################################


         ra = self.RA[iObj]
         dec = self.DEC[iObj]
         # coordinates in [rad]
         sourcecoord = np.array([dec, ra]) * np.pi/180.
         
         # find pixel indices (float) corresponding to ra, dec
         iY, iX = enmap.sky2pix(countDirac.shape, countDirac.wcs, sourcecoord, safe=True, corner=False)

         # nearest pixel
         jY = np.int(round(iY))
         jX = np.int(round(iX))

         # Check that the object is within the map boundaries
         if jX>=0 and jX<countDirac.shape[1] and jY>=0 and jY<countDirac.shape[0]:
             # fill the pixel
             countDirac[jY, jX] = 1.
             velDirac[jY, jX] = - self.vR[iObj] / 3.e5   # v_r/c  [dimless]

             # check that I filled the right pixel
             if countDirac.at(sourcecoord, prefilter=False, mask_nan=False, order=0)<>1:
                print "Filled the wrong pixel for  object", iObj, ", ra, dec=", ra, dec
                print iX, jX, countDirac.shape[1]
                print iY, jY, countDirac.shape[0]

             # normalize to integrate to 1 over angles in [muK*arcmin^2]
             countDirac[jY, jX] /= pixSizeMap[jY, jX] * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
             velDirac[jY, jX] /= pixSizeMap[jY, jX] * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
             


###################################################################################

   def extractStamp(self, ra, dec, dxDeg=0.3, dyDeg=0.3, resArcmin=0.25, proj='cea', test=False):
      """Extracts a small CEA or CAR map around the given position, with the given angular size and resolution.
      ra, dec in degrees.
      Does it for the map, the mask and the hit count.
      Makes sure that dxDeg = (2*n+1) * resArcmin / 60.,
      such that the object is exactly in the middle of the central pixel.
      """

      # replace dxDeg and dyDeg by the closest (2*n+1) * resArcmin / 60.
      # this ensures that the object is exactly in the middle of the central pixel
      nx = np.floor((dxDeg * 60. / resArcmin - 1.) / 2.) + 1.
      dxDeg = (2. * nx + 1.) * resArcmin / 60.
      ny = np.floor((dyDeg * 60. / resArcmin - 1.) / 2.) + 1.
      dyDeg = (2. * ny + 1.) * resArcmin / 60.

      # define geometry of small square maps to be extracted
      # here dxDeg * dyDeg, with 0.25arcmin pixel
      # the center of the small square map will be specified later
      # car: not equal area, but equally-spaced coordinates
      # cea: equal area pixels, but coordinates are not equally spaced
      # previous kSZ paper: went out to 5arcmin
      # Planck paper CV: went out to 18 arcmin
      shape, wcs = enmap.geometry(np.array([[-0.5*dxDeg,-0.5*dyDeg],[0.5*dxDeg,0.5*dyDeg]])*utils.degree, res=resArcmin*utils.arcmin, proj=proj)
      stampMap = enmap.zeros(shape, wcs)
      stampMask = enmap.zeros(shape, wcs)
      stampHit = enmap.zeros(shape, wcs)

      # coordinates of the square map (between -1 and 1 deg)
      # output map position [{dec,ra},ny,nx]
      opos = stampMap.posmap()

      # coordinate of the center of the square map we want to extract
      # ra, dec in this order
      sourcecoord = np.array([ra, dec])*utils.degree  # convert from degrees to radians

      # corresponding true coordinates on the big healpy map
      ipos = rotfuncs.recenter(opos[::-1], [0,0,sourcecoord[0],sourcecoord[1]])[::-1]

      # extract the small square map by interpolating the big map
      # these are now numpy arrays: the wcs info is gone

#      # Here, I use nearest neighbor interpolation (order=0)
#      stampMap[:,:] = self.cmbMap.at(ipos, prefilter=False, mask_nan=False, order=0)
#      stampMask[:,:] = self.cmbMask.at(ipos, prefilter=False, mask_nan=False, order=0)
#      stampHit[:,:] = self.cmbHit.at(ipos, prefilter=False, mask_nan=False, order=0)
      
      # Here, I use bilinear interpolation
      stampMap[:,:] = self.cmbMap.at(ipos, prefilter=True, mask_nan=False, order=1)
      stampMask[:,:] = self.cmbMask.at(ipos, prefilter=True, mask_nan=False, order=1)
      stampHit[:,:] = self.cmbHit.at(ipos, prefilter=True, mask_nan=False, order=1)

#      # Here, I use bicubic spline interpolation
#      stampMap[:,:] = self.cmbMap.at(ipos, prefilter=True, mask_nan=False, order=3)
#      stampMask[:,:] = self.cmbMask.at(ipos, prefilter=True, mask_nan=False, order=3)
#      stampHit[:,:] = self.cmbHit.at(ipos, prefilter=True, mask_nan=False, order=3)


      # re-threshold the mask map, to keep 0 and 1 only
      stampMask[:,:] = 1.*(stampMask[:,:]>0.5)



###################################################################################




