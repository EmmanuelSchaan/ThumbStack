# Last full run for VelDiraacVShuffle and VelGaussVShuffle: 2020/01/01

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

import flat_map
reload(flat_map)
from flat_map import *



# running on cori
# 68 cores per knl node, 32 cores per haswell node
#salloc -N 1 --qos=interactive -C haswell -t 04:00:00 -L SCRATCH
#plt.switch_backend('Agg')

##################################################################################

pathFig = "./figures/cmb_map/visu/" 
if not os.path.exists(pathFig):
    os.makedirs(pathFig)
pathOut = "./output/cmb_map/visu/" 
if not os.path.exists(pathOut):
    os.makedirs(pathOut)



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
# Generate mock maps
'''
# Hit count
pathHit = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits"
pactHit = enmap.read_map(pathHit)[0]

# Point sources and Gaussian profiles with sigma=1.5'
cmassMariana.generateMockMaps(pactHit, sigma=1.5)

# Same for the catalog with shuffled velocities
cmassMarianaVShuffle.generateMockMaps(pactHit, sigma=1.5)
'''

###################################################################################
###################################################################################
# Generate empty square map

# Generate an empty square map with RA in [0., 5.] and DEC in [-5., 0.] 
box = np.array([[-5., 5.], [0., 0.]]) * utils.degree
#box = np.array([[-5., 0.], [0., 5.]]) * utils.degree
resArcmin = 0.1 #1. #0.5  # 0.1   # map pixel size [arcmin]
shape,wcs = enmap.geometry(pos=box, res=resArcmin * utils.arcmin, proj='car')
# create a mask that keeps the whole area
#boxMask = enmap.zeros((1,) + shape, wcs=wcs)
boxMask = enmap.ones(shape, wcs=wcs)

# create flat map object for plotting
baseMap = FlatMap(nX=shape[0], nY=shape[1], sizeX=5.*np.pi/180., sizeY=5.*np.pi/180., name="test")



###################################################################################
###################################################################################
# Prepare cutouts for visualization

'''
###################################################################################
# create GRF map with CMB only, beamed

cmb1_4 = StageIVCMB(beam=1.4, noise=30., lMin=1., lMaxT=1.e5, lMaxP=1.e5, atm=False)
flensedTT = lambda l: cmb1_4.flensedTT(l) * cmb1_4.fbeam(l)**2
fdetectorNoise = lambda l: cmb1_4.fdetectorNoise(l) * cmb1_4.fbeam(l)**2
ftotal = lambda l: cmb1_4.ftotal(l) * cmb1_4.fbeam(l)**2

cutGRFCMBFourier = baseMap.genGRF(flensedTT)
cutGRFCMB = baseMap.inverseFourier(cutGRFCMBFourier)
#
np.savetxt(pathOut + "grf_cmb.txt", cutGRFCMB)


###################################################################################
# Read GRF map with full noise, and extract cutout

iMock = 0
pathMockGRF = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/code/thumbstack/output/cmb_map/mocks_grf_planck_act_coadd_2019_03_11/"
pathMockGRF += "mock_"+str(iMock)+"_grf_f150_daynight.fits"
pactMap = enmap.read_map(pathMockGRF)[0]
# extract a cutout
cutGRFFull = enmap.project(pactMap, boxMask.shape, boxMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
#
np.savetxt(pathOut + "grf_full.txt", cutGRFFull)


###################################################################################
# Read mock maps and extract cutout

pathMap = cmassMariana.pathOut + "mock_count_dirac_car.fits"
pactMap = enmap.read_map(pathMap)
# extract a cutout
cutCountDirac = enmap.project(pactMap, boxMask.shape, boxMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
#
np.savetxt(pathOut + "count_dirac.txt", cutCountDirac)


pathMap = cmassMariana.pathOut + "mock_count_gauss_car.fits"
pactMap = enmap.read_map(pathMap)
# extract a cutout
cutCountGauss = enmap.project(pactMap, boxMask.shape, boxMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
#
np.savetxt(pathOut + "count_gauss.txt", cutCountGauss)


pathMap = cmassMariana.pathOut + "mock_vel_dirac_car.fits"
pactMap = enmap.read_map(pathMap)
# extract a cutout
cutVelDirac = enmap.project(pactMap, boxMask.shape, boxMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
#
np.savetxt(pathOut + "vel_dirac.txt", cutVelDirac)


pathMap = cmassMariana.pathOut + "mock_vel_gauss_car.fits"
pactMap = enmap.read_map(pathMap)
# extract a cutout
cutVelGauss = enmap.project(pactMap, boxMask.shape, boxMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
#
np.savetxt(pathOut + "vel_gauss.txt", cutVelGauss)


###################################################################################
# Same on mocks with shuffled velocities

pathMap = cmassMarianaVShuffle.pathOut + "mock_vel_dirac_car.fits"
pactMap = enmap.read_map(pathMap)
# extract a cutout
cutVelDiracShuffled = enmap.project(pactMap, boxMask.shape, boxMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
#
np.savetxt(pathOut + "vel_dirac_vshuffled.txt", cutVelDiracShuffled)


pathMap = cmassMarianaVShuffle.pathOut + "mock_vel_gauss_car.fits"
pactMap = enmap.read_map(pathMap)
# extract a cutout
cutVelGaussShuffled = enmap.project(pactMap, boxMask.shape, boxMask.wcs, order=1, mode='constant', cval=0.0, force=False, prefilter=True, mask_nan=True, safe=True)
#
np.savetxt(pathOut + "vel_gauss_vshuffled.txt", cutVelGaussShuffled)
'''

###################################################################################
###################################################################################
# Read the cutout maps

cutGRFCMB = np.genfromtxt(pathOut + "grf_cmb.txt")
cutGRFFull = np.genfromtxt(pathOut + "grf_full.txt")
cutCountDirac = np.genfromtxt(pathOut + "count_dirac.txt")
cutCountGauss = np.genfromtxt(pathOut + "count_gauss.txt")
cutVelDirac = np.genfromtxt(pathOut + "vel_dirac.txt")
cutVelGauss = np.genfromtxt(pathOut + "vel_gauss.txt")
cutVelDiracShuffled = np.genfromtxt(pathOut + "vel_dirac_vshuffled.txt")
cutVelGaussShuffled = np.genfromtxt(pathOut + "vel_gauss_vshuffled.txt")



###################################################################################
# Visualize small cutouts of the maps



path = pathFig + "grf_cmb.pdf"
baseMap.plot(cutGRFCMB, save=True, path=path, cmap='seismic', vMin=-300., vMax=300., figSize=None, title=r'Lensed CMB', cbTitle=r'$\mu$K')
#
path = pathFig + "grf_full.pdf"
baseMap.plot(cutGRFFull, save=True, path=path, cmap='seismic', vMin=-300., vMax=300., figSize=None, title=r'Mock ACT', cbTitle=r'$\mu$K')
#
path = pathFig + "count_dirac.pdf"
baseMap.plot(cutCountDirac, save=True, path=path, cmap='Reds', vMin=0., vMax=2., figSize=None, title=r'tSZ, Pointlike', cbTitle=r'$\mu$K')
#
path = pathFig + "count_gauss.pdf"
baseMap.plot(cutCountGauss, save=True, path=path, cmap='Reds', vMin=0., vMax=0.5, figSize=None, title=r'tSZ, Gaussian', cbTitle=r'$\mu$K')
#
path = pathFig + "vel_dirac.pdf"
baseMap.plot(cutVelDirac, save=True, path=path, cmap='seismic', vMin=-0.2, vMax=0.2, figSize=None, title=r'kSZ, pointlike', cbTitle=r'$\mu$K')
#
path = pathFig + "vel_gauss.pdf"
baseMap.plot(cutVelGauss, save=True, path=path, cmap='seismic', vMin=-0.2, vMax=0.2, figSize=None, title=r'kSZ, Gaussian', cbTitle=r'$\mu$K')
#
path = pathFig + "vel_dirac_shuffled.pdf"
baseMap.plot(cutVelDiracShuffled, save=True, path=path, cmap='seismic', vMin=-0.2, vMax=0.2, figSize=None, title=r'kSZ, Pointlike, v-shuffled', cbTitle=r'$\mu$K')
#
path = pathFig + "vel_gauss_shuffled.pdf"
baseMap.plot(cutVelGaussShuffled, save=True, path=path, cmap='seismic', vMin=-0.2, vMax=0.2, figSize=None, title=r'kSZ, Gaussian, v-shuffled', cbTitle=r'$\mu$K')





