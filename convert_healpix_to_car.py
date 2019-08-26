from pixell import enmap, utils, powspec, enplot, reproject


# path to true hit count map and mask: Planck + ACT 150GHz day and night
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/"
pathHit = pathIn + "act_planck_f150_prelim_div_mono.fits"
# read maps in common for all mocks
pact150Hit = enmap.read_map(pathHit)

pathMock = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/mock_maps_CMASS_DR12_mariana/"


# Mock 1: counts
pathMap = pathMock + "CMASS_COUNTS_SMOOTH1.5_8192.fits"
hMap = hp.read_map(pathMap)
carMap = reproject.enmap_from_healpix(hMap, pact150Hit.shape, pact150Hit.wcs, rot=None)
enmap.write_map(pathMock + "CMASS_COUNTS_SMOOTH1.5_8192_car.fits", carMap)


# Mock 2: vel
pathMap = pathMock + "CMASS_VEL_SMOOTH1.5_8192.fits"
hMap = hp.read_map(pathMap)
carMap = reproject.enmap_from_healpix(hMap, pact150Hit.shape, pact150Hit.wcs, rot=None)
enmap.write_map(pathMock + "CMASS_VEL_SMOOTH1.5_8192_car.fits", carMap)
