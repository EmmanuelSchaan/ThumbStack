from pixell import enmap, utils, powspec, enplot

# Planck + ACT 150GHz day and night
pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/"
pathHit = pathIn + "f150_daynight_all_div_mono.fits"

# load a healpy map
bigmap = enmap.read_map(pathHit)

# setup the interpolation algorithm,
# done once for all, to speed up subsequent calls
bigmap = utils.interpol_prefilter(bigmap, inplace=True)

# Get map value at specific location on the sky
ra = 22. # deg
dec = 0. # deg
sourcecoord = np.array([ra, dec])*utils.degree
bigmap.at(sourcecoord, prefilter=False, mask_nan=False)

