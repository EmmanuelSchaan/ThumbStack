import numpy as np, time
import matplotlib.pyplot as plt
from pixell import enmap, utils, powspec, enplot, reproject
import healpy as hp

#########################################################################

pathIn = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2019_03_11/"

print "read CAR hit map to get the desired pixellation properties"
pathHit = pathIn + "act_planck_f150_prelim_div_mono.fits"
hitMap = enmap.read_map(pathHit)

print "Generate GRF map with enmap"
cov = np.ones(100000)
grfMap = enmap.rand_map(hitMap.shape, hitMap.wcs, cov[None,None])


print "Convert enmap to healpy map"
hpGrfMap = enmap.to_healpix(grfMap)
print "measure the unbinned power spectrum"
power = hp.anafast(hpGrfMap)
power = np.nan_to_num(power)
power /= grfMap.area()/(4.*np.pi)

print "plot power spectrum"
plt.loglog(power, label=r'mock')
plt.axhline(1., label=r'theory')
#
plt.legend(loc=2)
plt.xlabel(r'$\ell$')
plt.ylabel(r'$C_\ell / (2\pi)$')

plt.show()

