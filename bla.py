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
      #"lowz_kendrick": Catalog(u, massConversion, name="lowz_kendrick", nameLong="LOWZ K", save=False),
      #"boss_kendrick": Catalog(u, massConversion, name="boss_kendrick", nameLong="BOSS K", save=False)
      #"cmass_mk_diff": Catalog(u, massConversion, name="cmass_mk_diff", nameLong="CMASS M-K", save=False)
      }

tStop = time()
print("took "+str(round((tStop-tStart)/60., 2))+" min")


def intersectCatalog(cat, newCat, save=False, vDiff=False, nProc=3):
   '''Take the intersection of the two catalogs.
   Keep the galaxy properties from the first catalog (self).
   If vDiff is True, use the difference of the two velocities, 
   for a null test.
   '''

   # find the intersection
   hasMatch = np.zeros(cat.nObj, dtype=bool)

   def matchObj(iObj):
      if iObj%10000==0:
         print "matching object", iObj
      ra = cat.RA[iObj]
      dec = cat.DEC[iObj]
      z = cat.Z[iObj]

      diff = (newCat.RA - ra)**2 / (1.e-3)**2   # accuracy of Mariana's RA (3.6 arcsec)
      diff += (newCat.DEC - dec)**2 / (1.e-4)**2   # accuracy of Mariana's DEC (0.36 arcsec)
      diff += (newCat.Z - z) **2 / (1.e-4)**2   # accuracy of Mariana's redshifts 
      diff = np.sqrt(diff)

      minDiff = np.min(diff)
      #print "min diff", minDiff
      if (minDiff<1.):
         IMatch = np.where(diff==minDiff)[0]
         if len(IMatch) > 1:
            print "Problem: got", len(IMatch), "matches"
         hasMatch[iObj] = True
         iMatch = IMatch[0]
         #print iObj, minDiff
      else:
         iMatch = -1
      return iMatch



   with sharedmem.MapReduce(np=nProc) as pool:
      #IMatch = np.array(pool.map(matchObj, range(cat.nObj)))
      IMatch = np.array(pool.map(matchObj, range(500)))


   I0Match = np.where(IMatch<>-1.)[0]
   print "First catalog has", cat.nObj, "objects"
   print "Second catalog has", newCat.nObj, "objects"
   print "Intersection has", len(I0Match), "objects"


   cat.nObj = len(I0Match)
   #
   # sky coordinates and redshift
   cat.RA = cat.RA[I0Match]
   cat.DEC = cat.DEC[I0Match]
   cat.Z = cat.Z[I0Match]
   #
   # observed cartesian coordinates
   cat.coordX = cat.coordX[I0Match]
   cat.coordY = cat.coordY[I0Match]
   cat.coordZ = cat.coordZ[I0Match]
   #
   # displacement from difference,
   # not including the Kaiser displacement,
   # from differences of the observed and reconstructed fields
   cat.dX = cat.dX[I0Match]
   cat.dY = cat.dY[I0Match]
   cat.dZ = cat.dZ[I0Match]
   #
   # Kaiser-only displacement
   # originally from differences of the observed and reconstructed fields
   cat.dXKaiser = cat.dXKaiser[I0Match]
   cat.dYKaiser = cat.dYKaiser[I0Match]
   cat.dZKaiser = cat.dZKaiser[I0Match]
   #
   # velocity in cartesian coordinates
   cat.vX = cat.vX[I0Match]
   cat.vY = cat.vY[I0Match]
   cat.vZ = cat.vZ[I0Match]
   #
   # velocity in spherical coordinates,
   # from catalog of spherical displacements
   cat.vR = cat.vR[I0Match] - (vDiff==True) * newCat.vR[IMatch[I0Match]]
   cat.vTheta = cat.vTheta[I0Match] - (vDiff==True) * newCat.vTheta[IMatch[I0Match]]
   cat.vPhi = cat.vPhi[I0Match] - (vDiff==True) * newCat.vPhi[IMatch[I0Match]]
   #
   # Stellar masses
   cat.Mstellar = cat.Mstellar[I0Match]
   #
   # Halo mass
   cat.hasM = cat.hasM[I0Match]
   cat.Mvir = cat.Mvir[I0Match]
   #
   # Integrated optical depth [dimless]: int d^2theta n_e^2d sigma_T = (total nb of electrons) * sigma_T / (a chi)^2
   cat.integratedTau = cat.integratedTau[I0Match]
   #
   # Integrated kSZ signal [muK * sr]: int d^2theta n_e sigma_T (-v/c) Tcmb
   cat.integratedKSZ = cat.integratedKSZ[I0Match]
   #
   # Integrated Y signal [sr]: int d^2theta n_e sigma_T (kB Te / me c^2)
   # needs to be multiplied by Tcmb * f(nu) to get muK
   cat.integratedY = cat.integratedY[I0Match]


   # Write the full catalog to the output path, if needed
   if save:
      cat.writeCatalog()



cmassMKDiff = catalogs['cmass_mariana'].copy(name="cmass_mk_diff", nameLong="CMASS M-K")
#cmassMKDiff.intersectCatalog(cmassKendrick, vDiff=True, save=True, nProc=32)
intersectCatalog(cmassMKDiff, catalogs['cmass_kendrick'], vDiff=True, save=True, nProc=32)
cmassMKDiff = Catalog(u, massConversion, name="cmass_mk_diff", nameLong="CMASS M-K", save=False)




