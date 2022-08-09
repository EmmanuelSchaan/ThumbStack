from headers import *

##################################################################################
##################################################################################

class ThumbStack(object):

#   def __init__(self, U, Catalog, pathMap="", pathMask="", pathHit="", name="test", nameLong=None, save=False, nProc=1):
   def __init__(self, U, Catalog, cmbMap, cmbMask, cmbHit=None, name="test", nameLong=None, save=False, nProc=1, filterTypes='diskring', doStackedMap=False, doMBins=False, doVShuffle=False, doBootstrap=False, cmbNu=150.e9, cmbUnitLatex=r'$\mu$K'):
      
      self.nProc = nProc
      self.U = U
      self.Catalog = Catalog
      self.name = name
      if nameLong is None:
         self.nameLong = self.name
      else:
         self.nameLong = nameLong
      self.cmbMap = cmbMap
      self.cmbMask = cmbMask
      self.cmbHit = cmbHit
      self.doMBins = doMBins
      self.doVShuffle = doVShuffle
      self.doBootstrap = doBootstrap
      self.cmbNu = cmbNu
      self.cmbUnitLatex = cmbUnitLatex

      # aperture photometry filters to implement
      if filterTypes=='diskring':
         self.filterTypes = np.array(['diskring']) 
      elif filterTypes=='disk':
         self.filterTypes = np.array(['disk'])
      elif filterTypes=='ring':
         self.filterTypes = np.array(['ring'])
      elif filterTypes=='all':
         self.filterTypes = np.array(['diskring', 'disk', 'ring'])

      # estimators (ksz, tsz) and weightings (uniform, hit, var, ...)
      # for stacked profiles, bootstrap cov and v-shuffle cov
      if self.cmbHit is not None:
         self.Est = ['tsz_uniformweight', 'tsz_varweight']   #['tsz_uniformweight', 'tsz_hitweight', 'tsz_varweight', 'ksz_uniformweight', 'ksz_hitweight', 'ksz_varweight', 'ksz_massvarweight']
         self.EstBootstrap = ['tsz_uniformweight', 'tsz_varweight']  #['tsz_varweight', 'ksz_varweight']
         self.EstVShuffle = []   #['ksz_varweight']
         self.EstMBins = ['tsz_uniformweight', 'tsz_varweight']# ['tsz_varweight', 'ksz_varweight']
      else:
         self.Est = ['tsz_uniformweight'] #['tsz_uniformweight', 'ksz_uniformweight', 'ksz_massvarweight']
         self.EstBootstrap = ['tsz_uniformweight'] #['tsz_uniformweight', 'ksz_uniformweight']
         self.EstVShuffle = []   #['ksz_uniformweight']
         self.EstMBins = ['tsz_uniformweight'] #['ksz_uniformweight'] #['tsz_uniformweight', 'ksz_uniformweight']

      # resolution of the cutout maps to be extracted
      self.resCutoutArcmin = 0.25   # [arcmin]
      # projection of the cutout maps
      self.projCutout = 'cea'

      # number of samples for bootstraps, shuffles
      self.nSamples = 10000

      # number of mMax cuts to test,
      # for tSZ contamination to kSZ
      self.nMMax = 20
      # fiducial mass cuts, to avoid eg tSZ contamination
      # from massive clusters
      self.mMin = 0. #1.e6
      self.mMax = np.inf #1.e14 # 1.e17

      
      # Output path
      self.pathOut = "./output/thumbstack/"+self.name
      if not os.path.exists(self.pathOut):
         os.makedirs(self.pathOut)

      # Figures path
      self.pathFig = "./figures/thumbstack/"+self.name
      if not os.path.exists(self.pathFig):
         os.makedirs(self.pathFig)
      # test figures path
      self.pathTestFig = self.pathFig+"/tests"
      if not os.path.exists(self.pathTestFig):
         os.makedirs(self.pathTestFig)

      print "- Thumbstack: "+str(self.name)
      
      self.loadAPRadii()
      self.loadMMaxBins()
      
      if save:
         self.saveOverlapFlag(nProc=self.nProc)
      self.loadOverlapFlag()
      
      if save:
         self.saveFiltering(nProc=self.nProc)
      self.loadFiltering()
      
      self.measureAllVarFromHitCount(plot=save)
      
#      self.measureAllMeanTZBins(plot=save, test=False)


      if save:
         self.saveAllStackedProfiles()
      self.loadAllStackedProfiles()

      if save:
         self.plotAllStackedProfiles()
         self.plotAllCov()
         self.computeAllSnr()

      #if save:
      if True:
         if doStackedMap:
            # save all stacked maps
            #self.saveAllStackedMaps()
            # save only the stacked maps for
            # the best tsz and ksz estimators,
            # and for the diskring weighting
            #self.saveAllStackedMaps(filterTypes=['diskring'], Est=['tsz_varweight', 'ksz_varweight'])
            self.saveAllStackedMaps(filterTypes=None, Est=None)



   ##################################################################################
   ##################################################################################
   
   def loadAPRadii(self):
      # radii to use for AP filter: comoving Mpc/h
      self.nRAp = 9 #30 #9  #4
      
      # Aperture radii in Mpc/h
      #self.rApMinMpch = 1.
      #self.rApMaxMpch = 5
      #self.RApMpch = np.linspace(self.rApMinMpch, self.rApMaxMpch, self.nRAp)
      
      # Aperture radii in arcmin
      self.rApMinArcmin = 1.  #0.1   #1.  # 1.
      self.rApMaxArcmin = 6.  #6.  #6.  # 4.
      self.RApArcmin = np.linspace(self.rApMinArcmin, self.rApMaxArcmin, self.nRAp)


   def cutoutGeometry(self, test=False):
      '''Create enmap for the cutouts to be extracted.
      Returns a null enmap object with the right shape and wcs.
      '''
      
      # choose postage stamp size to fit the largest ring
      dArcmin = np.ceil(2. * self.rApMaxArcmin * np.sqrt(2.))
      #
      nx = np.floor((dArcmin / self.resCutoutArcmin - 1.) / 2.) + 1.
      dxDeg = (2. * nx + 1.) * self.resCutoutArcmin / 60.
      ny = np.floor((dArcmin / self.resCutoutArcmin - 1.) / 2.) + 1.
      dyDeg = (2. * ny + 1.) * self.resCutoutArcmin / 60.

      # define geometry of small square maps to be extracted
      shape, wcs = enmap.geometry(np.array([[-0.5*dxDeg,-0.5*dyDeg],[0.5*dxDeg,0.5*dyDeg]])*utils.degree, res=self.resCutoutArcmin*utils.arcmin, proj=self.projCutout)
      cutoutMap = enmap.zeros(shape, wcs)

      if test:
         print "cutout sides are dx, dy =", dxDeg*60., ",", dyDeg*60. , "arcmin"
         print "cutout pixel dimensions are", shape
         print "hence a cutout resolution of", dxDeg*60./shape[0], ",", dyDeg*60./shape[1], "arcmin per pixel"
         print "(requested", self.resCutoutArcmin, "arcmin per pixel)"

      return cutoutMap


   def loadMMaxBins(self, test=False):
      '''Choose the mMax values to have the same number of galaxies
      added in the sample for each mMax increment.
      '''
#      self.MMax = np.interp(np.linspace(0, self.Catalog.nObj, self.nMMax+1),
#                           np.arange(self.Catalog.nObj),
#                           np.sort(self.Catalog.Mvir))[1:]
      
      self.MMax = np.logspace(np.log10(self.Catalog.Mvir.min()*1.1), np.log10(self.Catalog.Mvir.max()), self.nMMax)

      if test:
         print "Checking the mMax bins:"
         print "number of bins:", self.nMMax, len(self.MMax)
         print "values of mMax:", self.MMax

      

   ##################################################################################
   

   def sky2map(self, ra, dec, map):
      '''Gives the map value at coordinates (ra, dec).
      ra, dec in degrees.
      Uses nearest neighbor, no interpolation.
      Will return 0 if the coordinates requested are outside the map
      '''
      # interpolate the map to the given sky coordinates
      sourcecoord = np.array([dec, ra]) * utils.degree   # convert from degrees to radians
      # use nearest neighbor interpolation
      return map.at(sourcecoord, prefilter=False, mask_nan=False, order=0)
   
   
   
   #def saveOverlapFlag(self, thresh=1.e-5, nProc=1):
   def saveOverlapFlag(self, thresh=0.95, nProc=1):
      '''1 for objects that overlap with the hit map,
      0 for objects that don't.
      '''
      print "Create overlap flag"
      # find if a given object overlaps with the CMB hit map
      def foverlap(iObj):
         '''Returns 1. if the object overlaps with the hit map and 0. otherwise.
         '''
         if iObj%100000==0:
            print "-", iObj
         ra = self.Catalog.RA[iObj]
         dec = self.Catalog.DEC[iObj]
         #hit = self.sky2map(ra, dec, self.cmbHit)
         hit = self.sky2map(ra, dec, self.cmbMask)
         return np.float(hit>thresh)
      
#      overlapFlag = np.array(map(foverlap, range(self.Catalog.nObj)))
      tStart = time()
      with sharedmem.MapReduce(np=nProc) as pool:
         overlapFlag = np.array(pool.map(foverlap, range(self.Catalog.nObj)))
      tStop = time()
      print "took", (tStop-tStart)/60., "min"
      print "Out of", self.Catalog.nObj, "objects,", np.sum(overlapFlag), "overlap, ie a fraction", np.sum(overlapFlag)/self.Catalog.nObj
      np.savetxt(self.pathOut+"/overlap_flag.txt", overlapFlag)
   
   
   def loadOverlapFlag(self):
      self.overlapFlag = np.genfromtxt(self.pathOut+"/overlap_flag.txt")
   
   
   ##################################################################################
   
   def examineCmbMaps(self):
   
      # mask, before re-thresholding
      x = self.cmbMask.copy()
      path = self.pathFig+"/hist_cmbmask_prerethresh.pdf"
      myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'CMB mask value', semilogy=True)

      # rethreshold the mask
      mask = (self.cmbMask>0.5)[0]

      # mask, after re-thresholding
      x = 1. * mask.copy()
      path = self.pathFig+"/hist_cmbmask_postrethresh.pdf"
      myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'CMB mask value', semilogy=True)

      # masked map histogram
      x = self.cmbMap[mask]
      path = self.pathFig+"/hist_cmbmap.pdf"
      myHistogram(x, nBins=71, lim=(-10.*np.std(x), 10.*np.std(x)), path=path, nameLatex=r'CMB map value', semilogy=True, doGauss=True, S2Theory=[110.**2])

      # masked hit count histogram
      x = self.cmbHit[mask]
      path = self.pathFig+"/hist_cmbhit.pdf"
      myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'CMB hit count', semilogy=True)
   
   
   ##################################################################################
   

   def extractStamp(self, ra, dec, test=False):
      """Extracts a small CEA or CAR map around the given position, with the given angular size and resolution.
      ra, dec in degrees.
      Does it for the map, the mask and the hit count.
      """

      stampMap = self.cutoutGeometry()
      stampMask = stampMap.copy()
      stampHit = stampMap.copy()

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
      if self.cmbHit is not None:
         stampHit[:,:] = self.cmbHit.at(ipos, prefilter=True, mask_nan=False, order=1)

#      # Here, I use bicubic spline interpolation
#      stampMap[:,:] = self.cmbMap.at(ipos, prefilter=True, mask_nan=False, order=3)
#      stampMask[:,:] = self.cmbMask.at(ipos, prefilter=True, mask_nan=False, order=3)
#      stampHit[:,:] = self.cmbHit.at(ipos, prefilter=True, mask_nan=False, order=3)


      # re-threshold the mask map, to keep 0 and 1 only
      stampMask[:,:] = 1.*(stampMask[:,:]>0.5)

      if test:
         print "Extracted cutouts around ra=", ra, "dec=", dec
         
         print "Map:"
         print "- min, mean, max =", np.min(stampMap), np.mean(stampMap), np.max(stampMap)
         print "- plot"
         plots=enplot.plot(stampMap, grid=True)
         enplot.write(self.pathTestFig+"/stampmap_ra"+np.str(np.round(ra, 2))+"_dec"+np.str(np.round(dec, 2)), plots)

         print "Mask:"
         print "- min, mean, max =", np.min(stampMask), np.mean(stampMask), np.max(stampMask)
         print "- plot"
         plots=enplot.plot(stampMask, grid=True)
         enplot.write(self.pathTestFig+"/stampmask_ra"+np.str(np.round(ra, 2))+"_dec"+np.str(np.round(dec, 2)), plots)

         print "Hit count:"
         print "- min, mean, max =", np.min(stampHit), np.mean(stampHit), np.max(stampHit)
         print "- plot the hit"
         plots=enplot.plot(stampHit, grid=True)
         enplot.write(self.pathTestFig+"/stamphit_ra"+np.str(np.round(ra, 2))+"_dec"+np.str(np.round(dec, 2)), plots)

      return opos, stampMap, stampMask, stampHit



   ##################################################################################


   def aperturePhotometryFilter(self, opos, stampMap, stampMask, stampHit, r0, r1, filterType='diskring',  test=False):
      """Apply an AP filter (disk minus ring) to a stamp map:
      AP = int d^2theta * Filter * map.
      Unit is [map unit * sr]
      The filter function is dimensionless:
      Filter = 1 in the disk, - (disk area)/(ring area) in the ring, 0 outside.
      Hence:
      int d^2theta * Filter = 0.
      r0 and r1 are the radius of the disk and ring in radians.
      stampMask should have values 0 and 1 only.
      Output:
      filtMap: [map unit * sr]
      filtMask: [mask unit * sr]
      filtHitNoiseStdDev: [1/sqrt(hit unit) * sr], ie [std dev * sr] if [hit map] = inverse var
      diskArea: [sr]
      """
      # coordinates of the square map in radians
      # zero is at the center of the map
      # output map position [{dec,ra},ny,nx]
      dec = opos[0,:,:]
      ra = opos[1,:,:]
      radius = np.sqrt(ra**2 + dec**2)
      # exact angular area of a pixel [sr] (same for all pixels in CEA, not CAR)
      pixArea = ra.area() / len(ra.flatten())
      
      # detect point sources within the filter:
      # gives 0 in the absence of point sources/edges; gives >=1 in the presence of point sources/edges
      filtMask = np.sum((radius<=r1) * (1-stampMask))   # [dimensionless]

      # disk filter [dimensionless]
      inDisk = 1.*(radius<=r0)
      # exact angular area of disk [sr]
      diskArea = np.sum(inDisk) * pixArea
      # ring filter [dimensionless]
      inRing = 1.*(radius>r0)*(radius<=r1)

      if filterType=='diskring':
         # normalize the ring so that the disk-ring filter integrates exactly to zero
         inRing *= np.sum(inDisk) / np.sum(inRing)
         # disk minus ring filter [dimensionless]
         filterW = inDisk - inRing
         if np.isnan(np.sum(filterW)):
            print "filterW sums to nan", r0, r1, np.sum(radius), np.sum(1.*(radius>r0)), np.sum(1.*(radius>r0)*(radius<=r1))
      elif filterType=='disk':
         # disk filter [dimensionless]
         inDisk = 1.*(radius<=r0)
         filterW = inDisk
      elif filterType=='ring':
         filterW = inRing
      elif filterType=='meanring':
         filterW = inRing / np.sum(pixArea * inRing)

      # apply the filter: int_disk d^2theta map -  disk_area / ring_area * int_ring d^2theta map
      filtMap = np.sum(pixArea * filterW * stampMap)   # [map unit * sr]
      # quantify noise std dev in the filter
      if self.cmbHit is not None:
         filtHitNoiseStdDev = np.sqrt(np.sum((pixArea * filterW)**2 / (1.e-16 + stampHit))) # to get the std devs [sr / sqrt(hit unit)]
      else:
         filtHitNoiseStdDev = 0.
      

      #print "filtHitNoiseStdDev = ", filtHitNoiseStdDev
      if np.isnan(filtHitNoiseStdDev):
         print "filtHitNoiseStdDev is nan"

      if test:
         print "AP filter with disk radius =", r0 * (180.*60./np.pi), "arcmin"
         # count nb of pixels where filter is strictly positive
         nbPix = len(np.where(filterW>0.)[0])
         print "- nb of pixels in the cutout: "+str(filterW.shape[0] * filterW.shape[1])
         print "= nb of pixels where filter=0: "+str(len(np.where(filterW==0.)[0]))
         print "+ nb of pixels where filter>0: "+str(len(np.where(filterW>0.)[0]))
         print "+ nb of pixels where filter<0: "+str(len(np.where(filterW<0.)[0]))
         print "- disk area: "+str(diskArea)+" sr, ie "+str(diskArea * (180.*60./np.pi)**2)+"arcmin^2"
         print "  (from r0, expect "+str(np.pi*r0**2)+" sr, ie "+str(np.pi*r0**2 * (180.*60./np.pi)**2)+"arcmin^2)"
         print "- disk-ring filter sums over pixels to "+str(np.sum(filterW))
         print "  (should be 0; compared to "+str(len(filterW.flatten()))+")"
         print "- filter on unit disk: "+str(np.sum(pixArea * filterW * inDisk))
         print "  (should be disk area in sr: "+str(diskArea)+")"
         print "- filter on map: "+str(filtMap)
         print "- filter on mask: "+str(filtMask)
         print "- filter on inverse hit: "+str(filtHitNoiseStdDev)
         print "- plot the filter"
         filterMap = stampMap.copy()
         filterMap[:,:] = filterW.copy()

         vmax = np.max(np.abs(filterW))
         plots=enplot.plot(filterMap,grid=True, min=-vmax, max=vmax)#, color='hotcold')
         enplot.write(self.pathTestFig+"/stampfilter_r0"+floatExpForm(r0)+"_r1"+floatExpForm(r1), plots)

      return filtMap, filtMask, filtHitNoiseStdDev, diskArea




   ##################################################################################

   def analyzeObject(self, iObj, test=False):
      '''Analysis to be done for each object: extract cutout map once,
      then apply all aperture photometry filters requested on it.
      Returns:
      filtMap: [map unit * sr]
      filtMask: [mask unit * sr]
      filtHitNoiseStdDev: [1/sqrt(hit unit) * sr], ie [std dev * sr] if [hit map] = inverse var
      diskArea: [sr]
      '''
      
      if iObj%10000==0:
         print "- analyze object", iObj
      
      
      # create arrays of filter values for the given object
      filtMap = {}
      filtMask = {}
      filtHitNoiseStdDev = {} 
      filtArea = {}
      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]
         # create arrays of filter values for the given object
         filtMap[filterType] = np.zeros(self.nRAp)
         filtMask[filterType] = np.zeros(self.nRAp)
         filtHitNoiseStdDev[filterType] = np.zeros(self.nRAp)
         filtArea[filterType] = np.zeros(self.nRAp)
      
      # only do the analysis if the object overlaps with the CMB map
      if self.overlapFlag[iObj]:
         # Object coordinates
         ra = self.Catalog.RA[iObj]   # in deg
         dec = self.Catalog.DEC[iObj] # in deg
         z = self.Catalog.Z[iObj]
         # choose postage stamp size to fit the largest ring
         dArcmin = np.ceil(2. * self.rApMaxArcmin * np.sqrt(2.))
         dDeg = dArcmin / 60.
         # extract postage stamp around it
         opos, stampMap, stampMask, stampHit = self.extractStamp(ra, dec, test=test)
         
         for iFilterType in range(len(self.filterTypes)):
            filterType = self.filterTypes[iFilterType]

            # loop over the radii for the AP filter
            for iRAp in range(self.nRAp):
               ## disk radius in comoving Mpc/h
               #rApMpch = self.RApMpch[iRAp]
               ## convert to radians at the given redshift
               #r0 = rApMpch / self.U.bg.comoving_transverse_distance(z) # rad

               # Disk radius in rad
               r0 = self.RApArcmin[iRAp] / 60. * np.pi/180.
               # choose an equal area AP filter
               r1 = r0 * np.sqrt(2.)
               
               # perform the filtering
               filtMap[filterType][iRAp], filtMask[filterType][iRAp], filtHitNoiseStdDev[filterType][iRAp], filtArea[filterType][iRAp] = self.aperturePhotometryFilter(opos, stampMap, stampMask, stampHit, r0, r1, filterType=filterType, test=test)

      if test:
         print " plot the measured profile"
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         for iFilterType in range(len(self.filterTypes)):
            filterType = self.filterTypes[iFilterType]
            ax.plot(self.RApArcmin, filtMap[filterType])
         #
         plt.show()

      return filtMap, filtMask, filtHitNoiseStdDev, filtArea



   def saveFiltering(self, nProc=1):
      
      print("Evaluate all filters on all objects")
      # loop over all objects in catalog
#      result = np.array(map(self.analyzeObject, range(self.Catalog.nObj)))
      tStart = time()
      with sharedmem.MapReduce(np=nProc) as pool:
         f = lambda iObj: self.analyzeObject(iObj, test=False)
         result = np.array(pool.map(f, range(self.Catalog.nObj)))
      tStop = time()
      print "took", (tStop-tStart)/60., "min"


      # unpack and save to file
      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]

         filtMap = np.array([result[iObj,0][filterType][:] for iObj in range(self.Catalog.nObj)])
         filtMask = np.array([result[iObj,1][filterType][:] for iObj in range(self.Catalog.nObj)])
         filtHitNoiseStdDev = np.array([result[iObj,2][filterType][:] for iObj in range(self.Catalog.nObj)])
         filtArea = np.array([result[iObj,3][filterType][:] for iObj in range(self.Catalog.nObj)])
         
         np.savetxt(self.pathOut+"/"+filterType+"_filtmap.txt", filtMap)
         np.savetxt(self.pathOut+"/"+filterType+"_filtmask.txt", filtMask)
         np.savetxt(self.pathOut+"/"+filterType+"_filtnoisestddev.txt", filtHitNoiseStdDev)
         np.savetxt(self.pathOut+"/"+filterType+"_filtarea.txt", filtArea)


   def loadFiltering(self):
      self.filtMap = {}
      self.filtMask = {}
      self.filtHitNoiseStdDev = {}
      self.filtArea = {}

      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]
         self.filtMap[filterType] = np.genfromtxt(self.pathOut+"/"+filterType+"_filtmap.txt")
         self.filtMask[filterType] = np.genfromtxt(self.pathOut+"/"+filterType+"_filtmask.txt")
         self.filtHitNoiseStdDev[filterType] = np.genfromtxt(self.pathOut+"/"+filterType+"_filtnoisestddev.txt")
         self.filtArea[filterType] = np.genfromtxt(self.pathOut+"/"+filterType+"_filtarea.txt")


   ##################################################################################


   def catalogMask(self, overlap=True, psMask=True, mVir=None, z=[0., 100.], extraSelection=1., filterType=None, outlierReject=True):
      '''Returns catalog mask: 1 for objects to keep, 0 for objects to discard.
      Use as:
      maskedQuantity = Quantity[mask]
      '''
      if mVir is None:
         mVir = [self.mMin, self.mMax]

      # Here mask is 1 for objects we want to keep
      mask = np.ones_like(self.Catalog.RA)
      print "start with fraction", np.sum(mask)/len(mask), "of objects"
      if mVir is not None:
         mask *= (self.Catalog.Mvir>=mVir[0]) * (self.Catalog.Mvir<=mVir[1])
         print "keeping fraction", np.sum(mask)/len(mask), "of objects after mass cut"
      if z is not None:
         mask *= (self.Catalog.Z>=z[0]) * (self.Catalog.Z<=z[1])
         print "keeping fraction", np.sum(mask)/len(mask), "of objects after further z cut"
      if overlap:
         mask *= self.overlapFlag.copy()
         print "keeping fraction", np.sum(mask)/len(mask), "of objects after further overlap cut"
      # PS mask: look at largest aperture, and remove if any point within the disk or ring is masked
      if psMask:
         # The point source mask may vary from one filterType to another
         if filterType is None:
            filterType = self.filtMask.keys()[0]
         mask *= 1.*(np.abs(self.filtMask[filterType][:,-1])<1.)
         print "keeping fraction", np.sum(mask)/len(mask), "of objects after PS mask"
      mask *= extraSelection
      #print "keeping fraction", np.sum(mask)/len(mask), " of objects"
      if outlierReject:
         mask = mask.astype(bool)
         # we reject objects whose filter values are such
         # that the probability to have >=1 object with such a high absolute value
         # in a sample of this size
         # is equivalent to a 5sigma PTE = 5.73e-7.
         # Sample of 1: proba that the value is found outside of [-n*sigma, +n*sigma]
         # called PTE
         # is p = erf( n / sqrt(2))
         # Sample of N independent objects: the proba that at least one of them is outside [-n*sigma, +n*sigma]
         # is 1 - (1-p)^N
         # ie N * p  if p<<1
         # Our threshold: n*p should correspond to a 5 sigma PTE, ie
         # N * p = erf(5/sqrt(2)) = 5.733e-7
         nObj = np.sum(mask)
         if nObj<>0:
            f = lambda nSigmas: nObj * special.erfc(nSigmas / np.sqrt(2.)) - special.erfc(5. / np.sqrt(2.))
            nSigmasCut = optimize.brentq(f , 0., 1.e2)
            # ts has shape (nObj, nRAp)
            # sigmas has shape  nRAp
            sigmas = np.std(self.filtMap[filterType][mask,:], axis=0)
            # shape is (nObj, nRAp)
            newMask = (np.abs(self.filtMap[filterType][:,:]) <= nSigmasCut * sigmas[np.newaxis,:])
            # take the intersection of the masks
            mask *= np.prod(newMask, axis=1).astype(bool)
            print "keeping fraction", np.sum(1.*mask)/len(mask), "of objects after further outlier cut"
      # make sure the mask is boolean
      mask = mask.astype(bool)
      print "keeping fraction", np.sum(1.*mask)/len(mask), "in the end"
      return mask


   ##################################################################################


   def measureVarFromHitCount(self, filterType,  plot=False):
      """Returns a list of functions, one for each AP filter radius,
      where the function takes filtHitNoiseStdDev**2 \propto [(map var) * sr^2] as input and returns the
      actual measured filter variance [(map unit)^2 * sr^2].
      The functions are expected to be linear if the detector noise is the main source of noise,
      and if the hit counts indeed reflect the detector noise.
      To be used for noise weighting in the stacking.
      """
      # keep only objects that overlap, and mask point sources
      mask = self.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=(self.Catalog.Mvir.min(), self.Catalog.Mvir.max()), outlierReject=False)
      # This array contains the true variances for each object and aperture 
      filtVarTrue = np.zeros((self.Catalog.nObj, self.nRAp))

      if self.cmbHit is not None:
         print("Interpolate variance=f(hit count) for each aperture")
         fVarFromHitCount = np.empty(self.nRAp, dtype='object')
         for iRAp in range(self.nRAp):
            #print("Aperture number "+str(iRAp))
            x = self.filtHitNoiseStdDev[filterType][mask, iRAp]**2
            y = self.filtMap[filterType][mask, iRAp].copy()
            y = (y - np.mean(y))**2

            # define bins of hit count values
            nBins = 10  #21
            binEdges = np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), nBins, 10.)

#            # define bins of hit count values,
#            # with an equal number of objects in each bin
#            nBins = 10
#            binEdges = splitBins(x, nBins)
#            print "bin edges"
#            print binEdges
#            print np.min(x), np.max(x)
            
            # compute histograms
            binCenters, binEdges, binIndices = stats.binned_statistic(x, x, statistic='mean', bins=binEdges)
            binCounts, binEdges, binIndices = stats.binned_statistic(x, x, statistic='count', bins=binEdges)
            binnedVar, binEdges, binIndices = stats.binned_statistic(x, y, statistic=np.mean, bins=binEdges)
            sBinnedVar, binEdges, binIndices = stats.binned_statistic(x, y, statistic=np.std, bins=binEdges)
            sBinnedVar /= np.sqrt(binCounts)

            # exclude the empty bins, which did not contain data
            I = np.where(np.isfinite(binCenters * binCounts * binnedVar * sBinnedVar))[0]
            binCenters = binCenters[I]
            binCounts = binCounts[I]
            binnedVar = binnedVar[I]
            sBinnedVar = sBinnedVar[I]
            
            # interpolate, to use as noise weighting
            fVarFromHitCount[iRAp] = interp1d(binCenters, binnedVar, kind='linear', bounds_error=False, fill_value=(binnedVar[0],binnedVar[-1])) 

            # evaluate the variance for each object
            filtVarTrue[mask,iRAp] = fVarFromHitCount[iRAp](x)

#            # evaluate the variance for each object
#            # for the objects masked, still give them a weight just in case,
#            # to avoid null weights
#            xFull = self.filtHitNoiseStdDev[filterType][:, iRAp]**2
#            filtVarTrue[:,iRAp] = fVarFromHitCount[iRAp](xFull)

            
            if plot:
               # plot
               fig=plt.figure(0)
               ax=fig.add_subplot(111)
               #
               # measured
               ax.errorbar(binCenters, binnedVar*(180.*60./np.pi)**2, yerr=sBinnedVar*(180.*60./np.pi)**2, fmt='.', label=r'measured')
               # interpolated
               newX = np.logspace(np.log10(np.min(x)/2.), np.log10(np.max(x)*2.), 10.*nBins, 10.)
               newY = np.array(map(fVarFromHitCount[iRAp], newX))
               ax.plot(newX, newY*(180.*60./np.pi)**2, label=r'interpolated')
               #
               ax.set_xscale('log', nonposx='clip')
               ax.set_yscale('log', nonposy='clip')
               ax.set_xlabel(r'Det. noise var. from combined hit [arbitrary]')
               ax.set_ylabel(r'Measured var. [$\mu$K.arcmin$^2$]')
               #
               path = self.pathFig+"/binned_noise_vs_hit_"+filterType+"_"+str(iRAp)+".pdf"
               fig.savefig(path, bbox_inches='tight')
               fig.clf()

      else:
         print("Measure var for each aperture (no hit count)")
         meanVarAperture = np.var(self.filtMap[filterType][mask, :], axis=0)
         for iRAp in range(self.nRAp):
            filtVarTrue[mask,iRAp] = meanVarAperture[iRAp] * np.ones(np.sum(mask))

      return filtVarTrue


   def measureAllVarFromHitCount(self, plot=False):
      self.filtVarTrue = {}
      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]
         self.filtVarTrue[filterType] = self.measureVarFromHitCount(filterType, plot=plot)


   ##################################################################################


   def measureMeanTZBins(self, filterType, plot=False, test=False):

      '''Measure the mean filter temperatures in redshift bins,
      in order to subtract it to make the kSZ estimator robust to dust evolution * mean velocities.
      '''
      print("Measure mean T in z-bins (to subtract for kSZ)")

      # keep only objects that overlap, and mask point sources
      mask = self.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=(self.Catalog.Mvir.min(), self.Catalog.Mvir.max()))

      # redshift bins
      nZBins = 10
      zBinEdges = np.linspace(self.Catalog.Z.min(), self.Catalog.Z.max(), nZBins)

      # array of mean t to fill
      meanT = np.zeros((self.Catalog.nObj, self.nRAp))
      fMeanT = np.empty(self.nRAp, dtype='object')
      for iRAp in range(self.nRAp):
         # quantities to be binned
         z = self.Catalog.Z[mask].copy()
         t = self.filtMap[filterType][mask, iRAp].copy()
   #      s2 = self.filtVarTrue[filterType][mask, iRAp].copy()
   #      # inverse-variance weight the temperatures
   #      t  = t/s2 / np.mean(1./s2)
         # compute histograms
         zBinCenters, zBinEdges, zBinIndices = stats.binned_statistic(z, z, statistic='mean', bins=zBinEdges)
         zBinCounts, zBinEdges, zBinIndices = stats.binned_statistic(z, z, statistic='count', bins=zBinEdges)
         tMean, zBinEdges, zBinIndices = stats.binned_statistic(z, t, statistic=np.mean, bins=zBinEdges)
         sTMean, zBinEdges, zBinIndices = stats.binned_statistic(z, t, statistic=np.std, bins=zBinEdges)
         sTMean /= np.sqrt(zBinCounts)
         # piecewise-constant interpolation
         fMeanT[iRAp] = interp1d(zBinEdges[:-1], tMean, kind='previous', bounds_error=False, fill_value=(tMean[0], tMean[-1]))

         # evaluate the interpolation at each object
         meanT[mask,iRAp] = fMeanT[iRAp](z)

         if plot:
            fig=plt.figure(0)
            ax=fig.add_subplot(111)
            factor = (180.*60./np.pi)**2
            #
            ax.axhline(0., color='k')
            # Measured
            ax.errorbar(zBinCenters, factor*tMean, yerr=factor*sTMean, label=r'Measured')
            # Check interpolation
            z = np.linspace(self.Catalog.Z.min(), self.Catalog.Z.max(), 101)
            ax.plot(z, factor*fMeanT[iRAp](z), '--', label=r'Interpolated')
            #
            ax.legend(loc=1, fontsize='x-small', labelspacing=0.1)
            ax.set_xlabel(r'$z$ bins')
            ax.set_ylabel(r'Mean T per $z$-bin [$\mu$K$\cdot$arcmin$^2$]')
            #
            path = self.pathFig+"/binned_mean_t_"+filterType+"_"+str(iRAp)+".pdf"
            fig.savefig(path, bbox_inches='tight')
            if test:
               plt.show()
            fig.clf()

      return meanT



   def measureAllMeanTZBins(self, plot=False, test=False):
      self.meanT = {}
      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]
         print("For "+filterType+" filter:")
         self.meanT[filterType] = self.measureMeanTZBins(filterType, plot=plot, test=test)



   ##################################################################################


   def computeStackedProfile(self, filterType, est, iBootstrap=None, iVShuffle=None, tTh='', stackedMap=False, mVir=None, z=[0., 100.], ts=None, mask=None):
      """Returns the estimated profile and its uncertainty for each aperture.
      est: string to select the estimator
      iBootstrap: index for bootstrap resampling
      iVShuffle: index for shuffling velocities
      tTh: to replace measured temperatures by a theory expectation
      ts: option to specify another thumbstack object
      """

      #tStart = time()

      print("- Compute stacked profile: "+filterType+", "+est+", "+tTh)

      # compute stacked profile from another thumbstack object
      if ts is None:
         ts = self
      if mVir is None:
         mVir = [ts.mMin, ts.mMax]

      # select objects that overlap, and reject point sources
      if mask is None:
         mask = ts.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=mVir, z=z)

#      tMean = ts.meanT[filterType].copy()

      # temperatures [muK * sr]
      if tTh=='':
         t = ts.filtMap[filterType].copy() # [muK * sr]
      elif tTh=='tsz':
         # expected tSZ signal
         # AP profile shape, between 0 and 1
         sigma_cluster = 3.   #1.5  # arcmin
         shape = ts.ftheoryGaussianProfile(sigma_cluster) # between 0 and 1 [dimless]
         # multiply by integrated y to get y profile [sr]
         t = np.column_stack([ts.Catalog.integratedY[:] * shape[iAp] for iAp in range(ts.nRAp)])
         # convert from y profile to dT profile if needed
         if self.cmbUnitLatex==r'$\mu$K':
            nu = self.cmbNu   # Hz
            Tcmb = 2.726   # K
            h = 6.63e-34   # SI
            kB = 1.38e-23  # SI
            def f(nu):
               """frequency dependence for tSZ temperature
               """
               x = h*nu/(kB*Tcmb)
               return x*(np.exp(x)+1.)/(np.exp(x)-1.) -4.
            #t *= 2. * f(nu) * Tcmb * 1.e6  # [muK * sr]
            t *= f(nu) * Tcmb * 1.e6  # [muK * sr]
      elif tTh=='ksz':
         # expected kSZ signal
         # AP profile shape, between 0 and 1
         sigma_cluster = 1.5  # arcmin
         shape = ts.ftheoryGaussianProfile(sigma_cluster) # between 0 and 1 [dimless]
         # multiply by integrated kSZ to get kSZ profile [muK * sr]
         t = np.column_stack([ts.Catalog.integratedKSZ[:] * shape[iAp] for iAp in range(ts.nRAp)])   # [muK * sr]
         if self.cmbUnitLatex=='':
            t /= 2.726e6   # convert from [muK*sr] to [sr]
      t = t[mask, :]
 #     tMean = tMean[mask,:]
      # -v/c [dimless]
      v = -ts.Catalog.vR[mask] / 3.e5
      v -= np.mean(v)

#      # expected sigma_{v_{true}}, for the normalization
#      #print "computing v1d norm"
#      #tStartV = time()
#      z = ts.Catalog.Z[mask]
#      #f = lambda zGal: ts.U.v1dRms(0., zGal, W3d_sth)**2
#      #sVTrue = np.sqrt(np.mean(np.array(map(f, z))))
#      sVTrue = ts.U.v1dRms(0., np.mean(z), W3d_sth) / 3.e5  # (v^true_rms/c) [dimless]
#      #print "sigma_v_true =", sVTrue
#      #print "at z=0.57, expect", np.sqrt(f(0.57))
#      #tStopV = time()
#      #print "v1d norm took", tStopV - tStartV, "sec"
      
      #true filter variance for each object and aperture,
      # valid whether or not a hit count map is available
      s2Full = ts.filtVarTrue[filterType][mask, :]
      # Variance from hit count (if available)
      s2Hit = ts.filtHitNoiseStdDev[filterType][mask, :]**2
      #print "Shape of s2Hit = ", s2Hit.shape
      # halo masses
      m = ts.Catalog.Mvir[mask]

      if iBootstrap is not None:
         # make sure each resample is independent,
         # and make the resampling reproducible
         np.random.seed(iBootstrap)
         # list of overlapping objects
         nObj = np.sum(mask)
         #print "sample "iBootstrap, ";", nObj, "objects overlap with", ts.name
         I = np.arange(nObj)
         # choose with replacement from this list
         J = np.random.choice(I, size=nObj, replace=True)
         #
         t = t[J,:]
         #tMean = tMean[J,:]
         v = v[J]
         s2Hit = s2Hit[J,:]
         s2Full = s2Full[J,:]
         m = m[J]

      if iVShuffle is not None:
         # make sure each shuffling is independent,
         # and make the shuffling reproducible
         np.random.seed(iVShuffle)
         # list of overlapping objects
         nObj = np.sum(mask)
         I = np.arange(nObj)
         # shuffle the velocities
         J = np.random.permutation(I)
         #
         v = v[J]

      # tSZ: uniform weighting
      if est=='tsz_uniformweight':
         weights = np.ones_like(s2Hit)
         norm = 1./np.sum(weights, axis=0)
      # tSZ: detector-noise weighted (hit count)
      elif est=='tsz_hitweight':
         weights = 1./s2Hit
         norm = 1./np.sum(weights, axis=0)
      # tSZ: full noise weighted (detector noise + CMB)
      elif est=='tsz_varweight':
         weights = 1./s2Full
         norm = 1./np.sum(weights, axis=0)

      # kSZ: uniform weighting
      elif est=='ksz_uniformweight':
         # remove mean temperature
         #t -= np.mean(t, axis=0)
#         t -= tMean
         weights = v[:,np.newaxis] * np.ones_like(s2Hit)
         #norm = sVTrue / np.sum(v[:,np.newaxis]*weights, axis=0)
         norm = np.std(v) / ts.Catalog.rV / np.sum(v[:,np.newaxis]*weights, axis=0)
      # kSZ: detector-noise weighted (hit count)
      elif est=='ksz_hitweight':
         # remove mean temperature
         #t -= np.mean(t, axis=0)
#         t -= tMean
         weights = v[:,np.newaxis] / s2Hit
         #norm = sVTrue / np.sum(v[:,np.newaxis]*weights, axis=0)
         norm = np.std(v) / ts.Catalog.rV / np.sum(v[:,np.newaxis]*weights, axis=0)
      # kSZ: full noise weighted (detector noise + CMB)
      elif est=='ksz_varweight':
         # remove mean temperature
         #t -= np.mean(t, axis=0)
#         t -= tMean
         weights = v[:,np.newaxis] / s2Full
         #norm = sVTrue / np.sum(v[:,np.newaxis]*weights, axis=0)
         norm = np.std(v) / ts.Catalog.rV / np.sum(v[:,np.newaxis]*weights, axis=0)
      # kSZ: full noise weighted (detector noise + CMB)
      elif est=='ksz_massvarweight':
         # remove mean temperature
         #t -= np.mean(t, axis=0)
#         t -= tMean
         weights = m[:,np.newaxis] * v[:,np.newaxis] / s2Full
         #norm = np.mean(m) * sVTrue / np.sum(m[:,np.newaxis]**2 * v[:,np.newaxis]**2 / s2Full, axis=0)
         norm = np.mean(m) * np.std(v) / ts.Catalog.rV / np.sum(m[:,np.newaxis]**2 * v[:,np.newaxis]**2 / s2Full, axis=0)

      #tStop = time()
      #print "stacked profile took", tStop-tStart, "sec"


      # return the stacked profiles
      if not stackedMap:
         stack = norm * np.sum(t * weights, axis=0)
         sStack = norm * np.sqrt(np.sum(s2Full * weights**2, axis=0))
         return stack, sStack


      # or, if requested, compute and return the stacked cutout map
      else:
         # define chunks
         nChunk = ts.nProc
         chunkSize = ts.Catalog.nObj / nChunk
         # list of indices for each of the nChunk chunks
         chunkIndices = [range(iChunk*chunkSize, (iChunk+1)*chunkSize) for iChunk in range(nChunk)]
         # make sure not to miss the last few objects:
         # add them to the last chunk
         chunkIndices[-1] = range((nChunk-1)*chunkSize, ts.Catalog.nObj)

         # select weights for a typical aperture size (not the smallest, not the largest)
         #iRAp0 = ts.nRAp / 2
         iRAp0 = ts.nRAp / 4
         norm = norm[iRAp0]
         # need to link object number with weight,
         # despite the mask
         weightsLong = np.zeros(ts.Catalog.nObj)
         weightsLong[mask] = weights[:,iRAp0]

         def stackChunk(iChunk):
            # object indices to be processed
            chunk = chunkIndices[iChunk]

            # start with a null map for stacking
            resMap = ts.cutoutGeometry()
            for iObj in chunk:
               if iObj%10000==0:
                  print "- analyze object", iObj
               if ts.overlapFlag[iObj]:
                  # Object coordinates
                  ra = ts.Catalog.RA[iObj]   # in deg
                  dec = ts.Catalog.DEC[iObj] # in deg
                  z = ts.Catalog.Z[iObj]
                  # extract postage stamp around it
                  opos, stampMap, stampMask, stampHit = ts.extractStamp(ra, dec, test=False)
                  resMap += stampMap * weightsLong[iObj]
            return resMap

         # dispatch each chunk of objects to a different processor
         with sharedmem.MapReduce(np=ts.nProc) as pool:
            resMap = np.array(pool.map(stackChunk, range(nChunk)))

         # sum all the chunks
         resMap = np.sum(resMap, axis=0)
         # normalize by the proper sum of weights
         resMap *= norm

         return resMap








#   def computeStackedProfile(self, filterType, est, iBootstrap=None, iVShuffle=None, tTh=None, stackedMap=False, mVir=None, z=[0., 100.]):
#      """Returns the estimated profile and its uncertainty for each aperture.
#      est: string to select the estimator
#      iBootstrap: index for bootstrap resampling
#      iVShuffle: index for shuffling velocities
#      tTh: to replace measured temperatures by a theory expectation
#      """
#      if mVir is None:
#         mVir = [self.mMin, self.mMax]
#
#      # select objects that overlap, and reject point sources
#      mask = self.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=mVir, z=z)
#      
#      # temperatures [muK * sr]
#      if tTh is None:
#         t = self.filtMap[filterType].copy() # [muK * sr]
#      elif tTh=='tsz':
#         # expected tSZ signal
#         # AP profile shape, between 0 and 1
#         sigma_cluster = 1.5  # arcmin
#         shape = self.ftheoryGaussianProfile(sigma_cluster) # between 0 and 1 [dimless]
#         # multiply by integrated y to get y profile [sr]
#         t = np.column_stack([self.Catalog.integratedY[:] * shape[iAp] for iAp in range(self.nRAp)])
#         # convert from y profile to dT profile
#         Tcmb = 2.726   # K
#         h = 6.63e-34   # SI
#         kB = 1.38e-23  # SI
#         def f(nu):
#            """frequency dependence for tSZ temperature
#            """
#            x = h*nu/(kB*Tcmb)
#            return x*(np.exp(x)+1.)/(np.exp(x)-1.) -4.
#         t *= 2. * f(150.e9) * Tcmb * 1.e6  # [muK * sr]
#      elif tTh=='ksz':
#         # expected kSZ signal
#         # AP profile shape, between 0 and 1
#         sigma_cluster = 1.5  # arcmin
#         shape = self.ftheoryGaussianProfile(sigma_cluster) # between 0 and 1 [dimless]
#         # multiply by integrated kSZ to get kSZ profile [muK * sr]
#         t = np.column_stack([self.Catalog.integratedKSZ[:] * shape[iAp] for iAp in range(self.nRAp)])   # [muK * sr]
#      t = t[mask, :]
#      # -v/c [dimless]
#      v = -self.Catalog.vR[mask] / 3.e5
#      v -= np.mean(v)
#      #true filter variance for each object and aperture,
#      # valid whether or not a hit count map is available
#      s2Full = self.filtVarTrue[filterType][mask, :]
#      # Variance from hit count (if available)
#      s2Hit = self.filtHitNoiseStdDev[filterType][mask, :]**2
#      #print "Shape of s2Hit = ", s2Hit.shape
#      # halo masses
#      m = self.Catalog.Mvir[mask]
#      
#      if iBootstrap is not None:
#         # make sure each resample is independent,
#         # and make the resampling reproducible
#         np.random.seed(iBootstrap)
#         # list of overlapping objects
#         nObj = np.sum(mask)
#         I = np.arange(nObj)
#         # choose with replacement from this list
#         J = np.random.choice(I, size=nObj, replace=True)
#         #
#         t = t[J,:]
#         v = v[J]
#         s2Hit = s2Hit[J,:]
#         s2Full = s2Full[J,:]
#         m = m[J]
#      
#      if iVShuffle is not None:
#         # make sure each shuffling is independent,
#         # and make the shuffling reproducible
#         np.random.seed(iVShuffle)
#         # list of overlapping objects
#         nObj = np.sum(mask)
#         I = np.arange(nObj)
#         # shuffle the velocities
#         J = np.random.permutation(I)
#         #
#         v = v[J]
#      
#      # tSZ: uniform weighting
#      if est=='tsz_uniformweight':
#         weights = np.ones_like(s2Hit)
#         norm = 1./np.sum(weights, axis=0)
#      # tSZ: detector-noise weighted (hit count)
#      elif est=='tsz_hitweight':
#         weights = 1./s2Hit
#         norm = 1./np.sum(weights, axis=0)
#      # tSZ: full noise weighted (detector noise + CMB)
#      elif est=='tsz_varweight':
#         weights = 1./s2Full
#         norm = 1./np.sum(weights, axis=0)
#
#      # kSZ: uniform weighting
#      elif est=='ksz_uniformweight':
#         # remove mean temperature
#         t -= np.mean(t, axis=0)
#         weights = v[:,np.newaxis] * np.ones_like(s2Hit)
#         norm = np.std(v) / np.sum(v[:,np.newaxis]*weights, axis=0)
#      # kSZ: detector-noise weighted (hit count)
#      elif est=='ksz_hitweight':
#         # remove mean temperature
#         t -= np.mean(t, axis=0)
#         weights = v[:,np.newaxis] / s2Hit
#         norm = np.std(v) / np.sum(v[:,np.newaxis]*weights, axis=0)
#      # kSZ: full noise weighted (detector noise + CMB)
#      elif est=='ksz_varweight':
#         # remove mean temperature
#         t -= np.mean(t, axis=0)
#         weights = v[:,np.newaxis] / s2Full
#         norm = np.std(v) / np.sum(v[:,np.newaxis]*weights, axis=0)
#      # kSZ: full noise weighted (detector noise + CMB)
#      elif est=='ksz_massvarweight':
#         # remove mean temperature
#         t -= np.mean(t, axis=0)
#         weights = m[:,np.newaxis] * v[:,np.newaxis] / s2Full
#         norm = np.mean(m) * np.std(v) / np.sum(m[:,np.newaxis]**2 * v[:,np.newaxis]**2 / s2Full, axis=0)
#
#      # return the stacked profiles
#      if not stackedMap:
#         stack = norm * np.sum(t * weights, axis=0)
#         sStack = norm * np.sqrt(np.sum(s2Full * weights**2, axis=0))
#         return stack, sStack
#
#
#      # or, if requested, compute and return the stacked cutout map
#      else:
#         # define chunks
#         nChunk = self.nProc
#         chunkSize = self.Catalog.nObj / nChunk
#         # list of indices for each of the nChunk chunks
#         chunkIndices = [range(iChunk*chunkSize, (iChunk+1)*chunkSize) for iChunk in range(nChunk)]
#         # make sure not to miss the last few objects:
#         # add them to the last chunk
#         chunkIndices[-1] = range((nChunk-1)*chunkSize, self.Catalog.nObj)
#
#         # select weights for a typical aperture size (not the smallest, not the largest)
#         iRAp0 = self.nRAp / 2
#         norm = norm[iRAp0]
#         # need to link object number with weight,
#         # despite the mask
#         weightsLong = np.zeros(self.Catalog.nObj)
#         weightsLong[mask] = weights[:,iRAp0]
#
#         def stackChunk(iChunk):
#            # object indices to be processed
#            chunk = chunkIndices[iChunk]
#
#            # start with a null map for stacking 
#            resMap = self.cutoutGeometry()
#            for iObj in chunk:
#               if iObj%10000==0:
#                  print "- analyze object", iObj
#               if self.overlapFlag[iObj]:
#                  # Object coordinates
#                  ra = self.Catalog.RA[iObj]   # in deg
#                  dec = self.Catalog.DEC[iObj] # in deg
#                  z = self.Catalog.Z[iObj]
#                  # extract postage stamp around it
#                  opos, stampMap, stampMask, stampHit = self.extractStamp(ra, dec, test=False)
#                  resMap += stampMap * weightsLong[iObj]
#            return resMap
#
#         # dispatch each chunk of objects to a different processor
#         with sharedmem.MapReduce(np=self.nProc) as pool:
#            resMap = np.array(pool.map(stackChunk, range(nChunk)))
#
#         # sum all the chunks
#         resMap = np.sum(resMap, axis=0)
#         # normalize by the proper sum of weights
#         resMap *= norm
#         return resMap


   ##################################################################################



   def plotFilterHistograms(self, filterType, mVir=None, z=[0., 100.], ts=None):
      '''Plot histograms of the AP filter values,
      for each AP filter size,
      and compare with the Gaussian expectation.
      This enables checking for outliers.
      '''
      if ts is None:
         ts = self
      if mVir is None:
         mVir = [ts.mMin, ts.mMax]

      # select objects that overlap, and reject point sources
      mask = ts.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=mVir, z=z)

      # aperture photometry filters
      t = ts.filtMap[filterType].copy() # [muK * sr]
      t = t[mask, :]
      # convert from sr to arcmin^2
      factor = (180.*60./np.pi)**2
      t *= factor

      #true filter variance for each object and aperture,
      # valid whether or not a hit count map is available
      s2Full = ts.filtVarTrue[filterType][mask, :]
      # Variance from hit count (if available)
      s2Hit = ts.filtHitNoiseStdDev[filterType][mask, :]**2


      for iRAp in range(ts.nRAp):
         # Raw temperatures
         path = ts.pathFig + "/histogram_t_"+filterType+"_uniformweight_"+str(iRAp)+".pdf"
         myHistogram(t[:,iRAp], nBins=101, lim=None, S2Theory=[], path=path, plot=False, nameLatex=r'$T_i$ [$\mu$K$\cdot$arcmin$^2$]', semilogx=False, semilogy=True, doGauss=True)
         
         # inverse variance weighted temperatures
         path = ts.pathFig + "/histogram_t_"+filterType+"_varweight_"+str(iRAp)+".pdf"
         myHistogram(t[:,iRAp] / s2Full[:,iRAp], nBins=101, lim=None, S2Theory=[], path=path, plot=False, nameLatex=r'$T_i/\sigma_{T, i}^2$', semilogx=False, semilogy=True, doGauss=True)



   ##################################################################################

   def SaveCovBootstrapStackedProfile(self, filterType, est, mVir=None, z=[0., 100.], nSamples=100, nProc=1):
      """Estimate covariance matrix for the stacked profile from bootstrap resampling
      """
      print "Performing", nSamples, "bootstrap resamples"
      if mVir is None:
         mVir = [self.mMin, self.mMax]
      tStart = time()
      with sharedmem.MapReduce(np=nProc) as pool:
         f = lambda iSample: self.computeStackedProfile(filterType, est, iBootstrap=iSample, mVir=mVir, z=z)
         result = np.array(pool.map(f, range(nSamples)))
         #result = np.array(map(f, range(nSamples)))
      tStop = time()
      print "took", (tStop-tStart)/60., "min"
      # unpack results
      stackSamples = result[:,0,:] # shape (nSamples, nRAp)
      #sStackSamples = result[:,1,:]
      # estimate cov
      covStack = np.cov(stackSamples, rowvar=False)
      # save it to file
      np.savetxt(self.pathOut+"/cov_"+filterType+"_"+est+"_bootstrap.txt", covStack)
      


   def SaveCovBootstrapTwoStackedProfiles(self, ts2, filterType, est, mVir=None, z=[0., 100.], nSamples=100, nProc=1):
      """Estimate the full covariance for two stacked profiles.
      These need to use the same galaxy catalog. The temperature maps can be different.
      Resamples only the objects in common, then rescales the cov assuming
      that the objects that are not in common are not statistically different
      from the objects in common.
      """
      print "Performing", nSamples, "bootstrap resamples"
      # for each resample, compute both profiles, and concatenate, before taking the cov
      if mVir is None:
         mVir = [self.mMin, self.mMax]
      tStart = time()

      # find the objects that overlap with both maps 
      # bootstrap resample only these objects,
      # then rescale the cov mat
      mask1 = self.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=mVir, z=z)
      mask2 = ts2.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=mVir, z=z)
      mask12 = mask1 * mask2
      n1 = np.sum(mask1)
      n2 = np.sum(mask2)
      n12 = np.sum(mask12)
      print "Objects that overlap with 1, 2, 1&2:", n1, n2, n12
      # build the rescaling matrix for the cov
      block = np.ones((self.nRAp, self.nRAp))
      f11 = 1. * n12 / n1 * block
      f12 = 1. * n12**2 / n1 / n2 * block
      f22 = 1. * n12 / n2 * block
      rescaleMat = np.block([[f11, f12], [f12, f22]])
      
      def f(iSample):
         #print "Joint Bootstrap", iSample
         prof1 = self.computeStackedProfile(filterType, est, iBootstrap=iSample, mVir=mVir, z=z, mask=mask12)
         prof2 = self.computeStackedProfile(filterType, est, iBootstrap=iSample, mVir=mVir, z=z, ts=ts2, mask=mask12)
         # concatenate the 2 profiles
         jointProf = np.concatenate((prof1[0], prof2[0]))
         return jointProf
         
      with sharedmem.MapReduce(np=nProc) as pool:
         stackSamples = np.array(pool.map(f, range(nSamples)))
         #result = np.array(map(f, range(nSamples)))
      tStop = time()
      print "took", (tStop-tStart)/60., "min"
      # estimate cov
      covStack = np.cov(stackSamples, rowvar=False)
      # rescale the cov mat, since only the galaxies in common for the two maps
      # were sampled
      covStack *= rescaleMat
      # save it to file
      path = self.pathOut+"/cov_"+filterType+"_"+est+"_joint_"+self.name+"_"+ts2.name+"_bootstrap.txt"
      print "saving to:"
      print path
      np.savetxt(path, covStack)
      

   def saveAllCovBootstrapTwoStackedProfiles(self, ts2):
      print "- compute full joint cov between stacked profiles from:"
      print self.name
      print ts2.name
      # Compute all filter types
      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]
         # covariance matrices from bootstrap,
         # only for a few select estimators
         for iEst in range(len(self.EstBootstrap)):
            est = self.EstBootstrap[iEst]
            self.SaveCovBootstrapTwoStackedProfiles(ts2, filterType, est, nSamples=self.nSamples, nProc=min(8,self.nProc))



   def SaveCovVShuffleStackedProfile(self, filterType, est, mVir=None, z=[0., 100.], nSamples=100, nProc=1):
      """Estimate covariance matrix for the stacked profile from shuffling velocities
      """
      if mVir is None:
         mVir = [self.mMin, self.mMax]
      tStart = time()
      with sharedmem.MapReduce(np=nProc) as pool:
         f = lambda iSample: self.computeStackedProfile(filterType, est, iVShuffle=iSample, mVir=mVir, z=z)
         result = np.array(pool.map(f, range(nSamples)))
      tStop = time()
      print "took", (tStop-tStart)/60., "min"
      # unpack results
      stackSamples = result[:,0,:] # shape (nSamples, nRAp)
      #sStackSamples = result[:,1,:]
      # estimate cov
      covStack = np.cov(stackSamples, rowvar=False)
      # save it to file
      np.savetxt(self.pathOut+"/cov_"+filterType+"_"+est+"_vshuffle.txt", covStack)

      # Null test: measure mean
      data = np.zeros((self.nRAp, 3))
      data[:,0] = self.RApArcmin
      data[:,1] = np.mean(stackSamples, axis=0) # mean of the shuffled stacks
      data[:,2] = np.sqrt(np.diag(covStack)) / np.sqrt(nSamples) # variance of the mean
      np.savetxt(self.pathOut+"/"+filterType+"_"+est+"_vshufflemean.txt", data)


   ##################################################################################


   def saveAllStackedMaps(self, filterTypes=None, Est=None):
      print "- compute all stacked maps"
      if filterTypes is None:
         filterTypes = self.filterTypes
      if Est is None:
         Est = self.Est

      # Get cutout geometry, for plotting
      cutoutMap = self.cutoutGeometry()
      size = cutoutMap.posmap()[0,:,:].max() - cutoutMap.posmap()[0,:,:].min()
      baseMap = FlatMap(nX=cutoutMap.shape[0], nY=cutoutMap.shape[1], sizeX=size, sizeY=size)


      # loop over filter types: only matter
      # because they determine the weights in the stacked map
      for iFilterType in range(len(filterTypes)):
         filterType = filterTypes[iFilterType]
         # Estimators (tSZ, kSZ, various weightings...)
         for iEst in range(len(Est)):
            est = Est[iEst]
            print "compute stacked map:", filterType, est
            stackedMap = self.computeStackedProfile(filterType, est, iBootstrap=None, iVShuffle=None, tTh='', stackedMap=True)

            # save the stacked cutout
            path = self.pathOut + "/stackedmap_"+filterType+"_"+est+".txt"
            np.savetxt(path, stackedMap)

            ## for low SNR signals, downgrade the map
            #lowBaseMap = FlatMap(nX=cutoutMap.shape[0]//2, nY=cutoutMap.shape[1]//2, sizeX=size, sizeY=size)
            #nXNew = stackedMap.shape[0] // 2
            #nYNew = stackedMap.shape[1] // 2
            #stackedMap = baseMap.downResolution(nXNew, nYNew, data=stackedMap)

            # plot the stacked map and save it
            path = self.pathFig + "/stackedmap_"+filterType+"_"+est+".pdf"
            baseMap.plot(data=stackedMap, save=True, path=path)
            #lowBaseMap.plot(data=stackedMap, save=True, path=path)


   ##################################################################################


   def saveAllStackedProfiles(self):
      print "- compute all stacked profiles and their cov"
      tStart = time()
      data = np.zeros((self.nRAp, 3))
      data[:,0] = self.RApArcmin # [arcmin]
      
      # Compute all filter types and estimators
      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]

         # check AP filter histograms
#         self.plotFilterHistograms(filterType)


         # Estimators (tSZ, kSZ, various weightings...)
         for iEst in range(len(self.Est)):
            est = self.Est[iEst]
            # measured stacked profile
            data[:,1], data[:,2] = self.computeStackedProfile(filterType, est) # [map unit * sr]
            np.savetxt(self.pathOut+"/"+filterType+"_"+est+"_measured.txt", data)
            # expected stacked profile from tSZ
            data[:,1], data[:,2] = self.computeStackedProfile(filterType, est, tTh='tsz') # [map unit * sr]
            np.savetxt(self.pathOut+"/"+filterType+"_"+est+"_theory_tsz.txt", data)
            # expected stacked profile from kSZ
            data[:,1], data[:,2] = self.computeStackedProfile(filterType, est, tTh='ksz') # [map unit * sr]
            np.savetxt(self.pathOut+"/"+filterType+"_"+est+"_theory_ksz.txt", data)

         # covariance matrices from bootstrap,
         # only for a few select estimators
         if self.doBootstrap: 
            for iEst in range(len(self.EstBootstrap)):
               est = self.EstBootstrap[iEst]
               self.SaveCovBootstrapStackedProfile(filterType, est, nSamples=self.nSamples, nProc=self.nProc)

         # covariance matrices from shuffling velocities,
         # for ksz only
         if self.doVShuffle:
            for iEst in range(len(self.EstVShuffle)):
               est = self.EstVShuffle[iEst]
               self.SaveCovVShuffleStackedProfile(filterType, est, nSamples=self.nSamples, nProc=self.nProc)


      # Stacked profiles in mass bins, to check for contamination
      if self.doMBins:
         for iFilterType in range(len(self.filterTypes)):
            filterType = self.filterTypes[iFilterType]
            for iEst in range(len(self.EstMBins)):
               est = self.EstMBins[iEst]
               data = np.zeros((self.nRAp, 2*self.nMMax+1))
               data[:,0] = self.RApArcmin # [arcmin]
               dataTsz = data.copy()
               dataKsz = data.copy()
               for iMMax in range(self.nMMax):
                  mMax = self.MMax[iMMax]

                  # measured stacked profile
                  data[:,1+2*iMMax], data[:,1+2*iMMax+1] = self.computeStackedProfile(filterType, est, mVir=[1.e6, mMax]) # [map unit * sr]
                  # expcted from tSZ
                  dataTsz[:,1+2*iMMax], dataTsz[:,1+2*iMMax+1] = self.computeStackedProfile(filterType, est, mVir=[1.e6, mMax], tTh='tsz') # [map unit * sr]
                  # expected from kSZ
                  dataKsz[:,1+2*iMMax], dataKsz[:,1+2*iMMax+1] = self.computeStackedProfile(filterType, est, mVir=[1.e6, mMax], tTh='ksz') # [map unit * sr]

               # Save all stacked profiles
               np.savetxt(self.pathOut+"/"+filterType+"_"+est+"_mmax_measured.txt", data)
               np.savetxt(self.pathOut+"/"+filterType+"_"+est+"_mmax_theory_tsz.txt", dataTsz)
               np.savetxt(self.pathOut+"/"+filterType+"_"+est+"_mmax_theory_ksz.txt", dataKsz)

      tStop = time()
      print "Computing all stacked profiles (and cov) took", (tStop-tStart)/60., "min"


   def loadAllStackedProfiles(self):
      print "- load stacked profiles and their cov"
      self.stackedProfile = {}
      self.sStackedProfile = {}
      self.covBootstrap = {}
      self.covVShuffle = {}
      
      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]

         # all stacked profiles
         for iEst in range(len(self.Est)):
            est = self.Est[iEst]
            # measured stacked profile
            data = np.genfromtxt(self.pathOut+"/"+filterType+"_"+est+"_measured.txt")
            self.stackedProfile[filterType+"_"+est] = data[:,1]
            self.sStackedProfile[filterType+"_"+est] = data[:,2]
            # expected stacked profile from tSZ
            data = np.genfromtxt(self.pathOut+"/"+filterType+"_"+est+"_theory_tsz.txt")
            self.stackedProfile[filterType+"_"+est+"_theory_tsz"] = data[:,1]
            self.sStackedProfile[filterType+"_"+est+"_theory_tsz"] = data[:,2]
            # expected stacked profile from kSZ
            data = np.genfromtxt(self.pathOut+"/"+filterType+"_"+est+"_theory_ksz.txt")
            self.stackedProfile[filterType+"_"+est+"_theory_ksz"] = data[:,1]
            self.sStackedProfile[filterType+"_"+est+"_theory_ksz"] = data[:,2]

         # Null tests from shuffling velocities,
         # for ksz only
         if self.doVShuffle:
            for iEst in range(len(self.EstVShuffle)):
               est = self.EstVShuffle[iEst]
               # null test from shuffling the velocities
               data = np.genfromtxt(self.pathOut+"/"+filterType+"_"+est+"_vshufflemean.txt")
               self.stackedProfile[filterType+"_"+est+"_vshufflemean"] = data[:,1]
               self.sStackedProfile[filterType+"_"+est+"_vshufflemean"] = data[:,2]



         # stacked profiles in mass bins
         if self.doMBins:
            for iEst in range(len(self.EstMBins)):
               est = self.EstMBins[iEst]

               # measured stacked profile
               data = np.genfromtxt(self.pathOut+"/"+filterType+"_"+est+"_mmax_measured.txt")
               for iMMax in range(self.nMMax):
                  self.stackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)] = data[:,1+2*iMMax]
                  self.sStackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)] = data[:,1+2*iMMax+1]
               # expected stacked profile from tSZ
               data = np.genfromtxt(self.pathOut+"/"+filterType+"_"+est+"_mmax_theory_tsz.txt")
               for iMMax in range(self.nMMax):
                  self.stackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)+"_theory_tsz"] = data[:,1+2*iMMax]
                  self.sStackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)+"_theory_tsz"] = data[:,1+2*iMMax+1]
               # expected stacked profile from kSZ
               data = np.genfromtxt(self.pathOut+"/"+filterType+"_"+est+"_mmax_theory_ksz.txt")
               for iMMax in range(self.nMMax):
                  self.stackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)+"_theory_ksz"] = data[:,1+2*iMMax]
                  self.sStackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)+"_theory_ksz"] = data[:,1+2*iMMax+1]

         # covariance matrices from bootstrap,
         # only for a few select estimators
         if self.doBootstrap:
            for iEst in range(len(self.EstBootstrap)):
               est = self.EstBootstrap[iEst]
               self.covBootstrap[filterType+"_"+est] = np.genfromtxt(self.pathOut+"/cov_"+filterType+"_"+est+"_bootstrap.txt")
         
         # covariance matrices from shuffling velocities,
         # for ksz only
         if self.doVShuffle:
            for iEst in range(len(self.EstVShuffle)):
               est = self.EstVShuffle[iEst]
               self.covVShuffle[filterType+"_"+est] = np.genfromtxt(self.pathOut+"/cov_"+filterType+"_"+est+"_vshuffle.txt")



   ##################################################################################

   def plotStackedProfile(self, filterType, Est, name=None, pathDir=None, theory=True, tsArr=None, plot=False, legend=True):
      """Compares stacked profiles, and their uncertainties.
      If pathDir is not specified, save to local figure folder.
      """
      if name is None:
         name = Est[0]
      if tsArr is None:
         tsArr = [self]
      if pathDir is None:
         pathDir = self.pathFig
      
      # stacked profile
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # convert from sr to arcmin^2
      factor = (180.*60./np.pi)**2
      #
      ax.axhline(0., c='k', lw=1)
      #
      #colors = ['r', 'g', 'b', 'm', 'c']
      lineStyles = ['-', '--', '-.', ':']
      for iEst in range(len(Est)):
         est = Est[iEst]
         #c = colors[iEst%len(colors)]

         for iTs in range(len(tsArr)):
            ts = tsArr[iTs]
            ls = lineStyles[iTs%len(tsArr)]

            ax.errorbar(ts.RApArcmin+iTs*0.05, factor * ts.stackedProfile[filterType+"_"+est], factor * ts.sStackedProfile[filterType+"_"+est], fmt=ls, label=filterType.replace('_',' ')+' '+est.replace('_', ' ')+' '+ts.name.replace('_',' '))
            if theory:
               ax.plot(ts.RApArcmin+iTs*0.05, factor * ts.stackedProfile[filterType+"_"+est+"_theory_tsz"], ls='--', label="theory tsz, "+filterType.replace('_',' ')+' '+est.replace('_', ' ')+' '+ts.name.replace('_',' '))
               ax.plot(ts.RApArcmin+iTs*0.05, factor * ts.stackedProfile[filterType+"_"+est+"_theory_ksz"], ls='-.', label="theory ksz, "+filterType.replace('_',' ')+' '+est.replace('_', ' ')+' '+ts.name.replace('_',' '))
      #
      if legend:
         ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'$T$ [$\mu K\cdot\text{arcmin}^2$]')
      #ax.set_ylim((0., 2.))
      #
      path = pathDir+"/"+name+".pdf"
      fig.savefig(path, bbox_inches='tight')
      if plot:
         plt.show()
      fig.clf()



      # uncertainty on stacked profile
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # convert from sr to arcmin^2
      factor = (180.*60./np.pi)**2
      #
      colors = ['r', 'g', 'b', 'm', 'c']
      lineStyles = ['-', '--', '-.', ':']
      for iEst in range(len(Est)):
         est = Est[iEst]
         c = colors[iEst%len(colors)]

         for iTs in range(len(tsArr)):
            ts = tsArr[iTs]
            ls = lineStyles[iTs%len(tsArr)]

            ax.plot(ts.RApArcmin+iTs*0.05, factor * ts.sStackedProfile[filterType+"_"+est], c=c, ls=ls, lw=2, label='analytic, '+filterType.replace('_',' ')+' '+est.replace('_', ' ')+' '+ts.name.replace('_',' '))
            if est in ts.covBootstrap:
               ax.plot(ts.RApArcmin+iTs*0.05, factor * np.sqrt(np.diag(ts.covBootstrap[filterType+"_"+est])), c=c, ls=ls, lw=1.5, label="bootstrap, "+filterType.replace('_',' ')+' '+est.replace('_', ' ')+' '+ts.name.replace('_',' '))
            if est in ts.covVShuffle:
               ax.plot(ts.RApArcmin+iTs*0.05, factor * np.sqrt(np.diag(ts.covVShuffle[filterType+"_"+est])), c=c, ls=ls, lw=1, label="v shuffle, "+filterType.replace('_',' ')+' '+est.replace('_', ' ')+' '+ts.name.replace('_',' '))
      #
      if legend:
         ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'$\sigma(T)$ [$\mu K\cdot\text{arcmin}^2$]')
      #ax.set_ylim((0., 2.))
      #
      path = pathDir+"/s_"+name+".pdf"
      fig.savefig(path, bbox_inches='tight')
      if plot:
         plt.show()
      fig.clf()


   def plotAllStackedProfiles(self):
      print "- plot all stacked profiles"
      
      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]
         # all stacked profiles
         for iEst in range(len(self.Est)):
            est = self.Est[iEst]
            self.plotStackedProfile(filterType, [est], name=filterType+"_"+est)

         # stacked profiles in mass bins
         if self.doMBins:
            self.plotTszKszContaminationMMax()

            for iEst in range(len(self.EstMBins)):
               est = self.EstMBins[iEst]
               # measured stacked profiles
               estArr = [est+"_mmax"+str(iMMax) for iMMax in range(self.nMMax)] + [est]
               self.plotStackedProfile(filterType, estArr, name=filterType+"_"+est+"_mmax", theory=False, legend=False)
#               # expected from tSZ
#               estArr = [est+"_mmax"+str(iMMax)+"_theory_tsz" for iMMax in range(self.nMMax)] + [est]
#               self.plotStackedProfile(filterType, estArr, name=filterType+"_"+est+"_mmax_theory_tsz", theory=False, legend=False)
#               # expected from kSZ
#               estArr = [est+"_mmax"+str(iMMax)+"_theory_ksz" for iMMax in range(self.nMMax)] + [est]
#               self.plotStackedProfile(filterType, estArr, name=filterType+"_"+est+"_mmax_theory_ksz", theory=False, legend=False)


   ##################################################################################

   def plotTszKszContaminationMMax(self):
      '''Bias from tSZ to kSZ and vice versa:
      compare bias to signal and error bar,
      as a function of maximum mass in the galaxy sample.
      '''
      print "Plotting contamination as a function of MMax"

      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]

         for iEst in range(len(self.EstMBins)):
            est = self.EstMBins[iEst]

            TszToKsz = np.zeros(self.nMMax)
            KszToTsz = np.zeros(self.nMMax)
            # checks that the estimators are unbiased for all Mmax
            KszToKsz = np.zeros(self.nMMax)
            TszToTsz = np.zeros(self.nMMax)
            # error bars as a function of Mmax
            sKsz = np.zeros(self.nMMax)
            sTsz = np.zeros(self.nMMax)
            for iMMax in range(self.nMMax):
               # ratio of expected tSZ for this Mmax
               # to the full expected kSZ
               ratio = self.stackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)+"_theory_tsz"].copy()
               ratio /= self.stackedProfile[filterType+"_"+est+"_theory_ksz"]
               # here the value at each aperture is equal to the mean
               TszToKsz[iMMax] = np.mean(ratio)

               # ratio of expected tSZ for this Mmax
               # to the full expected tSZ, as a check
               ratio = self.stackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)+"_theory_tsz"].copy()
               ratio /= self.stackedProfile[filterType+"_"+est+"_theory_tsz"]
               # here the value at each aperture is equal to the mean
               TszToTsz[iMMax] = np.mean(ratio)

               # expected kSZ error bar for the corresponding Mmax cut
               ratio = self.sStackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)].copy()
               ratio /= self.stackedProfile[filterType+"_"+est+"_theory_ksz"]
               #ratio /= self.stackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)+"_theory_ksz"]
               # here the value at each aperture is equal to the mean
               sKsz[iMMax] = np.abs(np.mean(ratio))

                  
               # ratio of expected kSZ for this Mmax
               # to the full expected tSZ
               ratio = self.stackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)+"_theory_ksz"].copy()
               ratio /= self.stackedProfile[filterType+"_"+est+"_theory_tsz"]
               # here the value at each aperture is equal to the mean
               KszToTsz[iMMax] = np.mean(ratio)

               # ratio of expected kSZ for this Mmax
               # to the full expected kSZ, as a check
               ratio = self.stackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)+"_theory_ksz"].copy()
               ratio /= self.stackedProfile[filterType+"_"+est+"_theory_ksz"]
               # here the value at each aperture is equal to the mean
               KszToKsz[iMMax] = np.mean(ratio)

               # expected error tSZ bar for the corresponding Mmax cut
               ratio = self.sStackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)].copy()
               ratio /= self.stackedProfile[filterType+"_"+est+"_theory_tsz"]
               #ratio /= self.stackedProfile[filterType+"_"+est+"_mmax"+str(iMMax)+"_theory_tsz"]
               # here the value at each aperture is equal to the mean
               sTsz[iMMax] = np.abs(np.mean(ratio))
   

            fig=plt.figure(0)
            ax=fig.add_subplot(111)
            #
            # convert from sr to arcmin^2
            factor = (180.*60./np.pi)**2
            #
            ax.axhline(0., c='k', lw=1)
            #
            #ax.axhline(1., lw=1, label=r'Expected kSZ')
            ax.plot(self.MMax, KszToKsz, 'b', label='kSZ signal')
            ax.plot(self.MMax, -KszToKsz, 'b--')
            ax.plot(self.MMax, 0.1*KszToKsz, 'b', alpha=0.3)
            ax.plot(self.MMax, -0.1*KszToKsz, 'b--', alpha=0.3)
            ax.plot([self.mMax], [1.], 'bo')
            #
            ax.plot(self.MMax, TszToKsz, 'r', label='tSZ bias to kSZ')
            ax.plot(self.MMax, -TszToKsz, 'r--')
            #
            ax.fill_between(self.MMax, sKsz, edgecolor='', facecolor='gray', alpha=0.5, label='kSZ error bar')
            ax.fill_between(self.MMax, 0.1*sKsz, edgecolor='', facecolor='gray', alpha=0.3)
            #
            # ignore the first mMax values for which there is no object in the sample
            xMin = np.min(self.MMax[np.isfinite(sKsz)])
            xMax = np.max(self.MMax[np.isfinite(sKsz)])
            ax.set_xlim((xMin, xMax))
            ax.set_ylim((1.e-3, 10.))
            ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
            ax.set_xscale('log', nonposx='clip')
            ax.set_yscale('log', nonposy='clip')
            ax.set_xlabel(r'$M_\text{vir max}$ [$M_\odot$]')
            ax.set_ylabel(r'Fraction of expected kSZ')
            #
            name = filterType+"_"+est+"_mmax_tsztoksz"
            path = self.pathFig+"/"+name+".pdf"
            fig.savefig(path, bbox_inches='tight')
            fig.clf()


            fig=plt.figure(1)
            ax=fig.add_subplot(111)
            #
            # convert from sr to arcmin^2
            factor = (180.*60./np.pi)**2
            #
            ax.axhline(0., c='k', lw=1)
            #ax.axhline(1., lw=1, label=r'Expected tSZ')
            ax.plot(self.MMax, TszToTsz, 'b', label='tSZ signal')
            ax.plot(self.MMax, -TszToTsz, 'b--')
            ax.plot(self.MMax, 0.1*TszToTsz, 'b', alpha=0.3)
            ax.plot(self.MMax, -0.1*TszToTsz, 'b--', alpha=0.3)
            ax.plot([self.mMax], [1.], 'bo')
            #
            ax.plot(self.MMax, KszToTsz, 'r', label='kSZ bias to tSZ')
            ax.plot(self.MMax, -KszToTsz, 'r--')
            #
            ax.fill_between(self.MMax, sTsz, edgecolor='', facecolor='gray', alpha=0.5, label='tSZ error bar')
            ax.fill_between(self.MMax, 0.1*sTsz, edgecolor='', facecolor='gray', alpha=0.3)
            #
            ax.set_xlim((self.MMax[1], self.MMax.max()))
            ax.set_ylim((1.e-3, 10.))
            ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
            ax.set_xscale('log', nonposx='clip')
            ax.set_yscale('log', nonposy='clip')
            ax.set_xlabel(r'$M_\text{vir max}$ [$M_\odot$]')
            ax.set_ylabel(r'Fraction of expected tSZ')
            #
            name = filterType+"_"+est+"_mmax_ksztotsz"
            path = self.pathFig+"/"+name+".pdf"
            fig.savefig(path, bbox_inches='tight')
            fig.clf()






   ##################################################################################
   
   def plotCov(self, cov, name="", show=False):
      # Covariance matrix
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # pcolor wants x and y to be edges of cell,
      # ie one more element, and offset by half a cell
      dR = (self.rApMaxArcmin - self.rApMinArcmin) / self.nRAp
      RApEdgesArcmin = np.linspace(self.rApMinArcmin-0.5*dR, self.rApMaxArcmin+0.5*dR, self.nRAp+1)
      # Covariance  matrix
      X, Y = np.meshgrid(RApEdgesArcmin, RApEdgesArcmin, indexing='ij')
      cp=ax.pcolormesh(X, Y, cov, cmap='YlOrRd')
      #
      ax.set_aspect('equal')
      plt.colorbar(cp)
      ax.set_xlim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_ylim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_xlabel(r'R [arcmin]')
      ax.set_ylabel(r'R [arcmin]')
      #
      path = self.pathFig+"/cov_"+name+".pdf"
      fig.savefig(path, bbox_inches='tight')
      if show:
         plt.show()
      else:
         fig.clf()
      

      # Correlation coefficient
      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      # pcolor wants x and y to be edges of cell,
      # ie one more element, and offset by half a cell
      dR = (self.rApMaxArcmin - self.rApMinArcmin) / self.nRAp
      RApEdgesArcmin = np.linspace(self.rApMinArcmin-0.5*dR, self.rApMaxArcmin+0.5*dR, self.nRAp+1)
      X, Y = np.meshgrid(RApEdgesArcmin, RApEdgesArcmin, indexing='ij')
      #
      sigma = np.sqrt(np.diag(cov))
      cor = np.array([[cov[i1, i2] / (sigma[i1]*sigma[i2]) for i2 in range(self.nRAp)] for i1 in range(self.nRAp)])
      cp=ax.pcolormesh(X, Y, cor, cmap='YlOrRd', vmin=0., vmax=1.)
      #
      ax.set_aspect('equal')
      plt.colorbar(cp)
      ax.set_xlim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_ylim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_xlabel(r'R [arcmin]')
      ax.set_ylabel(r'R [arcmin]')
      #
      path = self.pathFig+"/cor_"+name+".pdf"
      fig.savefig(path, bbox_inches='tight')
      if show:
         plt.show()
      else:
         fig.clf()


   def plotAllCov(self):
      print "- plot all covariances"
      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]

         # covariance matrices from bootstrap,
         # only for a few select estimators
         if self.doBootstrap:
            for iEst in range(len(self.EstBootstrap)):
               est = self.EstBootstrap[iEst]
               self.plotCov(self.covBootstrap[filterType+"_"+est], filterType+"_"+est+"_bootstrap")
         
         # covariance matrices from shuffling velocities,
         # for ksz only
         if self.doVShuffle:
            for iEst in range(len(self.EstVShuffle)):
               est = self.EstVShuffle[iEst]
               self.plotCov(self.covVShuffle[filterType+"_"+est], filterType+"_"+est+"_vshuffle")


   ##################################################################################


   def plotCovTwoStackedProfiles(self, cov, name="", show=False):
      '''For the figure captions, it is assumed that the first stack
      is at 150GHz, the second at 90GHz.
      '''

      RApArcmin = np.concatenate((self.RApArcmin, self.RApArcmin))
      X, Y = np.meshgrid(RApArcmin, RApArcmin, indexing='ij')
      #
      # corresponding indices (for masks)
      I = np.arange(len(RApArcmin))
      II, JJ = np.meshgrid(I, I, indexing='ij')


      # Covariance matrix
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # Covariance  matrix
      covLower = np.ma.masked_where(II<JJ, cov)
      cp=ax.imshow(covLower, cmap='YlOrRd')#vmin=0., vmax=1.
      plt.colorbar(cp)
      #
      ax.axhline(0.5 * I[-1], c='k')
      ax.axvline(0.5 * I[-1], c='k')
      #
      ax.text(0.25, 0.9, r'150 GHz', ha='center', va='center', transform=ax.transAxes, fontsize=14)
      ax.text(0.75, 0.75, r'150 GHz - 90 GHz', ha='center', va='center', transform=ax.transAxes, fontsize=14)
      ax.text(0.75, 0.4, r'90 GHz', ha='center', va='center', transform=ax.transAxes, fontsize=14)
      #
      ax.set_xticks(np.arange(len(RApArcmin)))
      ax.set_yticks(np.arange(len(RApArcmin)))
      ax.set_xticklabels(RApArcmin, fontdict={'fontsize': 10}, rotation=45)
      ax.set_yticklabels(RApArcmin, fontdict={'fontsize': 10})
      ax.set_aspect('equal')
      ax.set_xlabel(r'R [arcmin]')
      ax.set_ylabel(r'R [arcmin]')
      ax.set_title(r'Joint covariance matrix')
      #
      path = self.pathFig+"/cov_joint_"+name+".pdf"
      print "Saving cov plot to:"
      print path
      fig.savefig(path, bbox_inches='tight')
      if show:
         plt.show()
      else:
         fig.clf()

      # Correlation coefficients
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # Correlation coefficients
      sigma = np.sqrt(np.diag(cov))
      cor = np.array([[cov[i1, i2] / (sigma[i1]*sigma[i2]) for i2 in I] for i1 in I])
      corLower = np.ma.masked_where(II<JJ, cor)
      cp=ax.imshow(corLower, cmap='YlGnBu', vmin=0., vmax=1.)
      cb=plt.colorbar(cp)
      #
      ax.axhline(0.5 * I[-1], c='k')
      ax.axvline(0.5 * I[-1], c='k')
      #
      ax.text(0.25, 0.9, r'150 GHz', ha='center', va='center', transform=ax.transAxes, fontsize=14)
      ax.text(0.75, 0.75, r'150 GHz - 90 GHz', ha='center', va='center', transform=ax.transAxes, fontsize=14)
      ax.text(0.75, 0.4, r'90 GHz', ha='center', va='center', transform=ax.transAxes, fontsize=14)
      #
      ax.set_xticks(np.arange(len(RApArcmin)))
      ax.set_yticks(np.arange(len(RApArcmin)))
      ax.set_xticklabels(RApArcmin, fontdict={'fontsize': 10}, rotation=45)
      ax.set_yticklabels(RApArcmin, fontdict={'fontsize': 10})
      ax.set_aspect('equal')
      ax.set_xlabel(r'R [arcmin]')
      ax.set_ylabel(r'R [arcmin]')
      ax.set_title(r'Joint correlation matrix')
      #
      path = self.pathFig+"/cor_joint_"+name+".pdf"
      print "Saving cor plot to:"
      print path
      fig.savefig(path, bbox_inches='tight')
      if show:
         plt.show()
      else:
         fig.clf()



   def plotAllCovTwoStackedProfiles(self, ts2):
      print "- plot all full joint cov between stacked profiles from:"
      print self.name
      print ts2.name
      
      # Compute all filter types
      for iFilterType in range(len(self.filterTypes)):
         filterType = self.filterTypes[iFilterType]

         # covariance matrices from bootstrap,
         # only for a few select estimators
         for iEst in range(len(self.EstBootstrap)):
            est = self.EstBootstrap[iEst]

            # read cov from file
            pathCov = self.pathOut+"/cov_"+filterType+"_"+est+"_joint_"+self.name+"_"+ts2.name+"_bootstrap.txt"
            print "Reading joint cov from:"
            print pathCov
            cov = np.genfromtxt(pathCov)
            name = filterType+"_"+est+"_"+self.name+"_"+ts2.name+"_bootstrap"
            self.plotCovTwoStackedProfiles(cov, name=name, show=False)







   
   ##################################################################################

   def ftheoryGaussianProfile(self, sigma_cluster, filterType='diskring'):
      """Alpha_ksz signal, between 0 and 1.
      Assumes that the projected cluster profile is a 2d Gaussian,
      with sigma_cluster in arcmin
      Assumes equal area disk-ring filter
      """
      if filterType=='diskring':
         result = (1. - np.exp(-0.5*self.RApArcmin**2/sigma_cluster**2))**2
      elif filterType=='disk':
         result = 1. - np.exp(-0.5*self.RApArcmin**2/sigma_cluster**2)
      elif filterType=='ring':
         result = np.exp(-0.5*self.RApArcmin**2/sigma_cluster**2) - np.exp(-0.5*(self.RApArcmin*np.sqrt(2.))**2/sigma_cluster**2)
      return result



   def ftheoryGaussianProfilePixelated(self, sigma_cluster=1.5, filterType='diskring', dxDeg=0.3, dyDeg= 0.3, resArcmin=0.25, proj='cea', pixwin=0, test=False):
      """Alpha_ksz signal, between 0 and 1.
      Assumes that the projected cluster profile is a 2d Gaussian,
      with sigma_cluster in arcmin
      Assumes equal area disk-ring filter.
      This version is not analytical, but instead generates a mock cutout map
      with a Gaussian profile, and runs the AP filters on it.
      This takes into account the discreteness of the edges of the AP filters
      due to the pixelation of the map
      """

      ###########################################
      # generate pixellized cluster profile map

      # generate null cutout map
      stampMap = self.cutoutGeometry()
      opos = stampMap.posmap()

      # fill the central pixel
      ra = 0.
      dec = 0.
      # coordinates in [rad]
      sourcecoord = np.array([dec, ra]) * np.pi/180.
      # find pixel indices (float) corresponding to ra, dec
      iY, iX = enmap.sky2pix(shape, wcs, sourcecoord, safe=True, corner=False)
      # nearest pixel
      jY = np.int(round(iY))
      jX = np.int(round(iX))
      # fill in the central pixel
      # and normalize to integrate to 1 over angles in [muK*sr]
      pixSizeMap = stampMap.pixsizemap()
      stampMap[jY, jX] = 1. / pixSizeMap[jY, jX] # divide by pixel area in sr

      # convolve map with a Gaussian  profile of given sigma (not fwhm)
      stampMap = enmap.smooth_gauss(stampMap, sigma_cluster * np.pi/180./60.) # convert from arcmin to [rad]

      # apply the pixel window function if desired
      if pixwin<>0:
         stampMap = enmap.apply_window(stampMap, pow=pixwin)


      ###########################################
      # perform the AP filtering

      # create arrays of filter values for the given object
      filtMap = np.zeros(self.nRAp)
  
      # loop over the radii for the AP filter
      for iRAp in range(self.nRAp):
         # Disk radius in rad
         r0 = self.RApArcmin[iRAp] / 60. * np.pi/180.
         # choose an equal area AP filter
         r1 = r0 * np.sqrt(2.)
         
         # perform the filtering
         filtMap[iRAp],_,_,_ = self.aperturePhotometryFilter(opos, stampMap, stampMap, stampMap, r0, r1, filterType=filterType, test=False)
      
      
      if test:
         # compare to the non-pixelated theory profile
         nonPixelated = self.ftheoryGaussianProfile(sigma_cluster, filterType=filterType)
         
         fig=plt.figure(0)
         ax=fig.add_subplot(111)
         #
         ax.plot(self.RApArcmin, nonPixelated, 'k-', label=r'Analytical')
         ax.plot(self.RApArcmin, filtMap, 'b--', label=r'Pixelated')
         #
         ax.legend(loc=4, fontsize='x-small', labelspacing=0.1)
         
         plt.show()
      
      
      return filtMap


   ##################################################################################


   def computeSnrStack(self, filterType, est, tTh='', theory=None, name=None):
      """Compute null rejection, SNR (=detection significance)
      for the requested estimator.
      The estimator considered should have a bootstrap covariance.
      """

      print "- compute SNR and significances for "+filterType+" "+est+" "+tTh

      # replace data with theory if requested
      if tTh=='tsz':
         tTh = '_theory_tsz'
      elif tTh=='ksz':
         tTh = '_theory_ksz'
      else:
         tTh = ''
      
      if name is None:
         name = ''
      else:
         name = '_'+name
      
      if theory is None:
         sigma_cluster = 3. 
         theory = self.ftheoryGaussianProfile(sigma_cluster, filterType=filterType)



      path = self.pathFig+"/snr_"+filterType+"_"+est+tTh+name+".txt"
      with open(path, 'w') as f:
         f.write("*** "+est+" SNR ***\n")

         # data and covariance
         d = self.stackedProfile[filterType+"_"+est+tTh].copy()
         cov = self.covBootstrap[filterType+"_"+est].copy()
         dof = len(d)

         # Compute chi^2_null
         chi2Null = d.dot( np.linalg.inv(cov).dot(d) )
         # goodness of fit for null hypothesis
         f.write("number of dof:"+str(dof)+"\n")
         f.write("null chi2Null="+str(chi2Null)+"\n")
         pteNull = 1.- stats.chi2.cdf(chi2Null, dof)
         f.write("null pte="+str(pteNull)+"\n")
         # pte as a function of sigma, for a Gaussian random variable
         fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.)) - pteNull
         sigmaNull = optimize.brentq(fsigmaToPTE , 0., 1.e3)
         f.write("null pte significance="+str(sigmaNull)+"sigmas\n\n")

         # Gaussian model: find best fit amplitude
         sigma_cluster = 1.5  # arcmin
         def fdchi2(p):
            a = p[0]
            result = (d-a*theory).dot( np.linalg.inv(cov).dot(d-a*theory) )
            result -= chi2Null
            return result
         # Minimize the chi squared
         p0 = 1.
         res = optimize.minimize(fdchi2, p0)
         abest = res.x[0]
         #sbest= res.x[1]
         f.write("best-fit amplitude="+str(abest)+"\n")
         f.write("number of dof:"+str(dof - 1)+"\n\n")

         # goodness of fit for best fit
         chi2Best = fdchi2([abest])+chi2Null
         f.write("best-fit chi2="+str(chi2Best)+"\n")
         pteBest = 1.- stats.chi2.cdf(chi2Best, dof-1.)
         f.write("best-fit pte="+str(pteBest)+"\n")
         # pte as a function of sigma, for a Gaussian random variable
         fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.)) - pteBest
         sigma = optimize.brentq(fsigmaToPTE , 0., 1.e3)
         f.write("best-fit pte significance="+str(sigma)+"sigmas\n\n")

         # favour of best fit over null
         f.write("best-fit sqrt(delta chi2)="+str(np.sqrt(abs(fdchi2([abest]))))+"sigmas\n")
         fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.))
         pte = fsigmaToPTE( np.sqrt(abs(fdchi2([abest]))) )
         f.write("pte (if Gaussian)="+str(pte)+"\n")


   def computeAllSnr(self):
      print "- compute all SNR and significances"
      if self.doBootstrap:
         for filterType in self.filterTypes:
            for est in self.EstBootstrap:
               self.computeSnrStack(filterType, est)
               self.computeSnrStack(filterType, est, tTh='tsz')
               self.computeSnrStack(filterType, est, tTh='ksz')



   ##################################################################################







