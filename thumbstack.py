from headers import *

##################################################################################
##################################################################################

class ThumbStack(object):

   def __init__(self, U, Catalog, pathMap="", pathMask="", pathHit="", name="test", nameLong=None, save=False, nProc=1):
      
      self.nProc = nProc
      self.U = U
      self.Catalog = Catalog
      self.name = name
      if nameLong is None:
         self.nameLong = self.name
      else:
         self.nameLong = nameLong
      self.pathMap = pathMap
      self.pathMask = pathMask
      self.pathHit = pathHit
      
      # Output path
      self.pathOut = "./output/thumbstack/"+self.name
      if not os.path.exists(self.pathOut):
         os.makedirs(self.pathOut)

#      # catalog path
#      self.pathOutCatalog = self.pathOut + "/catalog.txt"

      # Figures path
      self.pathFig = "./figures/catalog/"+self.name
      if not os.path.exists(self.pathFig):
         os.makedirs(self.pathFig)
      
      self.loadMaps()
      self.loadAPRadii()
      
      if save:
         self.createOverlapFlag()
         self.doFiltering()
      
      






   ##################################################################################
   ##################################################################################
   
   def loadAPRadii(self):
   
      # radii to use for AP filter: comoving Mpc/h
      self.nRAp = 15
      self.rApMinMpch = 1. # arcmin
      self.rApMaxMpch = 10. # arcmin
      self.RApMpch = np.linspace(self.rApMinMpch, self.rApMaxMpch, self.nRAp)
   
   
   
   ##################################################################################

   def loadMaps(self):
      print "Read CMB map"
      self.cmbMap = enmap.read_map(self.pathMap)
      print "Read CMB mask"
      self.cmbMask = enmap.read_map(self.pathMask)
      print "Read hit count map"
      self.cmbHit = enmap.read_map(self.pathHit)

      # setup the interpolation algorithm,
      # done once for all, to speed up subsequent calls
      print "Set up map interpolations"
      tStart = time()
      self.cmbMap = utils.interpol_prefilter(self.cmbMap, inplace=True)
      self.cmbMask = utils.interpol_prefilter(self.cmbMask, inplace=True)
      self.cmbHit = utils.interpol_prefilter(self.cmbHit, inplace=True)
      tStop = time()
      print "took", tStop-tStart, "sec"


   ##################################################################################
   
   def createOverlapFlag(self):
   
   
      pass
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   ##################################################################################
   

   def extractStamp(self, ra, dec, dxDeg=1., dyDeg=1., resArcmin=0.25, proj='cea'):
      """Extracts a small CEA or CAR map around the given position, with the given angular size and resolution.
      """

      # define geometry of small square maps to be extracted
      # here 0.5deg * 0.5deg, with 0.25arcmin pixel
      # the center of the small square map will be specified later
      # car: not equal area, but equally-spaced coordinates
      # cea: equal area pixels, but coordinates are not equally spaced
      # previous kSZ paper: went out to 5arcmin
      # Planck paper CV: went out to 18 arcmin
      # --> let's go to 18'*sqrt(2) = 25' \simeq 0.5deg
      # probably good to extract bigger stamp than needed for now
      shape, wcs = enmap.geometry(np.array([[-0.5*dxDeg,-0.5*dyDeg],[0.5*dxDeg,0.5*dyDeg]])*utils.degree, res=resArcmin*utils.arcmin, proj=proj)
      stampMap = enmap.zeros(shape, wcs)

      # coordinates of the square map (between -1 and 1 deg)
      # output map position [{dec,ra},ny,nx]
      opos = stampMap.posmap()

      # coordinate of the center of the square map we want to extract
      # (known point source)
      # ra, dec in this order
      sourcecoord = np.array([ra, dec])*utils.degree

      # corresponding true coordinates on the big healpy map
      ipos = rotfuncs.recenter(opos[::-1], [0,0,sourcecoord[0],sourcecoord[1]])[::-1]

      # extract the small square map by interpolating the big healpy map
      stampMap = self.cmbMap.at(ipos, prefilter=False, mask_nan=False)
      stampMask = self.cmbMask.at(ipos, prefilter=False, mask_nan=False)
      stampHit = self.cmbHit.at(ipos, prefilter=False, mask_nan=False)

      return opos, stampMap, stampMask, stampHit



   ##################################################################################

   def diskRingFilter(self, opos, stampMap, stampMask, stampHit, r0, r1, test=False):
      """Apply an AP filter (disk minus ring) to a stamp map.
      The output is the mean pixel temperature among the pixels in the disk.
      r0 and r1 are the radius of the disk and ring in radians.
      """
      # coordinates of the square map (between -1 and 1 deg, or whatever the size is)
      # output map position [{dec,ra},ny,nx]
#      opos = stampMap.posmap()
      # local coordinates in rad.
      # zero is at the center of the map
      dec = opos[0,:,:]
      ra = opos[1,:,:]
      radius = np.sqrt(ra**2 + dec**2)
      
      # filter
      inDisk = 1.*(radius<=r0)
      inRing = 1.*(radius>r0)*(radius<=r1)
      filter = inDisk / np.sum(inDisk) - inRing / np.sum(inRing)
      # count nb of pixels where filter is strictly positive
      nbPix = len(np.where(filter>0.)[0])
      # estimate area of strictly positive part of filter
      pixArea = ra.area() / len(ra.flatten()) # in sr
      diskArea = np.sum(inDisk) * pixArea  # disk area in sr
      
      # apply the filter to the maps
      filtMap = np.sum(filter * stampMap)
      filtMask = np.sum(filter * stampMask)
      filtVar = np.sum(filter**2 / stampHit) # to get the variance (arbitrary units)

      if test:
         print "- nb of pixels where filter>0: "+str(nbPix)
         print "  ie area of "+str(diskArea)+" sr"
         
         print "- disk-ring filter sums over pixels to "+str(np.sum(filter))
         print "  (should be 0; to be compared with "+str(len(filter.flatten()))+")"

         print "- filter on map:"+str(filtMap)+" muK"
         print "- filter on mask:"+str(filtMask)+" muK"
         print "- filter on inverse hit:"+str(filtHit)+" muK"

      return filtMap, filtMask, filtVar, diskArea


   ##################################################################################


   def doFiltering(self):
      
      # initialize arrays
      self.filtMap = np.zeros((self.Catalog.nObj, self.nRAp))
      self.filtMask = np.zeros((self.Catalog.nObj, self.nRAp))
      self.filtVar = np.zeros((self.Catalog.nObj, self.nRAp))
      self.diskArea = np.zeros((self.Catalog.nObj, self.nRAp))
      
      
      # analysis to be done for each object
      def analyzeObject(iObj):
         print "- analyzing object" + str(iObj)
         
         # Object coordinates
         ra = self.Catalog.RA[iObj]   # in deg
         dec = self.Catalog.DEC[iObj] # in deg
         z = self.Catalog.Z[iObj]
         
         # extract postage stamp around it
         opos, stampMap, stampMask, stampHit = self.extractStamp(ra, dec, dxDeg=1., dyDeg=1., resArcmin=0.25, proj='cea')
         
         # create arrays of filter values for the given object
         filtMap = np.zeros(self.nRAp)
         filtMask = np.zeros(self.nRAp)
         filtVar = np.zeros(self.nRAp)
         diskArea = np.zeros(self.nRAp)
         
         # loop over the radii for the AP filter
         for iRAp in range(self.nRAp):
            # disk radius in comoving Mpc/h
            rApMpch = self.RApMpch[iRAp]
            # convert to radians at the given redshift
            r0 = rApMpch / self.U.bg.comoving_transverse_distance(z) # rad
            # choose an equal area AP filter
            r1 = r0 * np.sqrt(2.)
            # perform the filtering
            filtMap[iRAp], filtMask[iRAp], filtVar[iRAp], diskArea[iRAp] = self.diskRingFilter(opos, stampMap, stampMask, stampHit, r0, r1, test=False)
         return filtMap, filtMask, filtVar, diskArea


      # loop over objects in the catalog
#      pool = Pool(self.nProc)
#      result = np.array(pool.map(analyzeObject, range(self.Catalog.nObj)))
      result = np.array(map(analyzeObject, range(self.Catalog.nObj)))
      self.filtMap = result[:,0,:].copy()
      self.filtMask = result[:,1,:]
      self.filtVar = result[:,2,:]
      self.diskArea = result[:,3,:]


