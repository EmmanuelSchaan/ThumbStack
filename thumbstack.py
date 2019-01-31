from headers import *

##################################################################################
##################################################################################

class ThumbStack(object):

   def __init__(self, U, Catalog, pathMap="", pathMask="", pathHit="", name="test", nameLong=None, save=False):

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
      
      
      if save:
         pass
      
      self.loadMaps()






   ##################################################################################
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
      print "Set up map interpolation"
      tStart = time()
      self.cmbMap = utils.interpol_prefilter(self.cmbMap, inplace=True)
      tStop = time()
      print "took", tStop-tStart, "sec"


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

      return stampMap, stampMask, stampHit



   ##################################################################################

   def filterStamp(self, stampMap):
      pass
