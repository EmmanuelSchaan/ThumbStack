from headers import *

##################################################################################
##################################################################################

class ThumbStack(object):

#   def __init__(self, U, Catalog, pathMap="", pathMask="", pathHit="", name="test", nameLong=None, save=False, nProc=1):
   def __init__(self, U, Catalog, cmbMap, cmbMask, cmbHit, name="test", nameLong=None, save=False, nProc=1):
      
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
      
      # number of samples for bootstraps, shuffles
      self.nSamples = 1000
      
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
      
      if save:
         self.saveOverlapFlag(nProc=self.nProc)
      self.loadOverlapFlag()
      
      if save:
         self.saveFiltering(nProc=self.nProc)
      self.loadFiltering()
      
      self.measureVarFromHitCount(plot=save)

      if save:
         self.saveStackedProfiles()
      self.loadStackedProfiles()

      if False:
         self.plotAllStackedProfiles()
         self.plotAllCov()
         self.computeAllSnr()



   ##################################################################################
   ##################################################################################
   
   def loadAPRadii(self):
      # radii to use for AP filter: comoving Mpc/h
      self.nRAp = 9  #4
      
      # Aperture radii in Mpc/h
      self.rApMinMpch = 1.
      self.rApMaxMpch = 5
      self.RApMpch = np.linspace(self.rApMinMpch, self.rApMaxMpch, self.nRAp)
      
      # Aperture radii in arcmin
      self.rApMinArcmin = 1.  # 1.
      self.rApMaxArcmin = 6.  # 4.
      self.RApArcmin = np.linspace(self.rApMinArcmin, self.rApMaxArcmin, self.nRAp)


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
   
   
   
   def saveOverlapFlag(self, thresh=1.e-5, nProc=1):
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
         hit = self.sky2map(ra, dec, self.cmbHit)
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
   

   def extractStamp(self, ra, dec, dxDeg=0.3, dyDeg=0.3, resArcmin=0.25, proj='cea', test=False):
      """Extracts a small CEA or CAR map around the given position, with the given angular size and resolution.
      ra, dec in degrees.
      Does it for the map, the mask and the hit count.
      """

      # Make sure the number of pixels is (2*n+1),
      # so that the object is exactly in the middle of the central pixel
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


   def diskRingFilter(self, opos, stampMap, stampMask, stampHit, r0, r1, test=False):
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
      filtNoiseStdDev: [1/sqrt(hit unit) * sr], ie [std dev * sr] if [hit map] = inverse var
      diskArea: [sr]
      """
      # coordinates of the square map in radians
      # zero is at the center of the map
      # output map position [{dec,ra},ny,nx]
      dec = opos[0,:,:]
      ra = opos[1,:,:]
      radius = np.sqrt(ra**2 + dec**2)
      
      # disk filter [dimensionless]
      inDisk = 1.*(radius<=r0)
      # ring filter [dimensionless], normalized so that the disk-ring filter integrates exactly to zero
      inRing = 1.*(radius>r0)*(radius<=r1)
      inRing *= np.sum(inDisk) / np.sum(inRing)
      # disk minus ring filter [dimensionless]
      filter = inDisk - inRing
      
      # exact angular area of pixel [sr]
      pixArea = ra.area() / len(ra.flatten())
      # exact angular area of disk [sr]
      diskArea = np.sum(inDisk) * pixArea

      # apply the filter: int_disk d^2theta map -  disk_area / ring_area * int_ring d^2theta map
      filtMap = np.sum(pixArea * filter * stampMap)   # [map unit * sr]
      
      # detect point sources within the filter:
      # gives 0 in the absence of point sources/edges; gives >=1 in the presence of point sources/edges
      filtMask = np.sum((radius<=r1) * (1-stampMask))   # [dimensionless]
      
      # quantify noise std dev in the filter
      #!!! hardcoding a small number to avoid dividing by zero. Ugly.
      filtNoiseStdDev = np.sqrt(np.sum((pixArea * filter)**2 / (1.e-16 + stampHit))) # to get the std devs [sr / sqrt(hit unit)]
      
      #print "filtNoiseStdDev = ", filtNoiseStdDev
      if np.isnan(filtNoiseStdDev):
         print "filtNoiseStdDev is nan"
         print stampHit


      if test:
         print "AP filter with disk radius =", r0 * (180.*60./np.pi), "arcmin"

         # count nb of pixels where filter is strictly positive
         nbPix = len(np.where(filter>0.)[0])
         print "- nb of pixels in the cutout: "+str(filter.shape[0] * filter.shape[1])
         print "= nb of pixels where filter=0: "+str(len(np.where(filter==0.)[0]))
         print "+ nb of pixels where filter>0: "+str(len(np.where(filter>0.)[0]))
         print "+ nb of pixels where filter<0: "+str(len(np.where(filter<0.)[0]))
         print "- disk area: "+str(diskArea)+" sr, ie "+str(diskArea * (180.*60./np.pi)**2)+"arcmin^2"
         print "  (from r0, expect "+str(np.pi*r0**2)+" sr, ie "+str(np.pi*r0**2 * (180.*60./np.pi)**2)+"arcmin^2)"
         
         print "- disk-ring filter sums over pixels to "+str(np.sum(filter))
         print "  (should be 0; compared to "+str(len(filter.flatten()))+")"
         
         print "- filter on unit disk: "+str(np.sum(pixArea * filter * inDisk))
         print "  (should be disk area in sr: "+str(diskArea)+")"
         print "- filter on map: "+str(filtMap)
         print "- filter on mask: "+str(filtMask)
         print "- filter on inverse hit: "+str(filtNoiseStdDev)

         print "- plot the filter"
         filterMap = stampMap.copy()
         filterMap[:,:] = filter.copy()
         plots=enplot.plot(filterMap,grid=True)
         enplot.write(self.pathTestFig+"/stampfilter_r0"+floatExpForm(r0)+"_r1"+floatExpForm(r1), plots)

      return filtMap, filtMask, filtNoiseStdDev, diskArea





   ##################################################################################



   def analyzeObject(self, iObj, test=False):
      '''Analysis to be done for each object.
      Returns:
      filtMap: [map unit * sr]
      filtMask: [mask unit * sr]
      filtNoiseStdDev: [1/sqrt(hit unit) * sr], ie [std dev * sr] if [hit map] = inverse var
      diskArea: [sr]
      '''
      
      if iObj%1000==0:
         print "- analyze object", iObj
      
      # create arrays of filter values for the given object
      filtMap = np.zeros(self.nRAp)
      filtMask = np.zeros(self.nRAp)
      filtNoiseStdDev = np.zeros(self.nRAp)
      diskArea = np.zeros(self.nRAp)
      
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
         opos, stampMap, stampMask, stampHit = self.extractStamp(ra, dec, dxDeg=dDeg, dyDeg=dDeg, resArcmin=0.25, proj='cea', test=test)
         
         # loop over the radii for the AP filter
         for iRAp in range(self.nRAp):
#            # disk radius in comoving Mpc/h
#            rApMpch = self.RApMpch[iRAp]
#            # convert to radians at the given redshift
#            r0 = rApMpch / self.U.bg.comoving_transverse_distance(z) # rad

            # Disk radius in rad
            r0 = self.RApArcmin[iRAp] / 60. * np.pi/180.
            # choose an equal area AP filter
            r1 = r0 * np.sqrt(2.)
            
            # perform the filtering
            filtMap[iRAp], filtMask[iRAp], filtNoiseStdDev[iRAp], diskArea[iRAp] = self.diskRingFilter(opos, stampMap, stampMask, stampHit, r0, r1, test=test)

#      if test:
#         print " plot the measured profile"
#         fig=plt.figure(0)
#         ax=fig.add_subplot(111)
#         #
#         ax.plot(self.r, filtMap)

      return filtMap, filtMask, filtNoiseStdDev, diskArea



   def saveFiltering(self, nProc=1):
      
      # initialize arrays
      self.filtMap = np.zeros((self.Catalog.nObj, self.nRAp))
      self.filtMask = np.zeros((self.Catalog.nObj, self.nRAp))
      self.filtNoiseStdDev = np.zeros((self.Catalog.nObj, self.nRAp))
      self.diskArea = np.zeros((self.Catalog.nObj, self.nRAp))

      # loop over all objects in catalog
#      result = np.array(map(self.analyzeObject, range(self.Catalog.nObj)))
      tStart = time()
      with sharedmem.MapReduce(np=nProc) as pool:
         result = np.array(pool.map(self.analyzeObject, range(self.Catalog.nObj)))
      tStop = time()
      print "took", (tStop-tStart)/60., "min"
      
      # unpack and save to file
      self.filtMap = result[:,0,:].copy()
      np.savetxt(self.pathOut+"/filtmap.txt", self.filtMap)
      #
      self.filtMask = result[:,1,:].copy()
      np.savetxt(self.pathOut+"/filtmask.txt", self.filtMask)
      #
      self.filtNoiseStdDev = result[:,2,:].copy()
      np.savetxt(self.pathOut+"/filtnoisestddev.txt", self.filtNoiseStdDev)
      #
      self.diskArea = result[:,3,:].copy()
      np.savetxt(self.pathOut+"/diskarea.txt", self.diskArea)


   def loadFiltering(self):
      self.filtMap = np.genfromtxt(self.pathOut+"/filtmap.txt")
      self.filtMask = np.genfromtxt(self.pathOut+"/filtmask.txt")
      self.filtNoiseStdDev = np.genfromtxt(self.pathOut+"/filtnoisestddev.txt")
      self.diskArea = np.genfromtxt(self.pathOut+"/diskarea.txt")


   ##################################################################################


   def catalogMask(self, overlap=True, psMask=True, mVir=[1.e6, 1.e17], extraSelection=1.):
      '''Returns catalog mask: 1 for objects to keep, 0 for objects to discard.
      Use as:
      maskedQuantity = Quantity[mask]
      '''
      # Here mask is 1 for objects we want to keep
      mask = np.ones_like(self.Catalog.RA)
      #print "keeping fraction", np.sum(mask)/len(mask), " of objects"
      if mVir is not None:
         mask *= (self.Catalog.Mvir>=mVir[0]) * (self.Catalog.Mvir<=mVir[1])
         #print "keeping fraction", np.sum(mask)/len(mask), " of objects"
      if overlap:
         mask *= self.overlapFlag.copy()
         #print "keeping fraction", np.sum(mask)/len(mask), " of objects"
      # PS mask: look at largest aperture, and remove if any point within the disk or ring is masked
      if psMask:
         mask *= 1.*(np.abs(self.filtMask[:,-1])<1.)
         #print "keeping fraction", np.sum(mask)/len(mask), " of objects"
      mask *= extraSelection
      #print "keeping fraction", np.sum(mask)/len(mask), " of objects"

      mask = mask.astype(bool)
      #print "keeping fraction", np.sum(mask)/len(mask), " of objects"
      return mask


   ##################################################################################


   def measureVarFromHitCount(self, plot=False):
      """Returns a list of functions, one for each AP filter radius,
      where the function takes filtNoiseStdDev**2 \propto [(map var) * sr^2] as input and returns the
      actual measured filter variance [(map unit)^2 * sr^2].
      The functions are expected to be linear if the detector noise is the main source of noise,
      and if the hit counts indeed reflect the detector noise.
      To be used for noise weighting in the stacking.
      """
      print "- interpolate the relation hit count - noise"
      # keep only objects that overlap, and mask point sources
      mask = self.catalogMask(overlap=True, psMask=True)
   
      self.fVarFromHitCount = []
      for iRAp in range(self.nRAp):
         x = self.filtNoiseStdDev[mask, iRAp]**2
         y = self.filtMap[mask, iRAp].copy()
         y = (y - np.mean(y))**2

#         print "values in x"
#         print np.mean(x), np.std(x), np.max(x), np.min(x)
#         print "values in y"
#         print np.mean(y), np.std(y), np.max(y), np.min(y)

         
         # Check whether the hit count actually varies appreciably,
         # otherwise use uniform weighting
         if np.std(x) < 0.001 * np.mean(x):
            print "Using uniform weighting: hit count does not vary much."
            self.fVarFromHitCount.append(lambda x: np.mean(y)*np.ones_like(x))
         else:
            print "Using hit count weighting: hit count varies much."
            # define bins of hit count values
            nBins = 21
            BinsX = np.logspace(np.log10(np.min(x)), np.log10(np.max(x)), nBins, 10.)
            
            # compute histograms
            binCenters, binEdges, binIndices = stats.binned_statistic(x, x, statistic='mean', bins=BinsX)
            binCounts, binEdges, binIndices = stats.binned_statistic(x, x, statistic='count', bins=BinsX)
            binnedVar, binEdges, binIndices = stats.binned_statistic(x, y, statistic=np.mean, bins=BinsX)
            sBinnedVar, binEdges, binIndices = stats.binned_statistic(x, y, statistic=np.std, bins=BinsX)
            sBinnedVar /= np.sqrt(binCounts)
            
            # interpolate, to use as noise weighting
            self.fVarFromHitCount.append( interp1d(binCenters, binnedVar, kind='linear', bounds_error=False, fill_value=(binnedVar[0],binnedVar[-1])) )
            
            if plot:
               # plot
               fig=plt.figure(0)
               ax=fig.add_subplot(111)
               #
               # measured
               ax.errorbar(binCenters, binnedVar*(180.*60./np.pi)**2, yerr=sBinnedVar*(180.*60./np.pi)**2, fmt='.', label=r'measured')
               #
               # interpolated
               newX = np.logspace(np.log10(np.min(x)/2.), np.log10(np.max(x)*2.), 10.*nBins, 10.)
               newY = np.array(map(self.fVarFromHitCount[iRAp], newX))
               ax.plot(newX, newY*(180.*60./np.pi)**2, label=r'interpolated')
               #
               ax.set_xscale('log', nonposx='clip')
               ax.set_yscale('log', nonposy='clip')
               ax.set_xlabel(r'Det. noise var. from combined hit [arbitrary]')
               ax.set_ylabel(r'Measured var. [$\mu$K.arcmin$^2$]')
               #
               path = self.pathFig+"/binned_noise_vs_hit"+str(iRAp)+".pdf"
               fig.savefig(path, bbox_inches='tight')
               fig.clf()

      return


   ##################################################################################


   def stackedProfile(self, est, iBootstrap=None, iVShuffle=None, tTh=None):
      """Returns the estimated profile and its uncertainty for each aperture.
      est: string to select the estimator
      iBootstrap: index for bootstrap resampling
      iVShuffle: index for shuffling velocities
      tTh: to replace measured temperatures by a theory expectation
      """
      # select objects that overlap, and reject point sources
      mask = self.catalogMask(overlap=True, psMask=True)
      
      # temperatures [muK * sr]
      if tTh is None:
         t = self.filtMap.copy() # [muK * sr]
      elif tTh=='tsz':
         # expected tSZ signal
         # AP profile shape, between 0 and 1
         sigma_cluster = 1.5  # arcmin
         shape = self.ftheoryGaussianProfile(sigma_cluster) # between 0 and 1 [dimless]
         # multiply by integrated y to get y profile [sr]
         t = np.column_stack([self.Catalog.integratedY[:] * shape[iAp] for iAp in range(self.nRAp)])
         # convert from y profile to dT profile
         Tcmb = 2.726   # K
         h = 6.63e-34   # SI
         kB = 1.38e-23  # SI
         def f(nu):
            """frequency dependence for tSZ temperature
            """
            x = h*nu/(kB*Tcmb)
            return x*(np.exp(x)+1.)/(np.exp(x)-1.) -4.
         t *= 2. * f(150.e9) * Tcmb * 1.e6  # [muK * sr]
      elif tTh=='ksz':
         # expected kSZ signal
         # AP profile shape, between 0 and 1
         sigma_cluster = 1.5  # arcmin
         shape = self.ftheoryGaussianProfile(sigma_cluster) # between 0 and 1 [dimless]
         # multiply by integrated kSZ to get kSZ profile [muK * sr]
         t = np.column_stack([self.Catalog.integratedKSZ[:] * shape[iAp] for iAp in range(self.nRAp)])   # [muK * sr]
      t = t[mask, :]
      tNoMean = t - np.mean(t, axis=0)

      # v/c [dimless]
      v = -self.Catalog.vR[mask] / 3.e5
      v -= np.mean(v)
      # hit count and measured total noise (CMB + detector)
      s2Hit = self.filtNoiseStdDev[mask, :]**2
      #print "Shape of s2Hit = ", s2Hit.shape

#      s2Full = self.fVarFromHitCount[:](s2Hit)
      s2Full = np.column_stack([self.fVarFromHitCount[iAp](s2Hit[:, iAp]) for iAp in range(self.nRAp)])
      #print "Shape of s2Full = ", s2Full.shape

      # halo masses
      m = self.Catalog.Mvir[mask]
      
      if iBootstrap is not None:
         # make sure each resample is independent,
         # and make the resampling reproducible
         np.random.seed(iBootstrap)
         # list of overlapping objects
         nObj = np.sum(mask)
         I = np.arange(nObj)
         # choose with replacement from this list
         J = np.random.choice(I, size=nObj, replace=True)
         #
         t = t[J,:]
         tNoMean = tNoMean[J,:]
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
         stack = np.mean(t, axis=0)
         sStack = np.sqrt(np.sum(s2Full, axis=0)) / t.shape[0]
      # tSZ: detector-noise weighted (hit count)
      elif est=='tsz_hitweight':
         stack = np.sum(t/s2Hit, axis=0) / np.sum(1./s2Hit, axis=0)
         sStack = np.sum(s2Full/s2Hit**2, axis=0) / np.sum(1./s2Hit, axis=0)**2
         sStack = np.sqrt(sStack)
      # tSZ: full noise weighted (detector noise + CMB)
      elif est=='tsz_varweight':
         stack = np.sum(t/s2Full , axis=0) / np.sum(1./s2Full, axis=0)
         sStack = 1. / np.sum(1./s2Full, axis=0)
         sStack = np.sqrt(sStack)

      # kSZ: uniform weighting
      elif est=='ksz_uniformweight':
         stack = np.sum(tNoMean*v[:,np.newaxis], axis=0) / np.sum(v[:,np.newaxis]**2, axis=0)
         sStack = np.sum(s2Full*v[:,np.newaxis]**2, axis=0) / np.sum(v[:,np.newaxis]**2, axis=0)**2
         sStack = np.sqrt(sStack)
         # renormalize the kSZ estimator
         stack *= np.std(v)
         sStack *= np.std(v)
      # kSZ: detector-noise weighted (hit count)
      elif est=='ksz_hitweight':
         stack = np.sum(tNoMean*v[:,np.newaxis]/s2Hit, axis=0) / np.sum(v[:,np.newaxis]**2/s2Hit, axis=0)
         sStack = np.sum(s2Full*v[:,np.newaxis]**2/s2Hit**2, axis=0) / np.sum(v[:,np.newaxis]**2/s2Hit, axis=0)**2
         sStack = np.sqrt(sStack)
         # renormalize the kSZ estimator
         stack *= np.std(v)
         sStack *= np.std(v)
      # kSZ: full noise weighted (detector noise + CMB)
      elif est=='ksz_varweight':
         stack = np.sum(tNoMean*v[:,np.newaxis]/s2Full, axis=0) / np.sum(v[:,np.newaxis]**2/s2Full, axis=0)
         sStack = np.sum(s2Full*v[:,np.newaxis]**2/s2Full**2, axis=0) / np.sum(v[:,np.newaxis]**2/s2Full, axis=0)**2
         sStack = np.sqrt(sStack)
         # renormalize the kSZ estimator
         stack *= np.std(v)
         sStack *= np.std(v)
      # kSZ: full noise weighted (detector noise + CMB)
      elif est=='ksz_massvarweight':
         stack = np.sum(tNoMean*m[:,np.newaxis] * v[:,np.newaxis]/s2Full, axis=0)
         stack /= np.sum(m[:,np.newaxis]**2 * v[:,np.newaxis]**2/s2Full, axis=0)
         sStack = np.sum(s2Full * m[:,np.newaxis]**2 * v[:,np.newaxis]**2/s2Full**2, axis=0)
         sStack /= np.sum(m[:,np.newaxis]**2 * v[:,np.newaxis]**2/s2Full, axis=0)**2
         sStack = np.sqrt(sStack)
         # renormalize the kSZ estimator
         stack *= np.mean(m)
         sStack *= np.mean(m)
         stack *= np.std(v)
         sStack *= np.std(v)

      return stack, sStack


   ##################################################################################

   def SaveCovBootstrapStackedProfile(self, est, nSamples=1000, nProc=1):
      """Estimate covariance matrix for the stacked profile from bootstrap resampling
      """
      tStart = time()
      with sharedmem.MapReduce(np=nProc) as pool:
         f = lambda iSample: self.stackedProfile(est, iBootstrap=iSample)
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
      np.savetxt(self.pathOut+"/cov_"+est+"_bootstrap.txt", covStack)
      

   def SaveCovVShuffleStackedProfile(self, est, nSamples=1000, nProc=1):
      """Estimate covariance matrix for the stacked profile from shuffling velocities
      """
      tStart = time()
      with sharedmem.MapReduce(np=nProc) as pool:
         f = lambda iSample: self.stackedProfile(est, iVShuffle=iSample)
         result = np.array(pool.map(f, range(nSamples)))
      tStop = time()
      print "took", (tStop-tStart)/60., "min"
      # unpack results
      stackSamples = result[:,0,:] # shape (nSamples, nRAp)
      #sStackSamples = result[:,1,:]
      # estimate cov
      covStack = np.cov(stackSamples, rowvar=False)
      # save it to file
      np.savetxt(self.pathOut+"/cov_"+est+"_vshuffle.txt", covStack)


   ##################################################################################

   def saveStackedProfiles(self):
      print "- compute stacked profiles and their cov"
      data = np.zeros((self.nRAp, 3))
      data[:,0] = self.RApArcmin # [arcmin]
      
      # all stacked profiles
      Est = ['tsz_uniformweight', 'tsz_hitweight', 'tsz_varweight', 'ksz_uniformweight', 'ksz_hitweight', 'ksz_varweight', 'ksz_massvarweight']
      for iEst in range(len(Est)):
         est = Est[iEst]
         # measured stacked profile
         data[:,1], data[:,2] = self.stackedProfile(est) # [map unit * sr]
         np.savetxt(self.pathOut+"/"+est+"_measured.txt", data)
         # expected stacked profile from tSZ
         data[:,1], data[:,2] = self.stackedProfile(est, tTh='tsz') # [map unit * sr]
         np.savetxt(self.pathOut+"/"+est+"_theory_tsz.txt", data)
         # expected stacked profile from kSZ
         data[:,1], data[:,2] = self.stackedProfile(est, tTh='ksz') # [map unit * sr]
         np.savetxt(self.pathOut+"/"+est+"_theory_ksz.txt", data)

      # covariance matrices from bootstrap,
      # only for a few select estimators
      Est = ['tsz_varweight', 'ksz_varweight']
      for iEst in range(len(Est)):
         est = Est[iEst]
         self.SaveCovBootstrapStackedProfile(est, nSamples=1000, nProc=self.nProc)

      # covariance matrices from shuffling velocities,
      # for ksz only
      Est = ['ksz_varweight']
      for iEst in range(len(Est)):
         est = Est[iEst]
         self.SaveCovVShuffleStackedProfile(est, nSamples=1000, nProc=self.nProc)


   def loadStackedProfiles(self):
      print "- load stacked profiles and their cov"
      self.stackedProfile = {}
      self.sStackedProfile = {}
      self.covBootstrap = {}
      self.covVShuffle = {}
      
      # all stacked profiles
      Est = ['tsz_uniformweight', 'tsz_hitweight', 'tsz_varweight', 'ksz_uniformweight', 'ksz_hitweight', 'ksz_varweight', 'ksz_massvarweight']
      for iEst in range(len(Est)):
         est = Est[iEst]
         # measured stacked profile
         data = np.genfromtxt(self.pathOut+"/"+est+"_measured.txt")
         self.stackedProfile[est] = data[:,1]
         self.sStackedProfile[est] = data[:,2]
         # expected stacked profile from tSZ
         data = np.genfromtxt(self.pathOut+"/"+est+"_theory_tsz.txt")
         self.stackedProfile[est+"_theory_tsz"] = data[:,1]
         self.sStackedProfile[est+"_theory_tsz"] = data[:,2]
         # expected stacked profile from kSZ
         data = np.genfromtxt(self.pathOut+"/"+est+"_theory_ksz.txt")
         self.stackedProfile[est+"_theory_ksz"] = data[:,1]
         self.sStackedProfile[est+"_theory_ksz"] = data[:,2]

      # covariance matrices from bootstrap,
      # only for a few select estimators
      Est = ['tsz_varweight', 'ksz_varweight']
      for iEst in range(len(Est)):
         est = Est[iEst]
         self.covBootstrap[est] = np.genfromtxt(self.pathOut+"/cov_"+est+"_bootstrap.txt")
      
      # covariance matrices from shuffling velocities,
      # for ksz only
      Est = ['ksz_varweight']
      for iEst in range(len(Est)):
         est = Est[iEst]
         self.covVShuffle[est] = np.genfromtxt(self.pathOut+"/cov_"+est+"_vshuffle.txt")


   ##################################################################################

   def plotStackedProfile(self, Est, name=None):
      """Compares stacked profiles, and their uncertainties.
      """
      if name is None:
         name = Est[0]
      
      # stacked profile
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # convert from sr to arcmin^2
      factor = (180.*60./np.pi)**2
      #
      ax.axhline(0., c='k', lw=1)
      #
      colors = ['r', 'g', 'b', 'm', 'c']
      for iEst in range(len(Est)):
         est = Est[iEst]
         c = colors[iEst]
         ax.errorbar(self.RApArcmin, factor * self.stackedProfile[est], factor * self.sStackedProfile[est], fmt='-', c=c, label="measured")
         ax.plot(self.RApArcmin, factor * self.stackedProfile[est+"_theory_tsz"], ls='--', c=c, label="theory tsz")
         ax.plot(self.RApArcmin, factor * self.stackedProfile[est+"_theory_ksz"], ls='-.', c=c, label="theory ksz")
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'$T$ [$\mu K\cdot\text{arcmin}^2$]')
      #ax.set_ylim((0., 2.))
      #
      path = self.pathFig+"/"+name+".pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()

      # uncertainty on stacked profile
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # convert from sr to arcmin^2
      factor = (180.*60./np.pi)**2
      #
      for est in Est:
         ax.plot(self.RApArcmin, factor * self.sStackedProfile[est], ls='-', label="analytic")
         if est in self.covBootstrap:
            ax.plot(self.RApArcmin, factor * np.sqrt(np.diag(self.covBootstrap[est])), ls='--', label="bootstrap")
         if est in self.covVShuffle:
            ax.plot(self.RApArcmin, factor * np.sqrt(np.diag(self.covVShuffle[est])), ls='-.', label="v shuffle")
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'$\sigma(T)$ [$\mu K\cdot\text{arcmin}^2$]')
      #ax.set_ylim((0., 2.))
      #
      path = self.pathFig+"/s_"+name+".pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()

   def plotAllStackedProfiles(self):
      print "- plot all stacked profiles"
#      # tSZ
#      Est = ['tsz_uniformweight', 'tsz_hitweight', 'tsz_varweight']
#      self.plotStackedProfile(Est, name="tsz")
#      # kSZ
#      Est = ['ksz_uniformweight', 'ksz_hitweight', 'ksz_varweight', 'ksz_massvarweight']
#      self.plotStackedProfile(Est, name="ksz")
      
      # all stacked profiles
      Est = ['tsz_uniformweight', 'tsz_hitweight', 'tsz_varweight', 'ksz_uniformweight', 'ksz_hitweight', 'ksz_varweight', 'ksz_massvarweight']
      for iEst in range(len(Est)):
         est = Est[iEst]
         self.plotStackedProfile([est], name=est)

   ##################################################################################
   
   def plotCov(self, cov, name=""):
      # Covariance matrix
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # pcolor wants x and y to be edges of cell,
      # ie one more element, and offset by half a cell
      dR = (self.rApMaxArcmin - self.rApMinArcmin) / self.nRAp
      RApEdgesArcmin = np.linspace(self.rApMinArcmin-0.5*dR, self.rApMaxArcmin+0.5*dR, self.nRAp+1)
      #
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
      fig.clf()

      # Correlation coefficient
      fig=plt.figure(0)
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
      fig.clf()

   def plotAllCov(self):
      print "- plot all covariances"
      # covariance matrices from bootstrap,
      # only for a few select estimators
      Est = ['tsz_varweight', 'ksz_varweight']
      for iEst in range(len(Est)):
         est = Est[iEst]
         self.plotCov(self.covBootstrap[est], est+"_bootstrap")
      
      # covariance matrices from shuffling velocities,
      # for ksz only
      Est = ['ksz_varweight']
      for iEst in range(len(Est)):
         est = Est[iEst]
         self.plotCov(self.covVShuffle[est], est+"_vshuffle")


   ##################################################################################

   def ftheoryGaussianProfile(self, sigma_cluster):
      """Alpha_ksz signal, between 0 and 1.
      Assumes that the projected cluster profile is a 2d Gaussian,
      with sigma_cluster in arcmin
      Assumes equal area disk-ring filter
      """
      return (1. - np.exp(-0.5*self.RApArcmin**2/sigma_cluster**2))**2


   def ftheoryGaussianProfilePixelated(self, sigma_cluster=1.5, dxDeg=0.3, dyDeg= 0.3, resArcmin=0.25, proj='cea', pixwin=0, test=False):
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

         # Make sure the number of pixels is (2*n+1),
         # so that the object is exactly in the middle of the central pixel
         nx = np.floor((dxDeg * 60. / resArcmin - 1.) / 2.) + 1.
         dxDeg = (2. * nx + 1.) * resArcmin / 60.
         ny = np.floor((dyDeg * 60. / resArcmin - 1.) / 2.) + 1.
         dyDeg = (2. * ny + 1.) * resArcmin / 60.

         # generate null map
         shape, wcs = enmap.geometry(np.array([[-0.5*dxDeg,-0.5*dyDeg],[0.5*dxDeg,0.5*dyDeg]])*utils.degree, res=resArcmin*utils.arcmin, proj=proj)
         stampMap = enmap.zeros(shape, wcs)
         #
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
            filtMap[iRAp],_,_,_ = self.diskRingFilter(opos, stampMap, stampMap, stampMap, r0, r1, test=False)
         
         
         if test:
            # compare to the non-pixelated theory profile
            nonPixelated = self.ftheoryGaussianProfile(sigma_cluster)
            
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


   def computeSnrStack(self, est):
      """Compute null rejection, SNR (=detection significance)
      for the requested estimator.
      The estimator considered should have a bootstrap covariance.
      """

      path = self.pathFig+"/snr_"+est+".txt"
      with open(path, 'w') as f:
         f.write("*** "+est+" SNR ***\n")

         # data and covariance
         d = self.stackedProfile[est].copy()
         cov = self.covBootstrap[est].copy()
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
         theory = self.ftheoryGaussianProfile(sigma_cluster)
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
      Est = ['tsz_varweight', 'ksz_varweight']
      for iEst in range(len(Est)):
         est = Est[iEst]
         self.computeSnrStack(est)



   ##################################################################################





























##      myHistogram(DEC[mask], nBins=71, lim=(-90., 90.), path=path, nameLatex=r'$x$ [km/s]', semilogx=False, doGauss=False)
#
#      # first, keep all objects that overlap, even the masked ones
#      mask = self.catalogMask(overlap=True, psMask=False)
#
#      # check that the non-overlapping objects are the ones with the correct DEC
#      path = self.pathFig+"/hist_dec_overlap.pdf"
#      myHistogram(self.Catalog.DEC[mask], nBins=71, lim=(-30., 90.), path=path, nameLatex=r'Dec [deg]')
#
#      # check the values of the filters on the point source mask, to find a relevant cut
#      # look at the largest aperture
#      x = self.filtMask[mask,-1]
#      path = self.pathFig+"/hist_psmaskvalue_before.pdf"
#      myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'PS mask value', semilogy=True)
#
#      # then  remove the objects that overlap with point sources
#      mask = self.catalogMask(overlap=True, psMask=True)
#
#      # redo the mask histogram, to check
#      x = self.filtMask[mask,-1]
#      path = self.pathFig+"/hist_psmaskvalue_after.pdf"
#      myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'PS mask value', semilogy=True)
#
#      # Histograms of filter outputs
#      for iRAp in range(self.nRAp):
#         x = self.filtMap[mask, iRAp]
#         path = self.pathFig+"/hist_filtvalue"+str(iRAp)+".pdf"
#         S2Theory = []
#         for f in fsAp:
#            # theory assumes int d^2theta W = 1
#            s2Theory = f(self.RApArcmin[iRAp] * np.pi/180./60.)**2
#            # so multiply by disk area, but it can vary from object to object
#            # we neglect this source of scatter and just keep the mean
#            s2Theory *= np.mean(self.diskArea[mask, iRAp])**2
#            S2Theory.append(s2Theory)
#         myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'AP filter value', S2Theory=S2Theory, doGauss=True, semilogy=True)
#
#
#      # Histograms of noise std dev, from hit counts
#      for iRAp in range(self.nRAp):
#         x = self.filtNoiseStdDev[mask, iRAp]
#         path = self.pathFig+"/hist_noisestddevhit"+str(iRAp)+".pdf"
#         myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'Std dev value [arbitrary]', semilogy=True)
#
#
#      # tSZ / dust
#      for iRAp in range(self.nRAp):
#         weights = 1. / self.filtNoiseStdDev[mask, iRAp]**2   # need inverse variance, not std dev
#         weights /= np.mean(weights)   # to keep the size and units of the weighted AP outputs
#         x = self.filtMap[mask, iRAp] * weights
#         print "- mean tSZ= "+str(np.mean(x))+"; std on mean= "+str(np.std(x)/np.sqrt(len(x)))+"; SNR= "+str(np.mean(x)/np.std(x)*np.sqrt(len(x)))
#         path = self.pathFig+"/hist_tsz"+str(iRAp)+".pdf"
#         myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'Std dev value [arbitrary]', semilogy=True)


   ##################################################################################
   ##################################################################################
