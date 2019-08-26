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
   
#      self.loadMaps(nProc=self.nProc)
      self.loadAPRadii()
      
      if save:
         self.saveOverlapFlag(nProc=self.nProc)
      self.loadOverlapFlag()
      
      if save:
         self.saveFiltering(nProc=self.nProc)
      self.loadFiltering()
      
      self.measureVarFromHitCount(plot=save)

      if save:
         self.saveTszKsz(nSamples=self.nSamples, nProc=self.nProc)
      self.loadTszKsz()



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
      The value is interpolated between the nearest pixels.
      Will return 0 if the coordinates requested are outside the map
      '''
      # interpolate the map to the given sky coordinates
      sourcecoord = np.array([dec, ra]) * utils.degree
      # use nearest neighbor interpolation
      return map.at(sourcecoord, prefilter=False, mask_nan=False, order=0)
   
   
   ##################################################################################
   
   def saveOverlapFlag(self, nProc=1):
      print "Create overlap flag"
      # find if a given object overlaps with the CMB hit map
      def foverlap(iObj, thresh=1.e-5):
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
      stampMask = enmap.zeros(shape, wcs)
      stampHit = enmap.zeros(shape, wcs)

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
      # Here, I use nearest neighbor interpolation (order=0)
      # these are now numpy arrays: the wcs info is gone
      stampMap[:,:] = self.cmbMap.at(ipos, prefilter=False, mask_nan=False, order=0)
      stampMask[:,:] = self.cmbMask.at(ipos, prefilter=False, mask_nan=False, order=0)
      stampHit[:,:] = self.cmbHit.at(ipos, prefilter=False, mask_nan=False, order=0)
      
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
      The filter is dimensionless:
      Filter = 1 in the disk, - (disk area)/(ring area) in the ring, 0 outside.
      Hence:
      int d^2theta * Filter = 0.
      r0 and r1 are the radius of the disk and ring in radians.
      Output:
      filtMap: [map unit * sr]
      filtMask: [mask unit * sr]
      filtNoiseStdDev: [1/sqrt(hit unit) * sr], ie [std dev * sr] if [hit map] = inverse var
      diskArea: [sr]
      """
      # coordinates of the square map (between -1 and 1 deg, or whatever the size is)
      # local coordinates in rad.
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
      filtNoiseStdDev = np.sqrt(np.sum((pixArea * filter)**2 / stampHit)) # to get the std devs [sr / sqrt(hit unit)]


      if test:
         print "AP filter with disk radius =", r0 * (180.*60./np.pi), "arcmin"

         # count nb of pixels where filter is strictly positive
         nbPix = len(np.where(filter>0.)[0])
         print "- nb of pixels where filter>0: "+str(nbPix)
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
         
         # extract postage stamp around it
         opos, stampMap, stampMask, stampHit = self.extractStamp(ra, dec, dxDeg=0.3, dyDeg=0.3, resArcmin=0.25, proj='cea', test=test)
         
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

   def stack(self, quantity, mask, weights, norm=True):
      '''Stacks the quantity for the objects selected by mask,
      applying weights from weight.
      The quantity can be a single value per object, or an array for each object,
      like e.g. self.filtMap.
      Norm=True if you want to normalize by the sum of the weights.
      '''
      result = np.sum(quantity[mask,:] * weights[mask,np.newaxis], axis=0)
      if norm:
         result /= np.sum(weights[mask,np.newaxis], axis=0)
      return result


   def plotStack(self, stack, stackError=None, ylim=None, name='dT', nameLatex=r'$\delta T$ [$\mu$K]' ):
      """Generic stack plotter.
      """
      fig=plt.figure(0)
      ax=fig.subplot(111)
      #
      if stackError is None:
         ax.plot(self.RApMpch, stack, 'o-')
      else:
         ax.errorbar(self.RApMpch, stack, yerr=stackError, fmt='o-')
      #
      ax.set_xlabel(r'$R$ [cMpc/h]')
      ax.set_ylabel(nameLatex)
#      fig.savefig(self.pathFig+"/stack_"+name+".pdf")
#      fig.clf()
      plt.show()


   ##################################################################################


   def measureVarFromHitCount(self, plot=False):
      """Returns a list of functions, one for each AP filter radius,
      where the function takes filtNoiseStdDev**2 \propto [(map var) * sr^2] as input and returns the
      actual measured filter variance [(map unit)^2 * sr^2].
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

         # define bins
         nBins = 21
#         print "Binning for measuring var from hit count:", x
#         print "Binning for measuring var from hit count:", np.min(x), np.max(x), np.min(y), np.max(y)
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


   def examineHistograms(self, fsAp=[]):
      """fsAp is an optional list of functions of rAp in radians, which return the expected std def of the AP filter
      """
      
#      self.catalogMask(overlap=True, psMask=True, mVir=[1.e6, 1.e17], extraSelection=1.)
#      path = self.pathFig+"/hist_x.pdf"
#      myHistogram(DEC[mask], nBins=71, lim=(-90., 90.), path=path, nameLatex=r'$x$ [km/s]', semilogx=False, doGauss=False)

      # first, keep all objects that overlap, even the masked ones
      mask = self.catalogMask(overlap=True, psMask=False)

      # check that the non-overlapping objects are the ones with the correct DEC
      path = self.pathFig+"/hist_dec_overlap.pdf"
      myHistogram(self.Catalog.DEC[mask], nBins=71, lim=(-30., 90.), path=path, nameLatex=r'Dec [deg]')
      
      # check the values of the filters on the point source mask, to find a relevant cut
      # look at the largest aperture
      x = self.filtMask[mask,-1]
      path = self.pathFig+"/hist_psmaskvalue_before.pdf"
      myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'PS mask value', semilogy=True)
      
      # then  remove the objects that overlap with point sources
      mask = self.catalogMask(overlap=True, psMask=True)

      # redo the mask histogram, to check
      x = self.filtMask[mask,-1]
      path = self.pathFig+"/hist_psmaskvalue_after.pdf"
      myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'PS mask value', semilogy=True)
      
      # Histograms of filter outputs
      for iRAp in range(self.nRAp):
         x = self.filtMap[mask, iRAp]
         path = self.pathFig+"/hist_filtvalue"+str(iRAp)+".pdf"
         S2Theory = []
         for f in fsAp:
            # theory assumes int d^2theta W = 1
            s2Theory = f(self.RApArcmin[iRAp] * np.pi/180./60.)**2
            # so multiply by disk area, but it can vary from object to object
            # we neglect this source of scatter and just keep the mean
            s2Theory *= np.mean(self.diskArea[mask, iRAp])**2
            S2Theory.append(s2Theory)
         myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'AP filter value', S2Theory=S2Theory, doGauss=True, semilogy=True)

   
      # Histograms of noise std dev, from hit counts
      for iRAp in range(self.nRAp):
         x = self.filtNoiseStdDev[mask, iRAp]
         path = self.pathFig+"/hist_noisestddevhit"+str(iRAp)+".pdf"
         myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'Std dev value [arbitrary]', semilogy=True)


      # tSZ / dust
      for iRAp in range(self.nRAp):
         weights = 1. / self.filtNoiseStdDev[mask, iRAp]**2   # need inverse variance, not std dev
         weights /= np.mean(weights)   # to keep the size and units of the weighted AP outputs
         x = self.filtMap[mask, iRAp] * weights
         print "- mean tSZ= "+str(np.mean(x))+"; std on mean= "+str(np.std(x)/np.sqrt(len(x)))+"; SNR= "+str(np.mean(x)/np.std(x)*np.sqrt(len(x)))
         path = self.pathFig+"/hist_tsz"+str(iRAp)+".pdf"
         myHistogram(x, nBins=71, lim=(np.min(x), np.max(x)), path=path, nameLatex=r'Std dev value [arbitrary]', semilogy=True)


   ##################################################################################
   ##################################################################################

   def compareKszEstimators(self):

      # remove the objects that overlap with point sources
      mask = self.catalogMask(overlap=True, psMask=True)

      kSZ1 = np.zeros(self.nRAp)
      skSZ1 = np.zeros(self.nRAp)
      kSZ2 = np.zeros(self.nRAp)
      skSZ2 = np.zeros(self.nRAp)
      kSZ3 = np.zeros(self.nRAp)
      skSZ3 = np.zeros(self.nRAp)
      kSZ4 = np.zeros(self.nRAp)
      skSZ4 = np.zeros(self.nRAp)
      kSZ5 = np.zeros(self.nRAp)
      skSZ5 = np.zeros(self.nRAp)
      kSZ6 = np.zeros(self.nRAp)
      skSZ6 = np.zeros(self.nRAp)
      kSZ7 = np.zeros(self.nRAp)
      skSZ7 = np.zeros(self.nRAp)


      for iRAp in range(self.nRAp):
         t = self.filtMap[mask, iRAp]# - np.mean(self.filtMap[mask, iRAp])
         v = - self.Catalog.vR[mask]# + np.mean(self.Catalog.vR[mask])
         k = self.Catalog.integratedKSZ[mask]# - np.mean(self.Catalog.integratedKSZ[mask])
         s2True = self.fVarFromHitCount[iRAp](self.filtNoiseStdDev[mask, iRAp]**2)
         s2Hit = self.filtNoiseStdDev[mask, iRAp]**2


         # kSZ1: T * v / s2Hit
         num = np.sum(t * v / s2Hit)
         denom = np.sum(k * v / s2Hit)
         s2num = np.sum(s2True * (v / s2Hit)**2)
         #
         kSZ1[iRAp] = num / denom
         skSZ1[iRAp] = np.sqrt(s2num / denom**2)

         # kSZ6: T * v / s2True
         num = np.sum(t * v / s2True)
         denom = np.sum(k * v / s2True)
         s2num = np.sum(s2True * (v / s2True)**2)
         #
         kSZ6[iRAp] = num / denom
         skSZ6[iRAp] = np.sqrt(s2num / denom**2)


         # kSZ2: T * Mh v / s2Hit
         num = np.sum(t * k / s2Hit)
         denom = np.sum(k**2 / s2Hit)
         s2num = np.sum(s2True * (k / s2Hit)**2)
         #
         kSZ2[iRAp] = num / denom
         skSZ2[iRAp] = np.sqrt(s2num / denom**2)


         # subtracting mean
         tNoMean = self.filtMap[mask, iRAp] - np.mean(self.filtMap[mask, iRAp])
         vNoMean = - self.Catalog.vR[mask] + np.mean(self.Catalog.vR[mask])
         kNoMean = self.Catalog.integratedKSZ[mask] - np.mean(self.Catalog.integratedKSZ[mask])
         s2True = self.fVarFromHitCount[iRAp](self.filtNoiseStdDev[mask, iRAp]**2)
         s2Hit = self.filtNoiseStdDev[mask, iRAp]**2

         # kSZ3: T * Mh v / s2Hit, subtracting mean
         num = np.sum(tNoMean * kNoMean / s2Hit)
         denom = np.sum(kNoMean**2 / s2Hit)
         s2num = np.sum(s2True * (kNoMean / s2Hit)**2)
         #
         kSZ3[iRAp] = num / denom
         skSZ3[iRAp] = np.sqrt(s2num / denom**2)

         # kSZ4: T * Mh v / s2True, subtracting mean,
         # and using the measured noise weights
         num = np.sum(tNoMean * kNoMean / s2True)
         denom = np.sum(kNoMean**2 / s2True)
         #
         kSZ4[iRAp] = num / denom
         skSZ4[iRAp] = np.sqrt(1. / denom)

# Preferred estimator:
         # kSZ5: T * v / s2True, subtracting mean,
         # and using the measured noise weights
         num = np.sum(tNoMean * vNoMean / s2True)
         denom = np.sum(kNoMean * vNoMean / s2True)
         s2num = np.sum(s2True * (vNoMean / s2True)**2)
         #
         kSZ5[iRAp] = num / denom
         skSZ5[iRAp] = np.sqrt(s2num / denom**2)

         # kSZ7: T * v / s2Hit, subtracting mean,
         num = np.sum(tNoMean * vNoMean / s2Hit)
         denom = np.sum(kNoMean * vNoMean / s2Hit)
         s2num = np.sum(s2True * (vNoMean / s2Hit)**2)
         #
         kSZ7[iRAp] = num / denom
         skSZ7[iRAp] = np.sqrt(s2num / denom**2)


      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      #self.RApMpch
      ax.errorbar(self.RApArcmin, kSZ1, skSZ1, label=r'$Tv/\sigma^2$')
      ax.errorbar(self.RApArcmin+0.01, kSZ6, skSZ6, label=r'$Tv/\sigma^2$, better noise')
      ax.errorbar(self.RApArcmin+0.02, kSZ7, skSZ7, label=r'$Tv/\sigma^2$, sub. T,v')
      ax.errorbar(self.RApArcmin+0.03, kSZ5, skSZ5, label=r'$Tv/\sigma^2$, sub. T,v, better noise')
      ax.errorbar(self.RApArcmin+0.04, kSZ2, skSZ2, label=r'$TMv/\sigma^2$')
      ax.errorbar(self.RApArcmin+0.05, kSZ3, skSZ3, label=r'$TMv/\sigma^2$, sub. T,Mv')
      ax.errorbar(self.RApArcmin+0.06, kSZ4, skSZ4, label=r'$TMv/\sigma^2$, sub. T,Mv, better noise')
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'$\alpha_\text{kSZ}$')
      ax.set_ylim((0., 2.5))
      #
      # extra abscissa: disk radius in comoving Mpc/h
      ax2 = ax.twiny()
      ticks = ax.get_xticks()
      ax2.set_xticks(ticks)
      newticks = np.array(ticks) * np.pi/(180.*60.) * self.U.bg.comoving_distance(np.mean(self.Catalog.Z[mask]))  # disk radius in Mpc/h
      newticks = np.round(newticks, 2)
      ax2.set_xticklabels(newticks)
      ax2.set_xlim(ax.get_xlim())
      ax2.set_xlabel(r'$R$ [cMpc/h]', fontsize=20)
      ax2.xaxis.set_label_coords(0.5, 1.15)
      #
      # extra ordinate: signal in muK * arcmin^2
      ax3 = ax.twinx()
      ticks = ax.get_yticks()
      ax3.set_yticks(ticks)
      newticks = np.array(ticks) * np.mean(self.Catalog.integratedKSZ[mask]) * (180.*60./np.pi)**2 # [muK*arcmin^2]
      newticks = np.round(newticks, 3)
      ax3.set_yticklabels(newticks)
      ax3.set_ylim(ax.get_ylim())
      ax3.set_ylabel(r'$\langle T_\text{kSZ} \rangle$ [$\mu$K.arcmin$^2$]', fontsize=20)
      ax3.yaxis.set_label_coords(1.2, 0.5)
      #
      path = self.pathFig+"/ksz_comparison.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      #self.RApMpch
      ax.plot(self.RApArcmin, skSZ1, label=r'$Tv/\sigma^2$')
      ax.plot(self.RApArcmin+0.01, skSZ6, label=r'$Tv/\sigma^2$, better noise')
      ax.plot(self.RApArcmin+0.02, skSZ7, label=r'$Tv/\sigma^2$, sub. T,v')
      ax.plot(self.RApArcmin+0.03, skSZ5, label=r'$Tv/\sigma^2$, sub. T,v, better noise')
      ax.plot(self.RApArcmin+0.04, skSZ2, label=r'$TMv/\sigma^2$')
      ax.plot(self.RApArcmin+0.05, skSZ3, label=r'$TMv/\sigma^2$, sub. T,Mv')
      ax.plot(self.RApArcmin+0.06, skSZ4, label=r'$TMv/\sigma^2$, sub. T,Mv, better noise')
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'$\sigma \left(\alpha_\text{kSZ} \right)$')
      ax.set_ylim((0., 1.1))
      #
      # extra abscissa: disk radius in comoving Mpc/h
      ax2 = ax.twiny()
      ticks = ax.get_xticks()
      ax2.set_xticks(ticks)
      newticks = np.array(ticks) * np.pi/(180.*60.) * self.U.bg.comoving_distance(np.mean(self.Catalog.Z[mask]))  # disk radius in Mpc/h
      newticks = np.round(newticks, 2)
      ax2.set_xticklabels(newticks)
      ax2.set_xlim(ax.get_xlim())
      ax2.set_xlabel(r'$R$ [cMpc/h]', fontsize=20)
      ax2.xaxis.set_label_coords(0.5, 1.15)
      #
      # extra ordinate: signal in muK * arcmin^2
      ax3 = ax.twinx()
      ticks = ax.get_yticks()
      ax3.set_yticks(ticks)
      newticks = np.array(ticks) * np.mean(self.Catalog.integratedKSZ[mask]) * (180.*60./np.pi)**2 # [muK*arcmin^2]
      newticks = np.round(newticks, 3)
      ax3.set_yticklabels(newticks)
      ax3.set_ylim(ax.get_ylim())
      ax3.set_ylabel(r'$\langle T_\text{kSZ} \rangle$ [$\mu$K.arcmin$^2$]', fontsize=20)
      ax3.yaxis.set_label_coords(1.2, 0.5)
      #
      path = self.pathFig+"/sksz_comparison.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      #self.RApMpch
      ax.plot(self.RApArcmin, kSZ1/skSZ1, label=r'$Tv/\sigma^2$')
      ax.plot(self.RApArcmin+0.01, kSZ6/skSZ6, label=r'$Tv/\sigma^2$, better noise')
      ax.plot(self.RApArcmin+0.02, kSZ7/skSZ7, label=r'$Tv/\sigma^2$, sub. T,v')
      ax.plot(self.RApArcmin+0.03, kSZ5/skSZ5, label=r'$Tv/\sigma^2$, sub. T,v, better noise')
      ax.plot(self.RApArcmin+0.04, kSZ2/skSZ2, label=r'$TMv/\sigma^2$')
      ax.plot(self.RApArcmin+0.05, kSZ3/skSZ3, label=r'$TMv/\sigma^2$, sub. T,Mv')
      ax.plot(self.RApArcmin+0.06, kSZ4/skSZ4, label=r'$TMv/\sigma^2$, sub. T,Mv, better noise')
      #
      ax.legend(loc=4, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'kSZ SNR')
      ax.set_ylim((0., 6.))
      #
      # extra abscissa: disk radius in comoving Mpc/h
      ax2 = ax.twiny()
      ticks = ax.get_xticks()
      ax2.set_xticks(ticks)
      newticks = np.array(ticks) * np.pi/(180.*60.) * self.U.bg.comoving_distance(np.mean(self.Catalog.Z[mask]))  # disk radius in Mpc/h
      newticks = np.round(newticks, 2)
      ax2.set_xticklabels(newticks)
      ax2.set_xlim(ax.get_xlim())
      ax2.set_xlabel(r'$R$ [cMpc/h]', fontsize=20)
      ax2.xaxis.set_label_coords(0.5, 1.15)
      #
      path = self.pathFig+"/snr_ksz_comparison.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()


   ##################################################################################


   def tszKszEstimator(self, filtMap=None, v=None, k=None, filtNoiseStdDev=None, mask=None):
      """ Measures stacked tSZ and kSZ, as a function of aperture.
      Returns tSZ, stSZ, kSZ, skSZ in [map unit * sr]
      Input:
      filtMap: default is self.filtMap
      vR: default is - self.Catalog.vR / c
      kSZ: default is self.Catalog.integratedKSZ. Not used currently
      mask: default is self.catalogMask(overlap=True, psMask=True)
      """
      if mask is None:
         # remove the objects that overlap with point sources
         mask = self.catalogMask(overlap=True, psMask=True)
      if filtMap is None:
         filtMap = self.filtMap.copy()
      if v is None:
         v = -self.Catalog.vR.copy() / 3.e5  # v/c [dimless]
      if k is None:
         k = self.Catalog.integratedKSZ.copy()
      if filtNoiseStdDev is None:
         filtNoiseStdDev = self.filtNoiseStdDev.copy()
         
#      if np.any(~np.isfinite(v)) or np.any(~np.isfinite(k)):
#         print "problem, input infinite"

      tSZ = np.zeros(self.nRAp)
      stSZ = np.zeros(self.nRAp)
      kSZ = np.zeros(self.nRAp)
      skSZ = np.zeros(self.nRAp)
      for iRAp in range(self.nRAp):
         # AP filter values at the current radius
         t = filtMap[mask, iRAp]
         # subtracting mean
         tNoMean = t - np.mean(t)
         vNoMean = v[mask] - np.mean(v[mask])
         kNoMean = k[mask] - np.mean(k[mask])
         # hit count and measured total noise (CMB + detector)
         s2Hit = filtNoiseStdDev[mask, iRAp]**2
         s2True = self.fVarFromHitCount[iRAp](s2Hit)
         
#         if np.any(~np.isfinite(s2Hit)):
#            print "problem s2Hit", s2Hit[np.where(~np.isfinite(s2Hit))]
#         if np.any(~np.isfinite(s2True)):
#            print "problem s2True", s2True
#         if np.any(~np.isfinite(t)):
#            print "problem t", t
#         if np.any(~np.isfinite(tNoMean)):
#            print "problem tNoMean", tNoMean
#         if np.any(~np.isfinite(vNoMean)):
#            print "problem vNoMean", vNoMean
#         if np.any(~np.isfinite(kNoMean)):
#            print "problem kNoMean", kNoMean

         # tSZ
         # Subtracting mean,
         # and using the measured noise weights
         num = np.sum(t / s2True)
         denom = np.sum(1. / s2True)
         #
         tSZ[iRAp] = num / denom
         stSZ[iRAp] = np.sqrt(1. / denom)
         
         # kSZ
         # T * v / s2True, subtracting mean,
         # and using the measured noise weights
         num = np.sum(tNoMean * vNoMean / s2True)
         # show the result in terms of alpha[dimless]
         #denom = np.sum(kNoMean * vNoMean / s2True)
         # show the result in terms of AP filtered kSZ [muK * sr]
         denom = np.sum(vNoMean**2 / s2True)
         #
         kSZ[iRAp] = num / denom
         skSZ[iRAp] = np.sqrt(1. / denom)
         #
         # rescale by the RMS v/c, to get a fair quantity to compare to tSZ
         kSZ[iRAp] *= np.std(vNoMean)
         skSZ[iRAp] *= np.std(vNoMean)

      return tSZ, stSZ, kSZ, skSZ


   ##################################################################################
   

   def tszKszCovBootstrap(self, nSamples=1000, nProc=1):
      """Estimate tSZ/kSZ covariance matrix from bootstrap resampling.
      """
      # list of all objects, overlapping or not
      I = np.arange(self.Catalog.nObj)
      # remove the objects that overlap with point sources
      mask = self.catalogMask(overlap=True, psMask=True)
      
      def resample(iSample):
         print "starting sample number", iSample
         # make sure each random resample is really independent
         np.random.seed(iSample)
         # resample the overlapping objects from the overlapping objects, with replacement
         # leave the other objects untouched
         J = I.copy()
         J[mask] = np.random.choice(I[mask], size=np.sum(mask), replace=True)
         
#         print "test len(J) =", len(J)
#         print "do masks agree?", np.sum(mask-mask[J])
#         print "test s2Hit before:", np.any(~np.isfinite(self.filtNoiseStdDev[mask,:]))
#         x = self.filtNoiseStdDev[J,:]
#         print "test s2Hit after:", np.any(~np.isfinite(self.filtNoiseStdDev[mask[J],:]))

         # run kSZ estimator on the current resample
         tSZ, stSZ, kSZ, skSZ = self.tszKszEstimator(filtMap=self.filtMap[J,:], v=-self.Catalog.vR[J]/3.e5, k=self.Catalog.integratedKSZ[J], filtNoiseStdDev=self.filtNoiseStdDev[J,:], mask=mask[J])
         
#         if np.any(~np.isfinite(kSZ)):
#            print "problem: kSZ=", kSZ

         return tSZ, stSZ, kSZ, skSZ
      
      tStart = time()
      with sharedmem.MapReduce(np=nProc) as pool:
         result = np.array(pool.map(resample, range(nSamples)))
      tStop = time()
      print "took", (tStop-tStart)/60., "min"
      # unpack results
      tSZSamples = result[:,0,:] # shape (nObj, nRAp)
      stSZSamples = result[:,1,:]
      kSZSamples = result[:,2,:]
      skSZSamples = result[:,3,:]
      # estimate mean and cov
      meanTsz = np.mean(tSZSamples, axis=0)
      covTsz = np.cov(tSZSamples, rowvar=False)
      meanKsz = np.mean(kSZSamples, axis=0)
      covKsz = np.cov(kSZSamples, rowvar=False)
      return meanTsz, covTsz, meanKsz, covKsz


   def tszKszCovShuffleV(self, nSamples=1000, nProc=1):
      """Estimate tSZ/kSZ covariance matrix from shuffling the velocities.
      """
      # list of all objects, overlapping or not
      I = np.arange(self.Catalog.nObj)
      # remove the objects that overlap with point sources
      mask = self.catalogMask(overlap=True, psMask=True)
      
      def resample(iSample):
         print "starting sample number", iSample
         # make sure each random resample is really independent
         np.random.seed(iSample)
         # shuffle only the overlapping objects,
         # leave the other objects untouched
         J = I.copy()
         J[mask] = np.random.permutation(J[mask])
         
#         print "test len(J) =", len(J)
#         print "do masks agree?", np.sum(mask-mask[J])
#         print "test s2Hit before:", np.any(~np.isfinite(self.filtNoiseStdDev[mask,:]))
#         x = self.filtNoiseStdDev[J,:]
#         print "test s2Hit after:", np.any(~np.isfinite(self.filtNoiseStdDev[mask[J],:]))

         # run kSZ estimator on the current resample: shuffle only vr and ksz
         # no need to shuffle the mask, since masked objects are still masked
         tSZ, stSZ, kSZ, skSZ = self.tszKszEstimator(filtMap=self.filtMap, v=-self.Catalog.vR[J], k=self.Catalog.integratedKSZ[J], filtNoiseStdDev=self.filtNoiseStdDev, mask=mask)
         
#         if np.any(~np.isfinite(kSZ)):
#            print "problem: kSZ=", kSZ

         return tSZ, stSZ, kSZ, skSZ
      
      tStart = time()
      with sharedmem.MapReduce(np=nProc) as pool:
         result = np.array(pool.map(resample, range(nSamples)))
      tStop = time()
      print "took", (tStop-tStart)/60., "min"
      # unpack results
      tSZSamples = result[:,0,:] # shape (nObj, nRAp)
      stSZSamples = result[:,1,:]
      kSZSamples = result[:,2,:]
      skSZSamples = result[:,3,:]
      # estimate mean and cov
      meanTsz = np.mean(tSZSamples, axis=0)
      covTsz = np.cov(tSZSamples, rowvar=False)
      meanKsz = np.mean(kSZSamples, axis=0)
      covKsz = np.cov(kSZSamples, rowvar=False)
      return meanTsz, covTsz, meanKsz, covKsz


   ##################################################################################


   def saveTszKsz(self, nSamples=1000, nProc=1):
      # tSZ/kSZ signal and estimated variance
      tSZ, stSZ, kSZ, skSZ = self.tszKszEstimator()
      data = np.zeros((self.nRAp,5))
      data[:,0] = self.RApArcmin # [arcmin]
      data[:,1] = tSZ # [map unit * sr]
      data[:,2] = stSZ # [map unit * sr]
      data[:,3] = kSZ # [map unit * sr]
      data[:,4] = skSZ # [map unit * sr]
      np.savetxt(self.pathOut+"/tsz_ksz.txt", data)
      
      # null test and cov mat from bootstrap
      meanTsz, covTsz, meanKsz, covKsz = self.tszKszCovBootstrap(nSamples=nSamples, nProc=nProc)
      np.savetxt(self.pathOut+"/null_tsz_bootstrap.txt", meanTsz)
      np.savetxt(self.pathOut+"/cov_tsz_bootstrap.txt", covTsz)
      np.savetxt(self.pathOut+"/null_ksz_bootstrap.txt", meanKsz)
      np.savetxt(self.pathOut+"/cov_ksz_bootstrap.txt", covKsz)
   
      # null test and cov mat from shuffling the velocities
      meanTsz, covTsz, meanKsz, covKsz = self.tszKszCovShuffleV(nSamples=nSamples, nProc=nProc)
      np.savetxt(self.pathOut+"/null_tsz_shufflev.txt", meanTsz)
      np.savetxt(self.pathOut+"/cov_tsz_shufflev.txt", covTsz)
      np.savetxt(self.pathOut+"/null_ksz_shufflev.txt", meanKsz)
      np.savetxt(self.pathOut+"/cov_ksz_shufflev.txt", covKsz)
   
   
   def loadTszKsz(self, plot=False):
      data = np.genfromtxt(self.pathOut+"/tsz_ksz.txt")
      self.tSZ = data[:,1] # [map unit * sr]
      self.stSZ = data[:,2] # [map unit * sr]
      self.kSZ = data[:,3] # [map unit * sr]
      self.skSZ = data[:,4] # [map unit * sr]

      self.tSZNullBootstrap = np.genfromtxt(self.pathOut+"/null_tsz_bootstrap.txt")
      self.covTszBootstrap = np.genfromtxt(self.pathOut+"/cov_tsz_bootstrap.txt")
      #
      self.tSZNullShuffleV = np.genfromtxt(self.pathOut+"/null_tsz_shufflev.txt")
      self.covTszShuffleV = np.genfromtxt(self.pathOut+"/cov_tsz_shufflev.txt")
      #
      self.covTsz = self.covTszBootstrap.copy()

      self.kSZNullBootstrap = np.genfromtxt(self.pathOut+"/null_ksz_bootstrap.txt")
      self.covKszBootstrap = np.genfromtxt(self.pathOut+"/cov_ksz_bootstrap.txt")
      #
      self.kSZNullShuffleV = np.genfromtxt(self.pathOut+"/null_ksz_shufflev.txt")
      self.covKszShuffleV = np.genfromtxt(self.pathOut+"/cov_ksz_shufflev.txt")
      #
      self.covKsz = self.covKszBootstrap.copy()


   ##################################################################################


   def plotTszKsz(self):
      # tSZ
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # convert from sr to arcmin^2
      factor = (180.*60./np.pi)**2
      #
      ax.errorbar(self.RApArcmin, factor * self.tSZ, factor * np.sqrt(np.diag(self.covTsz)), fmt='-', c='r')
      #
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'$T_\text{tSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
      #ax.set_ylim((0., 2.))
      #
      path = self.pathFig+"/tsz.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      
      # kSZ
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # convert from sr to arcmin^2
      factor = (180.*60./np.pi)**2
      #
      ax.errorbar(self.RApArcmin, factor * self.kSZ, factor * np.sqrt(np.diag(self.covKsz)), fmt='-', c='r')
      #
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'$T_\text{kSZ}$ [$\mu K\cdot\text{arcmin}^2$]')
      #ax.set_ylim((0., 2.))
      #
      path = self.pathFig+"/ksz.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()



   def plotCovTszKsz(self):
      
      # Bootstrap tSZ
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # pcolor wants x and y to be edges of cell,
      # ie one more element, and offset by half a cell
      dR = (self.rApMaxArcmin - self.rApMinArcmin) / self.nRAp
      RApEdgesArcmin = np.linspace(self.rApMinArcmin-0.5*dR, self.rApMaxArcmin+0.5*dR, self.nRAp+1)
      #
      X, Y = np.meshgrid(RApEdgesArcmin, RApEdgesArcmin, indexing='ij')
      cp=ax.pcolormesh(X, Y, self.covTszBootstrap, cmap='YlOrRd')
      #
      ax.set_aspect('equal')
      plt.colorbar(cp)
      ax.set_xlim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_ylim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_xlabel(r'R [arcmin]')
      ax.set_ylabel(r'R [arcmin]')
      #
      path = self.pathFig+"/cov_tsz_bootstrap.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
      
      # Bootstrap kSZ
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # pcolor wants x and y to be edges of cell,
      # ie one more element, and offset by half a cell
      dR = (self.rApMaxArcmin - self.rApMinArcmin) / self.nRAp
      RApEdgesArcmin = np.linspace(self.rApMinArcmin-0.5*dR, self.rApMaxArcmin+0.5*dR, self.nRAp+1)
      #
      X, Y = np.meshgrid(RApEdgesArcmin, RApEdgesArcmin, indexing='ij')
      cp=ax.pcolormesh(X, Y, self.covKszBootstrap, cmap='YlOrRd')
      #
      ax.set_aspect('equal')
      plt.colorbar(cp)
      ax.set_xlim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_ylim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_xlabel(r'R [arcmin]')
      ax.set_ylabel(r'R [arcmin]')
      #
      path = self.pathFig+"/cov_ksz_bootstrap.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()

      # Bootstrap tSZ corr coeff
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # pcolor wants x and y to be edges of cell,
      # ie one more element, and offset by half a cell
      dR = (self.rApMaxArcmin - self.rApMinArcmin) / self.nRAp
      RApEdgesArcmin = np.linspace(self.rApMinArcmin-0.5*dR, self.rApMaxArcmin+0.5*dR, self.nRAp+1)
      X, Y = np.meshgrid(RApEdgesArcmin, RApEdgesArcmin, indexing='ij')
      #
      cov = self.covTszBootstrap.copy()
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
      path = self.pathFig+"/cor_tsz_bootstrap.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()

      # Bootstrap kSZ corr coeff
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # pcolor wants x and y to be edges of cell,
      # ie one more element, and offset by half a cell
      dR = (self.rApMaxArcmin - self.rApMinArcmin) / self.nRAp
      RApEdgesArcmin = np.linspace(self.rApMinArcmin-0.5*dR, self.rApMaxArcmin+0.5*dR, self.nRAp+1)
      X, Y = np.meshgrid(RApEdgesArcmin, RApEdgesArcmin, indexing='ij')
      #
      cov = self.covKszBootstrap.copy()
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
      path = self.pathFig+"/cor_ksz_bootstrap.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()





      # Shuffle v covariance
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # pcolor wants x and y to be edges of cell,
      # ie one more element, and offset by half a cell
      dR = (self.rApMaxArcmin - self.rApMinArcmin) / self.nRAp
      RApEdgesArcmin = np.linspace(self.rApMinArcmin-0.5*dR, self.rApMaxArcmin+0.5*dR, self.nRAp+1)
      #
      X, Y = np.meshgrid(RApEdgesArcmin, RApEdgesArcmin, indexing='ij')
      cp=ax.pcolormesh(X, Y, self.covKszShuffleV, cmap='YlOrRd')
      #
      ax.set_aspect('equal')
      plt.colorbar(cp)
      ax.set_xlim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_ylim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_xlabel(r'R [arcmin]')
      ax.set_ylabel(r'R [arcmin]')
      #
      path = self.pathFig+"/cov_ksz_shufflev.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()


      # Shuffle v corr coeff
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # pcolor wants x and y to be edges of cell,
      # ie one more element, and offset by half a cell
      dR = (self.rApMaxArcmin - self.rApMinArcmin) / self.nRAp
      RApEdgesArcmin = np.linspace(self.rApMinArcmin-0.5*dR, self.rApMaxArcmin+0.5*dR, self.nRAp+1)
      X, Y = np.meshgrid(RApEdgesArcmin, RApEdgesArcmin, indexing='ij')
      #
      cov = self.covKszShuffleV.copy()
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
      path = self.pathFig+"/cor_ksz_shufflev.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()



      # Compare bootstrap VS Shuffle v for kSZ
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # pcolor wants x and y to be edges of cell,
      # ie one more element, and offset by half a cell
      dR = (self.rApMaxArcmin - self.rApMinArcmin) / self.nRAp
      RApEdgesArcmin = np.linspace(self.rApMinArcmin-0.5*dR, self.rApMaxArcmin+0.5*dR, self.nRAp+1)
      #
      X, Y = np.meshgrid(RApEdgesArcmin, RApEdgesArcmin, indexing='ij')
      cp=ax.pcolormesh(X, Y, self.covKszBootstrap - self.covKszShuffleV, cmap='YlOrRd')
      #
      ax.set_aspect('equal')
      plt.colorbar(cp)
      ax.set_xlim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_ylim((np.min(RApEdgesArcmin), np.max(RApEdgesArcmin)))
      ax.set_xlabel(r'R [arcmin]')
      ax.set_ylabel(r'R [arcmin]')
      #
      path = self.pathFig+"/cov_ksz_bootstrap-shufflev.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
   
   
      # Diagonal: compare bootstrap, shuffled velocities
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.plot(self.RApArcmin, self.skSZ, label=r'semi-analytic')
      ax.plot(self.RApArcmin, np.sqrt(np.diag(self.covKszBootstrap)), label=r'bootstrap')
      ax.plot(self.RApArcmin, np.sqrt(np.diag(self.covKszShuffleV)), label=r'shuffled v')
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'$\sigma\left(\alpha_\text{kSZ} \right)$')
      #
      path = self.pathFig+"/sksz_bootstrap_shufflev_semianalytic.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()
   
   
   ##################################################################################


   def ftheoryGaussianProfile(self, sigma_cluster):
      """Alpha_ksz signal, between 0 and 1.
      Assumes that the projected cluster profile is a 2d Gaussian,
      with sigma_cluster in arcmin
      Assumes equal area disk-ring filter
      """
      return (1. - np.exp(-0.5*self.RApArcmin**2/sigma_cluster**2))**2


   def computeSnrKsz(self):

      path = self.pathFig+"/snr_ksz.txt"
      with open(path, 'w') as f:
         f.write("*** kSZ SNR ***\n")
      
         # Compute chi^2_null
         chi2Null = self.kSZ.dot( np.linalg.inv(self.covKsz).dot(self.kSZ) )
         # goodness of fit for null hypothesis
         f.write("number of dof:"+str(len(self.kSZ))+"\n")
         f.write("null chi2Null="+str(chi2Null)+"\n")
         pteNull = 1.- stats.chi2.cdf(chi2Null, len(self.kSZ))
         f.write("null pte="+str(pteNull)+"\n")
         # pte as a function of sigma, for a Gaussian random variable
         fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.)) - pteNull
         sigmaNull = optimize.brentq(fsigmaToPTE , 0., 50.)
         f.write("null sigma significance="+str(sigmaNull)+"sigmas\n\n")

         # Gaussian model: find best fit amplitude
         sigma_cluster = 1.5  # arcmin
         theory = self.ftheoryGaussianProfile(sigma_cluster)
         def fchi2(p):
            a = p[0]
            result = (self.kSZ-a*theory).dot( np.linalg.inv(self.covKsz).dot(self.kSZ-a*theory) )
            result -= chi2Null
            return result
         # Minimize the chi squared
         p0 = 1.
         res = optimize.minimize(fchi2, p0)
         #f.write( res
         abest = res.x[0]
         #sbest= res.x[1]
         f.write("best-fit amplitude="+str(abest)+"\n")
         f.write("number of dof:"+str(len(self.kSZ) - 1)+"\n\n")

         # goodness of fit for best fit
         chi2Best = fchi2([abest])+chi2Null
         f.write("best-fit chi2="+str(chi2Best)+"\n")
         pteBest = 1.- stats.chi2.cdf(chi2Best, len(self.kSZ)-1.)
         f.write("best-fit pte="+str(pteBest)+"\n")
         # pte as a function of sigma, for a Gaussian random variable
         fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.)) - pteBest
         sigma = optimize.brentq(fsigmaToPTE , 0., 50.)
         f.write("pte sigma significance="+str(sigma)+"sigmas\n\n")

         # favour of best fit over null
         f.write("fiducial delta chi2="+str(fchi2([1.]))+"\n")
         f.write("best-fit sqrt(delta chi2)="+str(np.sqrt(abs(fchi2([abest]))))+"sigmas\n")
         fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.))
         pte = fsigmaToPTE( np.sqrt(abs(fchi2([abest]))) )
         f.write("pte (if Gaussian)="+str(pte)+"\n")



   def computeSnrTsz(self):

      path = self.pathFig+"/snr_tsz.txt"
      with open(path, 'w') as f:
         f.write("*** tSZ SNR ***\n")
      
         # Compute chi^2_null
         chi2Null = self.tSZ.dot( np.linalg.inv(self.covTsz).dot(self.tSZ) )
         # goodness of fit for null hypothesis
         f.write("number of dof:"+str(len(self.tSZ))+"\n")
         f.write("null chi2Null="+str(chi2Null)+"\n")
         pteNull = 1.- stats.chi2.cdf(chi2Null, len(self.tSZ))
         f.write("null pte="+str(pteNull)+"\n")
         # pte as a function of sigma, for a Gaussian random variable
         fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.)) - pteNull
         sigmaNull = optimize.brentq(fsigmaToPTE , 0., 50.)
         f.write("null sigma significance="+str(sigmaNull)+"sigmas\n\n")

         # Gaussian model: find best fit amplitude
         sigma_cluster = 1.5  # arcmin
         theory = self.ftheoryGaussianProfile(sigma_cluster)
         def fchi2(p):
            a = p[0]
            result = (self.tSZ-a*theory).dot( np.linalg.inv(self.covTsz).dot(self.tSZ-a*theory) )
            result -= chi2Null
            return result
         # Minimize the chi squared
         p0 = 1.
         res = optimize.minimize(fchi2, p0)
         #f.write( res
         abest = res.x[0]
         #sbest= res.x[1]
         f.write("best-fit amplitude="+str(abest)+"\n")
         f.write("number of dof:"+str(len(self.tSZ) - 1)+"\n\n")

         # goodness of fit for best fit
         chi2Best = fchi2([abest])+chi2Null
         f.write("best-fit chi2="+str(chi2Best)+"\n")
         pteBest = 1.- stats.chi2.cdf(chi2Best, len(self.tSZ)-1.)
         f.write("best-fit pte="+str(pteBest)+"\n")
         # pte as a function of sigma, for a Gaussian random variable
         fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.)) - pteBest
         sigma = optimize.brentq(fsigmaToPTE , 0., 50.)
         f.write("pte sigma significance="+str(sigma)+"sigmas\n\n")

         # favour of best fit over null
         f.write("fiducial delta chi2="+str(fchi2([1.]))+"\n")
         f.write("best-fit sqrt(delta chi2)="+str(np.sqrt(abs(fchi2([abest]))))+"sigmas\n")
         fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.))
         pte = fsigmaToPTE( np.sqrt(abs(fchi2([abest]))) )
         f.write("pte (if Gaussian)="+str(pte)+"\n")



   ##################################################################################


   def kszNullTests(self):
      """Compute and plot null tests.
      """
      
      # remove the objects that overlap with point sources
      mask = self.catalogMask(overlap=True, psMask=True)
      
      # tSZ bias: tSZ contamination to the kSZ signal (ignore profile shape)
      # expected y signal for each object, based on mass
      t = np.column_stack((self.Catalog.integratedY for i in range(self.nRAp)))   # copy for each AP radius
      # convert from y to dT
      Tcmb = 2.726
      h = 6.63e-34
      kB = 1.38e-23
      def f(nu):
         """frequency dependence for tSZ temperature
         """
         x = h*nu/(kB*Tcmb)
         return x*(np.exp(x)+1.)/(np.exp(x)-1.) -4.
      t *= f(150.e9) * Tcmb * 1.e6  # muK
      print np.mean(t), np.std(t)
      # run into the estimator
      tSZ, stSZ, tSZBias, stSZBias = self.tszKszEstimator(filtMap=t)

      # basic estimator check: should give 1 for all apertures (ignore profile shape)
      t = np.column_stack((self.Catalog.integratedKSZ for i in range(self.nRAp)))   # copy for each AP radius
      tSZ, stSZ, kSZBias, skSZBias = self.tszKszEstimator(filtMap=t)

#      # Sign Null: replace v by |v|, should get zero
#      # This does null the kSZ, but no longer nulls tSZ and dust
#      # --> not a very useful null test
#      kSZNullAbsV, skSZNullAbsV = self.tszKszEstimator(v=np.abs(self.Catalog.vR), k=np.abs(self.Catalog.integratedKSZ))



      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # 1-sigma error band
      ax.fill_between(self.RApArcmin, -np.sqrt(np.diag(self.covKszBootstrap)), np.sqrt(np.diag(self.covKszBootstrap)), edgecolor='', facecolor='gray', alpha=0.2)
      #
      # Bootstrap convergence test?
      # Bootstrap mean should match measured signal
      ax.errorbar(self.RApArcmin+0.03, self.kSZNullBootstrap-self.kSZ, np.sqrt(np.diag(self.covKszBootstrap))/np.sqrt(self.nSamples), fmt='--', label=r'Bootstrap convergence')
      #
      # Estimator bias: should return 1 when run on the expected signal
      ax.errorbar(self.RApArcmin+0.01, kSZBias-1., fmt='--', label=r'Estimator bias')
      #
      # tSZ bias
      ax.errorbar(self.RApArcmin, tSZBias, fmt='--', label=r'tSZ bias')
      print tSZBias
      #
#      ax.errorbar(self.RApArcmin+0.02, kSZNullAbsV, label=r'$v \leftarrow |v|$')
      #
      # Shuffle v test: should trivially give zero
      ax.errorbar(self.RApArcmin+0.04, self.kSZNullShuffleV, np.sqrt(np.diag(self.covKszShuffleV))/np.sqrt(self.nSamples), label=str(self.nSamples)+r' v shuffles')
      #
      ax.legend(loc=2, fontsize='x-small', labelspacing=0.1)
      ax.set_xlabel(r'$R$ [arcmin]')
      ax.set_ylabel(r'$\alpha_\text{kSZ}$')
      ax.set_ylim((-0.45, 0.45))
      #
      # extra abscissa: disk radius in comoving Mpc/h
      ax2 = ax.twiny()
      ticks = ax.get_xticks()
      ax2.set_xticks(ticks)
      newticks = np.array(ticks) * np.pi/(180.*60.) * self.U.bg.comoving_distance(np.mean(self.Catalog.Z[mask]))  # disk radius in Mpc/h
      newticks = np.round(newticks, 2)
      ax2.set_xticklabels(newticks)
      ax2.set_xlim(ax.get_xlim())
      ax2.set_xlabel(r'$R$ [cMpc/h]', fontsize=20)
      ax2.xaxis.set_label_coords(0.5, 1.15)
      #
      # extra ordinate: signal in muK * arcmin^2
      ax3 = ax.twinx()
      ticks = ax.get_yticks()
      ax3.set_yticks(ticks)
      newticks = np.array(ticks) * np.mean(self.Catalog.integratedKSZ[mask]) * (180.*60./np.pi)**2 # [muK*arcmin^2]
      newticks = np.round(newticks, 3)
      ax3.set_yticklabels(newticks)
      ax3.set_ylim(ax.get_ylim())
      ax3.set_ylabel(r'$\langle T_\text{kSZ} \rangle$ [$\mu$K.arcmin$^2$]', fontsize=20)
      ax3.yaxis.set_label_coords(1.2, 0.5)
      #
      path = self.pathFig+"/null_ksz.pdf"
      fig.savefig(path, bbox_inches='tight')
      fig.clf()



   ##################################################################################
   ##################################################################################

