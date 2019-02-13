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

   
#      self.loadMaps(nProc=self.nProc)
      self.loadAPRadii()
      
      if save:
         self.saveOverlapFlag(nProc=self.nProc)
      self.loadOverlapFlag()
      
      if save:
         self.saveFiltering(nProc=self.nProc)
      self.loadFiltering()





   ##################################################################################
   ##################################################################################
   
   def loadAPRadii(self):
   
      # radii to use for AP filter: comoving Mpc/h
      self.nRAp = 4  #15
      
      # Aperture radii in Mpc/h
      self.rApMinMpch = 1.
      self.rApMaxMpch = 5
      self.RApMpch = np.linspace(self.rApMinMpch, self.rApMaxMpch, self.nRAp)
      
      # Aperture radii in arcmin
      self.rAnMinArcmin = 1.
      self.rAnMaxArcmin = 4.
      self.RApArcmin = np.linspace(self.rAnMinArcmin, self.rAnMaxArcmin, self.nRAp)


   ##################################################################################
   
   # check the hit map as a function of dec, at fixed ra = 0
   # this uses interpolation with the nearest pixels
   def sky2map(self, ra, dec, map):
      '''Gives the map value at coordinates (ra, dec).
      The value is interpolated between the nearest pixels.
      Will return 0 if the coordinates requested are outside the map
      '''
      # interpolate the map to the given sky coordinates
      sourcecoord = np.array([dec, ra]) * utils.degree
      return map.at(sourcecoord, prefilter=False, mask_nan=False)
   
   
   
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
   

   def extractStamp(self, ra, dec, dxDeg=0.25, dyDeg=0.25, resArcmin=0.25, proj='cea', test=False):
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
      # these are now numpy arrays: the wcs info is gone
      stampMap[:,:] = self.cmbMap.at(ipos, prefilter=False, mask_nan=False)
      stampMask[:,:] = self.cmbMask.at(ipos, prefilter=False, mask_nan=False)
      stampHit[:,:] = self.cmbHit.at(ipos, prefilter=False, mask_nan=False)

      if test:
         print "- plot the map"
         plots=enplot.get_plots(stampMap,grid=True)
         enplot.write(self.pathTestFig+"/stampmap_ra"+np.str(np.round(ra, 2))+"_dec"+np.str(np.round(dec, 2)), plots)
         print "- plot the mask"
         plots=enplot.get_plots(stampMask,grid=True)
         enplot.write(self.pathTestFig+"/stampmask_ra"+np.str(np.round(ra, 2))+"_dec"+np.str(np.round(dec, 2)), plots)
         print "- plot the hit"
         plots=enplot.get_plots(stampHit,grid=True)
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

      # apply the filter
      filtMap = np.sum(pixArea * filter * stampMap)   # [map unit * sr]
      
      # detect point sources within the filter:
      # gives 0 in the absence of point sources/edges; gives >=1 in the presence of point sources/edges
      filtMask = np.sum((radius<=r1) * (1-stampMask))   # [dimensionless]
      
      # quantify noise std dev in the filter
      filtNoiseStdDev = np.std(np.sum((pixArea * filter)**2 / stampHit)) # to get the variance [sr / hit unit]


      if test:
         '''
         print "- plot the map"
         plots=enplot.plot(stampMap,grid=True)
         plots.write(self.pathFig+"/stampmap_obj"+str(iObj)+".pdf")
         print "- plot the mask"
         plots=enplot.plot(stampMask,grid=True)
         plots.write(self.pathFig+"/stampmask_obj"+str(iObj)+".pdf")
         print "- plot the hit"
         plots = enplot.plot(stampHit,grid=True)
         plots.write(self.pathFig+"/stamphit_obj"+str(iObj)+".pdf")
         '''
         print "- plot the filter"
         filterMap = stampMap.copy()
         filterMap[:,:] = filter.copy()
         plots=enplot.get_plots(filterMap,grid=True)
         enplot.write(self.pathTestFig+"/stampfilter_r0"+floatExpForm(r0)+"_r1"+floatExpForm(r1), plots)
         

         # count nb of pixels where filter is strictly positive
         nbPix = len(np.where(filter>0.)[0])
         print "- nb of pixels where filter>0: "+str(nbPix)
         print "  ie area of "+str(diskArea)+" sr"
         
         print "- disk-ring filter sums over pixels to "+str(np.sum(filter))
         print "  (should be 0; to be compared with "+str(len(filter.flatten()))+")"

         print "- filter on map:"+str(filtMap)
         print "- filter on mask:"+str(filtMask)
         print "- filter on inverse hit:"+str(filtNoiseStdDev)

      return filtMap, filtMask, filtNoiseStdDev, diskArea


   ##################################################################################


   # analysis to be done for each object
   def analyzeObject(self, iObj, test=False):
      
      if iObj%1000==0:
         print "-", iObj
      
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
         opos, stampMap, stampMask, stampHit = self.extractStamp(ra, dec, dxDeg=0.25, dyDeg=0.25, resArcmin=0.25, proj='cea', test=test)
         
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
      
      
#      # analysis to be done for each object
#      def analyzeObject(iObj):
#
#         if iObj%1000==0:
#            print "-", iObj
#
#         # create arrays of filter values for the given object
#         filtMap = np.zeros(self.nRAp)
#         filtMask = np.zeros(self.nRAp)
#         filtNoiseStdDev = np.zeros(self.nRAp)
#         diskArea = np.zeros(self.nRAp)
#
#         # only do the analysis if the object overlaps with the CMB map
#         if self.overlapFlag[iObj]:
#            # Object coordinates
#            ra = self.Catalog.RA[iObj]   # in deg
#            dec = self.Catalog.DEC[iObj] # in deg
#            z = self.Catalog.Z[iObj]
#
#            # extract postage stamp around it
#            opos, stampMap, stampMask, stampHit = self.extractStamp(ra, dec, dxDeg=0.25, dyDeg=0.25, resArcmin=0.25, proj='cea')
#
#            # loop over the radii for the AP filter
#            for iRAp in range(self.nRAp):
#               # disk radius in comoving Mpc/h
#               rApMpch = self.RApMpch[iRAp]
#               # convert to radians at the given redshift
#               r0 = rApMpch / self.U.bg.comoving_transverse_distance(z) # rad
#               # choose an equal area AP filter
#               r1 = r0 * np.sqrt(2.)
#               # perform the filtering
#               filtMap[iRAp], filtMask[iRAp], filtNoiseStdDev[iRAp], diskArea[iRAp] = self.diskRingFilter(opos, stampMap, stampMask, stampHit, r0, r1, test=False)
#
#         return filtMap, filtMask, filtNoiseStdDev, diskArea


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
      '''Returns catalog mask: 1 for objects to discard, 0 for objects to keep.
      Use as:
      maskedQuantity = Quantity[~mask]
      '''
      # Here mask is 1 for objects we want to keep
      mask = np.ones_like(self.Catalog.RA)
      if mVir is not None:
         mask *= (self.Catalog.Mvir>=mVir[0]) * (self.Catalog.Mvir<=mVir[1])
      if overlap:
         mask *= self.overlapFlag.copy()
      # PS mask: look at largest aperture, and remove if any point within the disk or ring is masked
      if psMask:
         mask *= 1-self.filtMap[:, -1]
      mask *= extraSelection

      # Here mask is 1 for objects to discard (convenient for inversion)
      mask = 1 - mask
      mask = mask.astype(bool)
      return mask


   ##################################################################################

   def stack(self, quantity, mask, weights, norm=False):
      '''Stacks the quantity for the objects selected by mask,
      applying weights from weight.
      The quantity can be a single value per object, or an array for each object,
      like e.g. self.filtMap.
      Norm=True if you want to normalize by the sum of the weights.
      '''
      result = np.sum(quantity[~mask,:] * weights[~mask,np.newaxis], axis=0)
      if norm:
         result /= np.sum(weights[~mask,np.newaxis], axis=0)
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


   def histogram(self, X, nBins=71, lim=(-1000., 1000.), sigma2Theory=None, name='x', nameLatex=r'$x$ [km/s]', semilogx=False, semilogy=False, doGauss=False):
      """Generic histogram plotter.
      """
      # Bin edges
      if semilogx:
         Bins = np.logspace(np.log10(lim[0]), np.log10(lim[1]), nBins, 10.)
      else:
         Bins = np.linspace(lim[0], lim[1], nBins)
      binwidth = Bins[1:] - Bins[:-1]

      # Data histogram
      histX = np.histogram(X, Bins)[0]
      histX = histX.astype(float)

      # histogram for a Gaussian with the variance from the data
      if doGauss:
         mean = np.mean(X)
         std = np.std(X)
         sigma2 = std**2
         av = mean
         fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
         g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-3)[0]
         histGaussFit = np.array(map(g, range(nBins-1)))
         histGaussFit *= self.Catalog.nObj

      # Theory histogram
      if sigma2Theory is not None:
         av = 0.
         fPDF = lambda x: (2*np.pi*sigma2Theory)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2Theory))
         g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-3)[0]
         histTheory = np.array(map(g, range(nBins-1)))
         histTheory *= self.Catalog.nObj

      # Plot
      fig = plt.figure(0)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histX, binwidth, color='b', alpha=0.5, label=r'Data')
      if doGauss:
         ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      if sigma2Theory is not None:
         ax.step(Bins[:-1], histTheory, color='r', lw=3, where='post', label=r'Theory')
      #
      ax.legend(loc=1)
      ax.set_xlim((lim[0], lim[1]))
      if semilogx:
         ax.set_xscale('log', nonposx='clip')
      if semilogy:
         ax.set_yscale('log', nonposy='clip')
      ax.set_xlabel(nameLatex)
      ax.set_ylabel(r'number of objects')
      fig.savefig(self.pathFig+"/hist_"+name+".pdf", bbox_inches='tight')
      fig.clf()
#      plt.show()


   ##################################################################################


   def examineHistograms(self):
      
      
#      self.catalogMask(overlap=True, psMask=True, mVir=[1.e6, 1.e17], extraSelection=1.)
#      self.histogram(DEC[~mask], nBins=71, lim=(-90., 90.), sigma2Theory=None, name='x', nameLatex=r'$x$ [km/s]', semilogx=False, doGauss=False)

      # first, keep all objects that overlap, even the masked ones
      mask = self.catalogMask(overlap=True, psMask=False)

      # check that the non-overlapping objects are the ones with the correct DEC
      self.histogram(self.Catalog.DEC[~mask], nBins=71, lim=(-30., 90.), name='dec_overlap', nameLatex=r'Dec [deg]')
      
      # check the values of the filters on the point source mask, to find a relevant cut
      # look at the largest aperture
      x = self.filtMask[~mask,-1]
      self.histogram(x, nBins=71, lim=(np.min(x), np.max(x)), name='psmaskvalue_before', nameLatex=r'PS mask value', semilogy=True)
      
      # then  remove the objects that overlap with point sources
      mask = self.catalogMask(overlap=True, psMask=True)

      # redo the mask histogram, to check
      x = self.filtMask[~mask,-1]
      self.histogram(x, nBins=71, lim=(-1., 1.), name='psmaskvalue_after', nameLatex=r'PS mask value', semilogy=True)
      
      # is the scatter in T reasonable? Are there outliers?
      for iRAp in range(self.nRAp):
         x = self.filtMap[~mask, iRAp]
         self.histogram(x, nBins=71, lim=(np.min(x), np.max(x)), name='filtvalue'+str(iRAp), nameLatex=r'AP filter value', semilogy=True)
      
      
      # is the kSZ signal visible by eye from the histogram? Probably not.
      # what about the tSZ signal?
      # are there outliers?
      for iRAp in range(self.nRAp):
         x = self.filtMap[~mask, iRAp] * self.Catalog.vR[~mask]
         self.histogram(x, nBins=71, lim=(np.min(x), np.max(x)), name='tv'+str(iRAp), nameLatex=r'$T\times v_r$ [$\mu $K $\times$ km/s]', semilogy=True)
      
      pass




