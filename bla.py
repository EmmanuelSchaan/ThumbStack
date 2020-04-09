#from headers import *
#
#class Class(object):
#
#   def __init__(self, name):
#      self.name = name
#      return
#
#   def fprint(self):
#      print self.name
#
#   def func(self, ts=None):
#      if ts is None:
#         ts = self
#      ts.fprint()
#
#
#
#
#ts1 = Class("1")
#ts2 = Class("2")
#
#
#ts1.func()
#ts1.func(ts=ts2)
#










   def computeStackedProfile(self, filterType, est, iBootstrap=None, iVShuffle=None, tTh=None, stackedMap=False, mVir=None, z=[0., 100.]):
      """Returns the estimated profile and its uncertainty for each aperture.
      est: string to select the estimator
      iBootstrap: index for bootstrap resampling
      iVShuffle: index for shuffling velocities
      tTh: to replace measured temperatures by a theory expectation
      """
      if mVir is None:
         mVir = [self.mMin, self.mMax]

      # select objects that overlap, and reject point sources
      mask = self.catalogMask(overlap=True, psMask=True, filterType=filterType, mVir=mVir, z=z)

      # temperatures [muK * sr]
      if tTh is None:
         t = self.filtMap[filterType].copy() # [muK * sr]
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
      # -v/c [dimless]
      v = -self.Catalog.vR[mask] / 3.e5
      v -= np.mean(v)
      #true filter variance for each object and aperture,
      # valid whether or not a hit count map is available
      s2Full = self.filtVarTrue[filterType][mask, :]
      # Variance from hit count (if available)
      s2Hit = self.filtHitNoiseStdDev[filterType][mask, :]**2
      #print "Shape of s2Hit = ", s2Hit.shape
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
         t -= np.mean(t, axis=0)
         weights = v[:,np.newaxis] * np.ones_like(s2Hit)
         norm = np.std(v) / np.sum(v[:,np.newaxis]*weights, axis=0)
      # kSZ: detector-noise weighted (hit count)
      elif est=='ksz_hitweight':
         # remove mean temperature
         t -= np.mean(t, axis=0)
         weights = v[:,np.newaxis] / s2Hit
         norm = np.std(v) / np.sum(v[:,np.newaxis]*weights, axis=0)
      # kSZ: full noise weighted (detector noise + CMB)
      elif est=='ksz_varweight':
         # remove mean temperature
         t -= np.mean(t, axis=0)
         weights = v[:,np.newaxis] / s2Full
         norm = np.std(v) / np.sum(v[:,np.newaxis]*weights, axis=0)
      # kSZ: full noise weighted (detector noise + CMB)
      elif est=='ksz_massvarweight':
         # remove mean temperature
         t -= np.mean(t, axis=0)
         weights = m[:,np.newaxis] * v[:,np.newaxis] / s2Full
         norm = np.mean(m) * np.std(v) / np.sum(m[:,np.newaxis]**2 * v[:,np.newaxis]**2 / s2Full, axis=0)

      # return the stacked profiles
      if not stackedMap:
         stack = norm * np.sum(t * weights, axis=0)
         sStack = norm * np.sqrt(np.sum(s2Full * weights**2, axis=0))
         return stack, sStack


      # or, if requested, compute and return the stacked cutout map
      else:
         # define chunks
         nChunk = self.nProc
         chunkSize = self.Catalog.nObj / nChunk
         # list of indices for each of the nChunk chunks
         chunkIndices = [range(iChunk*chunkSize, (iChunk+1)*chunkSize) for iChunk in range(nChunk)]
         # make sure not to miss the last few objects:
         # add them to the last chunk
         chunkIndices[-1] = range((nChunk-1)*chunkSize, self.Catalog.nObj)

         # select weights for a typical aperture size (not the smallest, not the largest)
         iRAp0 = self.nRAp / 2
         norm = norm[iRAp0]
         # need to link object number with weight,
         # despite the mask
         weightsLong = np.zeros(self.Catalog.nObj)
         weightsLong[mask] = weights[:,iRAp0]

         def stackChunk(iChunk):
            # object indices to be processed
            chunk = chunkIndices[iChunk]

            # start with a null map for stacking
            resMap = self.cutoutGeometry()
            for iObj in chunk:
               if iObj%10000==0:
                  print "- analyze object", iObj
               if self.overlapFlag[iObj]:
                  # Object coordinates
                  ra = self.Catalog.RA[iObj]   # in deg
                  dec = self.Catalog.DEC[iObj] # in deg
                  z = self.Catalog.Z[iObj]
                  # extract postage stamp around it
                  opos, stampMap, stampMask, stampHit = self.extractStamp(ra, dec, test=False)
                  resMap += stampMap * weightsLong[iObj]
            return resMap

         # dispatch each chunk of objects to a different processor
         with sharedmem.MapReduce(np=self.nProc) as pool:
            resMap = np.array(pool.map(stackChunk, range(nChunk)))

         # sum all the chunks
         resMap = np.sum(resMap, axis=0)
         # normalize by the proper sum of weights
         resMap *= norm
         return resMap


