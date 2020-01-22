   
   def computeStackedProfile(self, filterType, est, iBootstrap=None, iVShuffle=None, tTh=None, stackedMap=False):
      """Returns the estimated profile and its uncertainty for each aperture.
      est: string to select the estimator
      iBootstrap: index for bootstrap resampling
      iVShuffle: index for shuffling velocities
      tTh: to replace measured temperatures by a theory expectation
      """
      # select objects that overlap, and reject point sources
      mask = self.catalogMask(overlap=True, psMask=True, filterType=filterType)
      
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
      tNoMean = t - np.mean(t, axis=0)

      # v/c [dimless]
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


         
         if stackedMap:










































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

