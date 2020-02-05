

def computeBiasSnrMMax(self, filterType, est):
   '''Compute the significance of the bias as a function
   of the Mmax cut in the sample.
   '''


      # values of mMax
      nMMax = 20
      # choose the mMax values to have equal number of galaxies
      MMax = np.interp(np.linspace(0, self.Catalog.nObj, nMMax+1),
            np.arange(self.Catalog.nObj),
            np.sort(self.Catalog.Mvir))[1:]
      print nMMax, len(MMax)
      print MMax

      # bias from tsz to ksz for each mMax cut
      snrTszBias = np.zeros_like(MMax)

      for mMax in MMax:
         stack, sStack = self.computeStackedProfile(filterType, est, tTh='tsz')
         
         
