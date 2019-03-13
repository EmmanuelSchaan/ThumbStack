      ###################################################################
      # compute significance

      # compute chi2_null
      chi2_0 = Alpha[imin:imax].dot( np.linalg.inv(Cov[imin:imax,imin:imax]).dot(Alpha[imin:imax]) )
      #print "chi^2_null", chi2_0

      # a is the gas fraction
      # sigma_cluster in arcmin
      def fchi2(p):
         a = p[0]
         result = (Alpha[imin:imax]-a*Theory[imin:imax]).dot( np.linalg.inv(Cov[imin:imax,imin:imax]).dot(Alpha[imin:imax]-a*Theory[imin:imax]) )
         result -= chi2_0
         return result

      p0 = 1.
      res = optimize.minimize(fchi2, p0)
      #print res
      abest = res.x[0]
      
      # number of dof of the fit
      #print "number of dof of fit =", len(Alpha[imin:imax])-1.

      # goodness of fit for null hypothesis
      #print "null chi20=", chi2_0
      pte0 = 1.- stats.chi2.cdf(chi2_0, len(Alpha[imin:imax])-1.)
      #print "null pte=", pte0
      # pte as a function of sigma, for a Gaussian random variable
      fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.)) - pte0
      sigma0 = optimize.brentq(fsigmaToPTE , 0., 50.)
      #print "null sigma significance=", sigma0

      # goodness of fit for best fit
      chi2_best = fchi2([abest, 2.])+chi2_0
      #print "best-fit chi2=", chi2_best
      pte = 1.- stats.chi2.cdf(chi2_best, len(Alpha[imin:imax])-1.)
      #print "best fit pte=", pte
      # pte as a function of sigma, for a Gaussian random variable
      fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.)) - pte
      sigma = optimize.brentq(fsigmaToPTE , 0., 50.)
      #print "pte sigma significance=", sigma

      # favour of best fit over null
      #print "best-fit chi2_0-chi2=", -fchi2([abest, 2.])
      sqrtDeltaChi2 = np.sqrt(abs(fchi2([abest, 2.])))
      #print "best-fit sqrt(delta chi2)=", sqrtDeltaChi2
