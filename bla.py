   ###############################################################################
   # optical depth relevant for kSZ

   # integrated optical depth to Thompson scattering
   # this is (total nb of electrons) * sigma_T / (a chi)^2
   # dimless
   # mVir is dark matter + baryon mass in M_sun
   # z redshift
   def fintegratedOpticalDepth(self, mVir, z):
      
      # convert from total mass to baryon mass
      # assume baryon mass = gas mass
      #tau = m * self.U.OmB/self.U.OmC
      tau = mVir * self.U.OmB/self.U.OmM
      
      # convert from total baryon mass to electron total number
      me = 9.10938291e-31  # electron mass (kg)
      mH = 1.67262178e-27  # proton mass (kg)
      mHe = 4.*mH # helium nucleus mass (kg)
      xH = 0.76   # primordial hydrogen fraction by mass
      nH_ne = 2.*xH/(xH+1.)
      nHe_ne = (1.-xH)/(2.*(1.+xH))
      msun = 1.989e30   # solar mass (kg)
      factor = (me + nH_ne*mH + nHe_ne*mHe) * (1./msun)   # in (Msun)
      tau /= factor
      #print "Ne=", tau
      
      # multiply by Thomson cross section (physical)
      mpc = 3.08567758e16*1.e6   # 1Mpc in m
      sigma_T = 6.6524e-29 # Thomson cross section in m^2
      tau *= sigma_T *(self.U.h/mpc)**2
      #print "sigmaT=", sigma_T *(self.U.h/mpc)**2
      
      # divide by (a chi)^2
      tau /= (self.U.ComovDist(1./(1.+z), 1.)/(1.+z))**2
      #print "(achi)^2=", (self.U.ComovDist(1./(1.+z), 1.)/(1.+z))**2
      
      return tau

#   # this is the expected signal, in muK * (Mpc/h)^2
#   # before division by the area of the filter
#   # assuming the cluster is entirely contained within the disk of the filter
#   # m is halo mass in Msun/h, ie mvir
#   def fkSZsignal(self, m, z, v):
#      c = 3.e5 # speed of light (km/s)
#      Tcmb = 2.726   # in Kelvin
#      result = self.fintegratedOpticalDepth(m, z) * v/c * Tcmb*1.e6
#      return result

   ###############################################################################
   # tSZ
   
   # frequency dependence for tSZ
   def ftSZFreqDpdce(self, nu):
      kB = 1.38e-23
      Tcmb = 2.725
      h = 6.63e-34
      x = h*nu/(kB*Tcmb)
      return x*(np.exp(x)+1.)/(np.exp(x)-1.) -4.
   
   
   # simple power-law fit to Greco et al 2014, fig4
   # giving the expected tSZ signal = int d^2theta delta(T), in muK*sr
   # before division by the area of the filter
   # as a function of mStellar
   def ftSZsignal(self, mStellar, z):
      
      # in arcmin^2
      yCcyltilda = (mStellar/1.e11)**3.2 * 1.e-6
      
#      # test the scaling
#      Mstellar = np.logspace(np.log10(1.e11), np.log10(1.e12), 1001, 10.)
#      YCcyltilda = (Mstellar/1.e11)**3.2 * 1.e-6
#      plt.figure(0)
#      plt.loglog(Mstellar, YCcyltilda)
#      plt.show()

      # in arcmin^2
      yCcyl = yCcyltilda * (self.U.Hubble(1./(1.+z))/self.U.Hubble(1.))**(2./3.)
      yCcyl /= (self.U.ComovDist(1./(1.+z), 1.)/(500.*self.U.h))**2
      
      # in muK * arcmin^2
      result = yCcyl * self.ftSZFreqDpdce(148.e9)
      Tcmb = 2.725
      result *= Tcmb*1.e6
      # in muK * sr
      result *= (np.pi/180./60.)**2
      return result
