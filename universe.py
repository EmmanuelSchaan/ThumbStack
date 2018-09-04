from headers import *
##################################################################################


class Universe(object):

   def __init__(self, params=None):


      # all cosmological parameters
      if params is None:
         # neutrino masses
         self.Mnu = 0.06 # eV, minimum possible masses
         self.normalHierarchy = True
         self.nuMasses = self.computeNuMasses(self.Mnu, normal=self.normalHierarchy)

         self.params = {
                  'output': 'dTk vTk mPk',#'lCl tCl pCl mPk',
#                  'l_max_scalars': 2000,
#                  'lensing': 'yes',
                  'A_s': 2.3e-9,
                  'n_s': 0.9624,
                  'h': 0.6711,
                  'N_ur': 3.046,
                  'omega_b': 0.022068,
                  'Omega_cdm': 0.32,
                  'Omega_k': 0.,
                  'P_k_max_1/Mpc': 10.,
#                  'N_ncdm': 3,
#                  'm_ncdm': str(self.nuMasses[0])+','+str(self.nuMasses[1])+','+str(self.nuMasses[2]),
#                  'deg_ncdm': '1, 1, 1',
                  'non linear': 'halofit',
                  'z_max_pk': 100.
                  }
      else:
         self.params = params
      
      # run CLASS
      self.engine = CLASS.ClassEngine(self.params)
      self.bg = CLASS.Background(self.engine)
      self.sp = CLASS.Spectra(self.engine)
      self.th = CLASS.Thermo(self.engine)
      self.pm = CLASS.Primordial(self.engine)


      # wave vectors computed for power spectrum (h/Mpc)
      self.kMin = self.sp.P_k_min
      self.kMax = self.sp.P_k_max
      self.nK = 1001
      self.K = np.logspace(np.log10(self.kMin), np.log10(self.kMax), 10001, 10.)


      # physical constants
      self.G = 6.67e-11   # Newton's constant in SI
      self.c_kms = 299792458 / 1.e3 # light celerity in km/s
      self.rhoCRIT = 2.7617e11 # critical density today, in (h^-1 solarM) (h^-1 Mpc)^-3

      
      # convert from physical to comving density
      # and to (h^-1 solarM) (h^-1 Mpc)^-3
      self.rho_crit = lambda z: self.bg.rho_crit(z)/(1.+z)**3 * 1.e10
      self.rho_m = lambda z: self.bg.rho_m(z)/(1.+z)**3 * 1.e10
      self.rho_cdm = lambda z: self.bg.rho_cdm(z)/(1.+z)**3 * 1.e10
      self.rho_b = lambda z: self.bg.rho_b(z)/(1.+z)**3 * 1.e10
      self.rho_g = lambda z: self.bg.rho_g(z)/(1.+z)**3 * 1.e10
      self.rho_k = lambda z: self.bg.rho_k(z)/(1.+z)**3 * 1.e10
      self.rho_lambda = lambda z: self.bg.rho_lambda(z)/(1.+z)**3 * 1.e10
      self.rho_fld = lambda z: self.bg.rho_fld(z)/(1.+z)**3 * 1.e10
      self.rho_ncdm = lambda z: self.bg.rho_ncdm(z)/(1.+z)**3 * 1.e10
      self.rho_r = lambda z: self.bg.rho_r(z)/(1.+z)**3 * 1.e10
      self.rho_ur = lambda z: self.bg.rho_ur(z)/(1.+z)**3 * 1.e10
      self.rho_tot = lambda z: self.bg.rho_tot(z)/(1.+z)**3 * 1.e10

      # convert to (km/s)/(Mpc/h)
      self.hubble = lambda z: self.bg.hubble_function(z) * self.c_kms / self.bg.h
      
      # age of universe, in Gyr/h
      self.time = lambda z: self.bg.time(z) * self.bg.h



   def computeNuMasses(self, mSum, normal=True):
      '''mSum: sum of neutrino masses in eV
      normal=True for normal hierarchy
      output: masses in eV
      '''
      dmsq_atm = 2.5e-3 # eV^2
      dmsq_solar = 7.6e-5 # eV^2
      if normal:
         f = lambda m0: m0 + np.sqrt(m0**2+dmsq_solar) + np.sqrt(m0**2+dmsq_solar+dmsq_atm) - mSum
         m0 = optimize.brentq(f , 0., mSum)
         result = np.array([m0, np.sqrt(m0**2+dmsq_solar), np.sqrt(m0**2+dmsq_solar+dmsq_atm)])
      else:
         f = lambda m0: m0 + np.sqrt(m0**2+dmsq_atm) + np.sqrt(m0**2+dmsq_atm+dmsq_solar) - mSum
         m0 = optimize.brentq(f , 0., mSum)
         result = np.array([m0, np.sqrt(m0**2+dmsq_atm), np.sqrt(m0**2+dmsq_atm+dmsq_solar)])
      return result



   def fP2hinterp(self, k, z):
      '''This is actually Plin.
      Used for my halo model code
      '''
      if k<self.kMin or k>self.kMax:
         return 0.
      else:
         return self.sp.get_pklin(k, z)
   
   def fP1hinterp(self, k, z):
      '''Used for my halo model code
      '''
      return 0.
   
   def fPinterp(self, k, z):
      '''This is actually Pnl.
      Used for my halo model code
      '''
      if k<self.kMin or k>self.kMax:
         return 0.
      else:
         return self.sp.get_pk(k, z)


   ##################################################################################
   # spherical collapse
   ##################################################################################


   def deltaC_z(self, z):
      """critical density for spherical collapse at redshift z
      from Henry 2000, from Nakamura & Suto 1997
      usual 3.*(12.*pi)**(2./3.) / 20. = 1.686 if OmM=1.
      """
      x = ( 1./self.OmM - 1. )**(1./3.)
      x /= 1.+z
      dc = 3.*(12.*np.pi)**(2./3.) / 20.
      dc *= 1. - 0.0123* np.log( 1. + x**3 )
      return dc


   '''
   def DeltaVir(self, z):
      """Overdensity wrt mean for virialized halo at z
      from Bullock et al 2001, from Bryan & Norman 1998
      usual 18*pi**2 if OmM=1.
      gives 337 at z=0 for OmM= 0.3
      Omega = rhocrit(z)/rho_matter(z)
      """
      Omega = self.OmM*(1.+z)**3
      Omega /= self.OmM*(1.+z)**3 + (1. - self.OmM)
      x = Omega - 1.
      Dvir = 18*np.pi**2 + 82.*x - 39.*x**2
      # convert between rho_m(z) and rho_crit(z)
      Dvir /= Omega
      return Dvir
   '''
   '''
   def DeltaVir(self, z):
      """Overdensity wrt mean for virialized halo at z
      from Henry 2000, from Nakamura & Suto 1997
      usual 18*pi**2 if OmM=1.
      """
      x = ( 1./self.OmM - 1. )**(1./3.)
      x /= 1.+z
      Dvir = 18*np.pi**2 * ( 1. + 0.4093* x**2.71572 )
      return Dvir
   '''


   def Deltacrit_z(self, z):
      """ratio of virialized density to critical density at collapse (dimless).
      from Bullock et al 2001, from Bryan & Norman 1998
      usual 18*pi**2 if OmM=1.
      Omega = rhocrit(z)/rho_matter(z).
      """
      f = self.OmM * (1.+z)**3 / ( self.OmM * (1.+z)**3 + (1. - self.OmM) )
      return 18.*np.pi**2 + 82.*(f-1.) - 39.*(f-1.)**2


   def frvir(self, m, z):
      """Comoving virial and scale radii (Mpc/h)
      input mass is mvir (Msun/h)
      """
      Rvir = ( 3.*m / (4*np.pi*self.rhocrit_z(z)*self.Deltacrit_z(z)) )**(1./3.)  # in h^-1 Mpc
      return Rvir
   
   
   ##################################################################################

   def fK(self, chi):
      """Comoving transverse distance, as a function of the comoving radial distance chi.
      Input and output in Mpc/h
      """
      if self.bg.Omega0_k == 0.:
         result = chi
      # negative Om0_k corresponds to positive curvature
      elif self.bg.Omega0_k < 0.:
         Rk = 3.e3 / np.sqrt(-self.bg.Omega0_k) # (c/H0) / sqrt(|Omega0_k|), in Mpc/h
         result = Rk * np.sin(chi / Rk)
      # positive Om0_k corresponds to negative curvature
      elif self.bg.Omega0_k > 0.:
         Rk = 3.e3 / np.sqrt(self.bg.Omega0_k) # (c/H0) / sqrt(|Omega0_k|), in Mpc/h
         result = Rk * np.sinh(chi / Rk)
      return result

   def plotDistances(self, zMax=2.):
      z = np.linspace(0., zMax, 512)
      
      # comoving radial distance
      plt.plot(z, self.bg.comoving_distance(z), label=r"$\chi$")
      # comoving angular diameter distance
      plt.plot(z, self.bg.comoving_transverse_distance(z), ls='-.', label=r"$d_A$")
      # comoving angular diameter distance 2
      plt.plot(z, self.bg.angular_diameter_distance(z)*(1.+z), ls='--', label=r"$d_A=D_A/a$")
      # comoving angular diameter distance 2
      plt.plot(z, self.fK(self.bg.comoving_distance(z)), '.', label=r"$d_A=f_K(\chi)$")
      # luminosity distance
      plt.plot(z, self.bg.luminosity_distance(z), label=r"$d_L$")
      # physical angular diameter distance
      plt.plot(z, self.bg.angular_diameter_distance(z), label=r"$D_A$")
      plt.legend()
      plt.xlabel(r"$z$")
      plt.ylabel(r"distance $[h^{-1} \mathrm{Mpc}]$")
      plt.show()
   
      # distance modulus
      DistMod = 5.*( np.log10(self.bg.luminosity_distance(z)*1.e6*self.bg.h) - 1. )
      plt.figure(1)
      ax=plt.subplot(111)
      ax.plot(z, DistMod)
      ax.grid()
      ax.set_title(r'Distance modulus $\mu=5 (log_{10}(d/pc) - 1)$')
      ax.set_xlabel('redshift z')
      ax.set_ylabel('distance modulus')
      plt.show()
   
   
   def printCosmoParams(self):
      print "h = ", self.bg.h
      print "Omega0_m = ", self.bg.Omega0_m
      print "Omega0_lambda = ", self.bg.Omega0_lambda
      print "Omega0_r = ", self.bg.Omega0_r
      print "Omega0_k = ", self.bg.Omega0_k

   def plotDensity(self):
      z = np.linspace(0., 300., 512)
      plt.loglog(1.+z, self.rho_m(z), label=r"$\rho_m$")
      plt.loglog(1.+z, self.rho_cdm(z), label=r"$\rho_{cdm}$")
      plt.loglog(1.+z, self.rho_b(z), label=r"$\rho_{b}$")
      plt.loglog(1.+z, self.rho_r(z), label=r"$\rho_r$")
      plt.loglog(1.+z, self.rho_g(z), label=r"$\rho_g$")
      plt.loglog(1.+z, self.rho_ur(z), label=r"$\rho_{ur}$")
      #plt.loglog(1.+z, self.rho_ncdm(z), label=r"$\rho_{ncdm}$")
      plt.loglog(1.+z, self.rho_lambda(z), label=r"$\rho_\Lambda$")
      plt.loglog(1.+z, self.rho_fld(z), label=r"$\rho_\text{fluid}$")
      plt.legend()
      plt.xlabel(r"$1.+z$")
      plt.ylabel(r"Comoving density [(M$_\odot$/h)/(Mpc/h)$^3$]")
      plt.show()
   
   def plotDensityParameters(self):
      z = np.linspace(0., 300., 512)
      plt.loglog(1.+z, self.bg.Omega_m(z), label=r"$\Omega_m$")
      plt.loglog(1.+z, self.bg.Omega_r(z), label=r"$\Omega_r$")
      plt.loglog(1.+z, self.bg.Omega_lambda(z), label=r"$\Omega_\Lambda$")
      plt.loglog(1.+z, self.bg.Omega_fld(z), label=r"$\Omega_\text{fluid}$")
      plt.legend()
      plt.xlabel(r"$1.+z$")
      plt.ylabel("Energy density parameters")
      plt.show()

   def plotP(self, z=0.):
      plt.loglog(self.K, self.sp.get_pk(k=self.K, z=z), label='nonlinear')
      plt.loglog(self.K, self.sp.get_pklin(k=self.K, z=z), label='linear')
      plt.legend()
      plt.xlabel(r"$k$ $[h / \mathrm{Mpc}]$")
      plt.ylabel(r"$P$ $[(\mathrm{Mpc}/h)^3]$")
      plt.show()

   def plotSigma8(self, zMax=2.):
      z = np.linspace(0., zMax, 512)
      plt.plot(1. + z, self.sp.sigma8_z(z))
      plt.xlabel(r"$1+z$")
      plt.ylabel(r"$\sigma_8(z)$")
      plt.show()

   def plotTransfer(self):
      transfer = self.sp.get_transfer(z=0)
      print(transfer.dtype.names)
      plt.subplot(211)
      plt.plot(transfer['k'], transfer['d_tot'])
      plt.ylabel("total density transfer")
      plt.subplot(212)
      plt.plot(transfer['k'], transfer['t_tot'])
      plt.xlabel(r"$k$ $[h\mathrm{Mpc}^{-1}]$")
      plt.ylabel("total velocity transfer")
      plt.show()

   def plotHubble(self, zMax=2.):
      z = np.linspace(0., zMax, 512)
      plt.plot(1. + z, self.hubble(z) * self.bg.h)
      plt.xlabel(r"$1+z$")
      plt.ylabel(r"$H(z)$ [km/s/Mpc]")
      plt.show()

   def plotTime(self, zMax=2.):
      z = np.linspace(0., zMax, 512)
      plt.plot(1. + z, self.time(z) / self.bg.h, '-')
      plt.plot(1. + z, self.bg.time(z), '--')
      plt.xlabel(r"$1+z$")
      plt.ylabel(r"$t(z)$ [Gyr]")
      plt.show()

   def plotLinearGrowth(self, zMax=10.):
      z = np.linspace(0., zMax, 512)
      plt.plot(1. + z, self.bg.scale_independent_growth_factor(z), '-', label=r'$D(a)$')
      plt.plot(1. + z, 1./(1.+z), '--', label=r'$a(z)$')
      plt.legend(loc=1)
      plt.xscale('log')
      plt.yscale('log')
      plt.xlabel(r"$1+z$")
      plt.ylabel(r"Linear growth rate")
      plt.show()

   def plotThermo(self):
      # recombination
      print "recombination redshift =", self.th.z_rec
      print "sound horizon at recombination =", self.th.rs_rec, "Mpc/h"
      print "sound horizon angle at recombination =", self.th.theta_s * 180./np.pi, "deg"
      # drag
      print "drag redshift =", self.th.z_drag
      print "sound horizon at z_drag =", self.th.rs_drag, "Mpc/h"
      # reionization
      print "reionization redshift =", self.th.z_reio
      print "reionization optical depth =", self.th.tau_reio
   
   def plotPrimordial(self):
      plt.loglog(self.K, self.pm.get_pkprim(self.K), label='from CLASS')
      plt.loglog(self.K, self.sp.A_s * (self.K / self.sp.k_pivot)**(self.sp.n_s-1.), ls='--', label='analytic')
      plt.legend()
      plt.xlabel(r"$k$ $[h\mathrm{Mpc}^{-1}]$")
      plt.ylabel(r"dimless prim. power $\Delta_\mathcal{R}(k)$")
      plt.show()
   


   ##################################################################################
   # density perturbations
   ##################################################################################


   def dlnPlindlnK(self, k, z):
      """derivative of linear power spectrum wrt k
      """
      e = 0.01
      kup = k*(1.+e)
      kdown = k*(1.-e)
      if kup>self.kMax or kdown<self.kMin:
         result = 0.
      else:
         result = self.sp.get_pklin(kup, z) / self.sp.get_pklin(kdown, z)
         result = np.log(result) / (2.*e)
      return result

   def dlnPnldlnK(self, k, z):
      """derivative of halofit power spectrum wrt k
      """
      e = 0.01
      lup = k*(1.+e)
      ldown = k*(1.-e)
      if kup>self.kMax or kdown<self.kMin:
         result = 0.
      else:
         result = self.sp.get_pk(kup, z) / self.sp.get_pk(kdown, z)
         result = np.log(result) / (2.*e)
      return result


   def Sigma2(self, R, z, W3d):
      """variance of delta on an isotropic 3d domain,
      defined by W3d
      R in h^-1 Mpc, comoving scale, output is dimless
      """
      f = lambda lnk: np.exp(lnk)**3 * self.sp.get_pklin(np.exp(lnk), z) * W3d(np.exp(lnk)*R)**2 / (2* np.pi**2) # dimensionless
      result = integrate.quad(f, np.log(self.kMin), np.log(self.kMax), epsabs=0., epsrel=1.e-3)[0]
      return result


   def nonLinMass(self, z):
      """nonlin mass at z, in Msun/h
      from Takada and Jain 2002/2003
      """
      # bounds for looking for m_nonlinin, in (h^-1 solarM) (h^-1 Mpc)^-3
      ma = 1.e11
      mb = 1.e13
      # solve for sigma(m)^2 = deltaC**2
      R = lambda m: ( 3.* m / (4*np.pi*self.rho_m(z)) )**(1./3.)   # in h^-1 Mpc
      f = lambda m: self.Sigma2(R(m), z, W3d_sth) - self.deltaC_z(z)**2
      # find mass such that nu(m, z) = 1
      result = optimize.brentq(f , ma, mb)
      return result


   def Sigma2_2d(self, r, z, W2d):
      """this is Dchi * sigma^2,
      where sigma^2 = <( delta averaged on bin Dchi )^2>
      r is comoving scale in h^-1 Mpc
      """
      f = lambda k: k * self.sp.get_pklin(k, z) * W2d(k*r)**2 / (2* np.pi)
      result = integrate.quad(f, self.kMin, self.kMax, epsabs=0., epsrel=1.e-3)[0]
      return result


   def Sigma2_cone(self, aMin, aMax, theta, W2d):
      """Variance of delta on a section of a cone,
      with fixed angular radius theta
      from aMin to aMax
      output dimless
      """
      z = lambda a: 1./a-1.
      r = lambda a: self.bg.comoving_distance(z(a)) * theta
      integrand = lambda a: 3.e5/(self.hubble(a) * a**2) * self.Sigma2_2d(r(a), z(a), W2d)
      result = integrate.quad(integrand, aMin, aMax, epsabs=0., epsrel=1.e-3)[0]
      chiMin = self.comoving_distance(1./aMax-1.)
      chiMax = self.comoving_distance(1./aMin-1.)
      result /= (chiMax-chiMin)**2
      return result


   def dlnSigma2_dlnR(self, R, z):
      """R in h^-1 Mpc, comoving scale, output is dimless
      dln(sigma2) / dln(R)
      """
      f = lambda lnk: np.exp(lnk)**3 * self.sp.get_pklin(np.exp(lnk),z) * 2. * W3d_sth(np.exp(lnk)*R) * dW3d_sth(np.exp(lnk)*R) * np.exp(lnk)*R / (2* np.pi**2)  # dimensionless
      result = integrate.quad(f, np.log(self.kMin), np.log(self.kMax), epsabs=0., epsrel=1.e-3)[0]
      result /= self.Sigma2(R, z, W3d_sth)
      return result


   def fdlnSigma_dlnM(self, m, z):
      """dln(sigma)/dln(m)
      """
      R = (3.*m / (4.*np.pi*self.rho_m(z)))**(1./3.)
      result = self.dlnSigma2_dlnR(R, z) /6.
      return result



   def fnu(self, m, z):
      """nu = dc**2/sigma2(m, z)
      """
      r = (3.*m / (4.*np.pi*self.rho_m(z)))**(1./3.)
      s2 = self.Sigma2(r, z, W3d_sth)
      nu = self.deltaC_z(z)**2 / s2
      return nu


   def fdlnnu_dlnm(self, m, z):
      """dln(nu)/dln(m)
      """
      return -2.* self.fdlnSigma_dlnM(m, z)


   ##################################################################################
   # halo mass conversion


   def massRadiusConversion(self, m, z, value=200, ref="m"):
      """converts virial mass (Msun/h)
      into M_value,ref (Msun/h),
      and corresponding comoving R_value,ref (Mpc/h).
      Assumes concentration from Duffy et al 2008.
      Used by Tinker mass function
      """

      # concentration params from Duffy et al 2008, used by Tinker
      cNFW0 = 5.71
      cNFWam = -0.084
      cNFWaz = -0.47

      # from Duffy et al 2008: different pivot mass
      cNFW = cNFW0 * (m/2.e12)**cNFWam * (1.+z)**cNFWaz
      # comoving virial radius and scale radius in h^-1 Mpc
      Rvir = ( 3.*m / (4*np.pi*self.rho_crit(z) * self.Deltacrit_z(z)) )**(1./3.)
      Rs = Rvir / cNFW
      # NFW scale density (comoving)
      rhoS = m / (4.*np.pi*Rs**3) / (np.log(1.+cNFW) - cNFW/(1.+cNFW))

      # comoving reference density
      if ref=="m":   # ie wrt mean density
         rhoRef = self.rho_m(z)
      elif ref=="c": # ie wrt critical density
         rhoRef = self.rho_crit(z)

      # get R200 and M200
      f = lambda x: -1. + 1./(1.+x) + np.log(1.+x) - value/3.*(rhoRef/rhoS)*x**3
      x = optimize.brentq(f , 0.1, 100.)
      Rnew = x * Rs
      Mnew = 4./3.*np.pi*Rnew**3 * rhoRef * value

      return Mnew, Rnew


   ##################################################################################
   # response of power spectrum to local overdensity

   def dlnPlindDelta(self, k, z):
      result = 68./21.
      result -= self.dlnPlindlnK(k, z) / 3.
      result -= 1.   # this is -1/3 * dlnk^3/dlnk
      return result



   def plotdlnPlindDelta(self, z=0.):
      """for the linear power spectrum, dlnP/ddelta is independent of z
      but dP/ddelta is not
      """
      #
      f = lambda k: self.sp.get_pklin(k, z)
      P = np.array(map(f, self.K))
      #
      f = lambda k: self.dlnPlindDelta(k, z)
      dP = np.array(map(f, self.K))

      # dlnP/ddelta
      fig = plt.figure(0)
      ax = plt.subplot(111)
      #
      ax.semilogx(self.K, dP, 'b', label=r'$\frac{68}{21} - \frac{1}{3}\frac{d\ln k^3P_\text{lin}}{d\ln k}$')
      ax.axhline(68./21., color='r', label=r'$\frac{68}{21}$')
      #
      ax.legend(loc=4)
      ax.grid()
      #ax.set_xlim((0., 0.4))
      #ax.set_ylim((2.2, 3.2))
      ax.set_xlabel(r'k [h/Mpc]')
      ax.set_ylabel(r'$d\ln P_\text{lin}/d\delta$ [(Mpc/h)$^3$]')
      #fig.savefig('./figures/response_powerspectrum/dLnPlindDelta.pdf', bbox_inches='tight')

      # dP/ddelta
      fig = plt.figure(1)
      ax = plt.subplot(111)
      #
      ax.loglog(self.K, P*dP, 'b', label=r'$P(k)\left[\frac{68}{21} - \frac{1}{3}\frac{d\ln k^3P_\text{lin}}{d\ln k} \right]$')
      ax.loglog(self.K, P*68./21., 'r', label=r'$P(k)\left[\frac{68}{21}\right]$')
      ax.loglog(self.K, P, 'k', label=r'$P(k)$')
      #
      ax.legend(loc=3)
      ax.grid()
      #ax.set_xlim((0., 0.4))
      #ax.set_ylim((2.2, 3.2))
      ax.set_xlabel(r'k [h/Mpc]')
      ax.set_ylabel(r'$dP_\text{lin}/d\delta$ [(Mpc/h)$^3$]')
      #fig.savefig('./figures/response_powerspectrum/dPlindDelta_loglog.pdf', bbox_inches='tight')

      plt.show()


   ##################################################################################
   # Velocity fluctuations

   def RMSVelocity(self, R, z, W3d):
      """Computes |v|_RMS in km/s.
      Input R in Mpc/h comoving.
      Assumes linear (Zel'dovich) relation between velocity and density,
      and uses the linear matter power spectrum.
      """
      def integrand(lnk):
         k = np.exp(lnk)
         result = k**3 / (2* np.pi**2) # d^3k/(2pi)^3 = dlnk*k^3/(2 pi^2) [(h/Mpc)^3]
         result *= np.abs(W3d(k*R))**2 # window function [dimless]
         result *= self.sp.get_pklin(k, z)  / k**2 # velocity power spectrum [(Mpc/h)^5]
         result *= self.bg.Omega_m(z)**(2.*5./9.) # f**2, with f = Omega_m(z)**5/9
         result *= (self.hubble(z) / (1.+z))**2 # (a*H(a))**2 [(km/s/(Mpc/h))^2]
         return result
      result = integrate.quad(integrand, np.log(self.kMin), np.log(self.kMax), epsabs=0., epsrel=1.e-3)[0]
      result = np.sqrt(result)
      return result


#   def RMSErrorRecVel(self, R, z, W3d):
#      """R in h^-1 Mpc, comoving scale, output is the rms error on reconstructed velocity in (km/s)
#      """
#      CorrCoeff = np.array(map(self.fr, self.K))
#      F = (self.K**3) * ( np.array(map(W3d, self.K * R))**2 ) / (2* np.pi**2)  # dimensionless
#      F *= self.Plin_z(z) / self.K**2
#      F *= 2.* (1. - CorrCoeff)
#      F *= ( self.OmM*(1.+z)**3 / (self.Hubble(1./(1.+z))/self.Hubble(1.))**2 )**(2.*5./9.)   # f**2, where f = Omega_m(z)**5/9
#      F *= ( self.Hubble(1./(1.+z))/(1.+z) )**2
#      dlnK = ( (self.K[1:]-self.K[:-1]) / self.K[:-1] )   # dimensionless
#      Itrap = np.sum( dlnK * ( F[:-1] + F[1:] ) ) * 0.5
#      return np.sqrt(Itrap)
#
#   def growthLogDerivativeF(self, z):
#      return ( self.OmM*(1.+z)**3 / (self.OmM*(1.+z)**3+(1. - self.OmM)) )**(5./9.)
#
#   def convertVelToDisp(self, z):
#      """multiply a velocity at redshift z by this factor to get the displacement at that redshift
#      v in km/s, displacement in Mpc/h
#      """
#      factor = self.Hubble(z)/(1.+z)*self.growthLogDerivativeF(z)
#      return 1./factor
#
#   def RMSDisplacement(self, R, z, W3d):
#      """R in h^-1 Mpc, comoving scale, output is the rms displacement in (Mpc/h)
#      """
#      result = self.RMSVelocity(R, z, W3d)
#      result *= self.convertVelToDisp(z)
#      return result
#
##   # variance of vRadial_RMS^2, the average vRadial^2 on a spherical volume
##   # with radius R in comoving Mpc/h
##   def varRMSVRadialSpherical(self, R, z, W3d):
##
##      def integrand(pars):
##         k = pars[0]
##         K = pars[1]
##         mu = pars[2]
##         x = K/k
##         # Put the density instead of velocity!!!
##         result = self.fPlin_z(k, z) / k**2
##         result *= self.fPlin_z(k*np.sqrt(1 - 2.*x*mu + x**2), z) / (k*np.sqrt(1 - 2.*x*mu + x**2))**2
##
##         factor = ( self.OmM*(1.+z)**3 / (self.Hubble(1./(1.+z))/self.Hubble(1.))**2 )**(2.*5./9.)   # f**2, where f = Omega_m(z)**5/9
##         factor *= ( self.Hubble(1./(1.+z))/(1.+z) )**2 # (a*H)**2
##
##         result *= factor**2
##         result *= W3d(k*x*R)**2
##         result *= x**2 * k**5
##         result *= k
##         result *= 4./ 3. / (2.*np.pi)**4
##         return result
##
##      integ = vegas.Integrator([[1.e-4, 1.e1], [1.e-4, 1.e1], [-1., 1.]])
##      integ(integrand, nitn=4, neval=1000)
##      result = integ(integrand, nitn=8, neval=1000)
##      print result.sdev/result.mean
##      return result.mean
###      return result.summary()
#
#
#   def varRMSVSpherical(self, R, z, W3d):
#      """variance of v_RMS^2, the average v^2 on a spherical volume
#      with radius R in comoving Mpc/h
#      """
#
#      def integrand(pars):
#         k = pars[0]
#         K = pars[1]
#         mu = pars[2]
##         x = K/k
#         # Put the density instead of velocity!!!
#         result = self.fPlin_z(k, z) / k**2
#         result *= self.fPlin_z(np.sqrt(k**2 - 2.*k*K*mu + K**2), z) / (k**2 - 2.*k*K*mu + K**2)
#
#         factor = ( self.OmM*(1.+z)**3 / (self.Hubble(1./(1.+z))/self.Hubble(1.))**2 )**(2.*5./9.)   # f**2, where f = Omega_m(z)**5/9
#         factor *= ( self.Hubble(1./(1.+z))/(1.+z) )**2 # (a*H)**2
#
#         result *= factor**2
#         result *= W3d(K*R)**2
#         result *= K**2 * k**2
#         result *= 4. / (2.*np.pi)**4
#         return result
#
#      integ = vegas.Integrator([[1.e-4, 1.e1], [1.e-4, 1.e1], [-1., 1.]])
#      integ(integrand, nitn=8, neval=1000)
#      result = integ(integrand, nitn=8, neval=1000)
#      print result.sdev/result.mean
##      print result.summary()
#      return result.mean
#
#
#   def plotVarRMSVSpherical(self, z=0.):
#      vMin = 1.e-2 * 1.e9  # (Mpc/h)^3
#      vMax = 1.e1 * 1.e9
#      V = np.logspace(np.log10(vMin), np.log10(vMax), 11, 10.)
#      R = ( 3.*V/(4.*np.pi) )**(1./3.)
#
#      # We want the v_RMS^2, not the square of the average v,
#      # this is why the scale is 0 Mpc/h
#      f = lambda r: self.RMSVelocity(r*0., z, W3d_sth)**2
#      v2 = np.array(map(f, R))
#      # We want the cosmic variance on v_RMS^2;
#      # now this quantity will depend on the scale r
#      f = lambda r: self.varRMSVSpherical(r, z, W3d_sth)
#      varV2 = np.array(map(f, R))
#
#      # volume of D56
#      volumeD56 = 700.*(np.pi/180.)**2 # area in sr
#      volumeD56 *= self.ComovDist(1./(1.+0.57), 1.)**2 # area in (Mpc/h)^2
#      volumeD56 *= self.ComovDist(1./(1.+0.7), 1./(1.+0.4)) # volume in (Mpc/h)^3
#      # volume of D56+BOSS N
#      volumeD56BN = volumeD56 * 2700./700.
#
#      fig=plt.figure(0)
#      ax=fig.add_subplot(111)
#      #
#      ax.plot(V / 1.e9, np.sqrt(varV2)/v2, 'b-')
#      #
#      ax.axvline(volumeD56 / 1.e9, c='k', linestyle='--')
#      ax.axvline(volumeD56BN / 1.e9, c='k', linestyle='--')
#      #
#      ax.set_xscale('log')
#      ax.set_xlabel(r'Volume [Gpc/h]$^3$')
#      ax.set_ylabel(r'$\sigma_{v_\text{RMS}^2} / v_\text{RMS}^2$')
#      #
##      fig.savefig("./figures/velocities/variance_of_vrms2_spherical.pdf", bbox_inches='tight')
#
#      fig=plt.figure(1)
#      ax=fig.add_subplot(111)
#      #
#      ax.plot(V / 1.e9, np.sqrt(v2), 'b-')
#      ax.plot(V / 1.e9, np.sqrt(np.sqrt(varV2)), 'b-')
#      #
#      ax.set_xscale('log')
#      ax.set_xlabel(r'Volume [Gpc/h]$^3$')
#      ax.set_ylabel(r'$v_\text{RMS}$')
#
#
#      plt.show()
#      return


   ##################################################################################

   def plotHSV_2d_scaling(self):
      # range of redshift
      Z = np.linspace(0.005, 4., 8.)
      Z = np.array([0.005, 0.01, 0.02, 0.1, 0.2, 1.,  2., 3., 4.])
      A = 1./(1.+Z)
      # values of fsky to compute
      Fsky = np.logspace(np.log10(1.e-2), np.log10(1.), 51, 10.)

      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      for a in A:
         fscaling = lambda fsky: self.Sigma2_2d(self.bg.comoving_distance(1./a-1.) * 2.*np.sqrt(fsky), 1./a-1., W2d_cth) * fsky
         HSV = np.array(map(fscaling, Fsky))
         ax.loglog(Fsky, HSV, label=r'$z=$'+str(round(1./a-1., 2)))
      #
      ax.legend(loc=1)
      ax.set_xlabel(r'$f_{sky}$', fontsize=18)
      ax.set_ylabel(r'$\sigma^2(z, f_{sky})  f_{sky}$', fontsize=18)
      #fig.savefig('./figures/hsv_2d_scaling.pdf')

      plt.show()



   def plotHSV_3d_scaling(self):
      # survey volumes
      Vmin = 1.e-3 * (1.e3)**3
      Vmax = 10. * (1.e3)**3
      NV = 51
      V = np.logspace(np.log10(Vmin), np.log10(Vmax), NV, 10.)

      # variance of the matter overdensity smoothed over the survey (only used for HSV for 3d cov)
      R = (3.*V/(4.*np.pi))**(1./3.)   # for a spherical survey
      HSV = np.array( map(lambda r: self.Sigma2(r, 0., W3d_sth), R) ) * V

      fig=plt.figure(0)
      ax=plt.subplot(111)
      #
      ax.loglog(V/(1.e3)**3, HSV)
      #ax.loglog(V/(1.e3)**3, V**(-1./3.)* 1.e7)
      #ax.loglog(V/(1.e3)**3, V)
      #
      ax.set_xlabel(r'survey volume [(Gpc/h)$^3$]', fontsize=18)
      ax.set_ylabel(r'$\sigma^2 V$', fontsize=18)
      #fig.savefig('./figures/hsv_3d_scaling.pdf')

      plt.show()



##################################################################################

class UnivHillPajer13(Universe):

   def __init__(self):
      params = {
               'output': 'dTk vTk mPk',#'lCl tCl pCl mPk',
#                  'l_max_scalars': 2000,
#                  'lensing': 'yes',
               'A_s': 2.3e-9,
               'n_s': 0.9624,
               'h': 0.697,
               'N_ur': 3.046,
               'Omega_b': 0.0461,
               'Omega_cdm': 0.236,
               'Omega_k': 0.,
               'P_k_max_1/Mpc': 10.,
               'non linear': 'halofit',
               'z_max_pk': 100.
               }
      super(UnivHillPajer13, self).__init__()


##################################################################################

class UnivTinkerEtAl08(Universe):

   def __init__(self):
      params = {
               'output': 'dTk vTk mPk',#'lCl tCl pCl mPk',
#                  'l_max_scalars': 2000,
#                  'lensing': 'yes',
               'A_s': 2.3e-9,
               'n_s': 0.9624,
               'h': 0.7,
               'N_ur': 3.046,
               'Omega_b': 0.04,
               'Omega_cdm': 0.26,
               'Omega_k': 0.,
               'P_k_max_1/Mpc': 10.,
               'non linear': 'halofit',
               'z_max_pk': 100.
               }
      super(UnivTinkerEtAl08, self).__init__(params=params)


##################################################################################

class UnivSchaanEtAl14(Universe):

   def __init__(self):
      params = {
               'output': 'dTk vTk mPk',#'lCl tCl pCl mPk',
#                  'l_max_scalars': 2000,
#                  'lensing': 'yes',
               'A_s': 2.3e-9,
               'n_s': 0.9624,
               'h': 0.732,
               'N_ur': 3.046,
               'Omega_b': 0.042,
               'Omega_cdm': 0.196,
               'Omega_k': 0.,
               'P_k_max_1/Mpc': 10.,
               'non linear': 'halofit',
               'z_max_pk': 100.
               }
      super(UnivSchaanEtAl14, self).__init__(params=params)


##################################################################################

# uses WMAP9 + eCMB parameters
class UnivHandEtAl13(Universe):

   def __init__(self):
      params = {
               'output': 'dTk vTk mPk',#'lCl tCl pCl mPk',
#                  'l_max_scalars': 2000,
#                  'lensing': 'yes',
               'A_s': 2.3e-9,
               'n_s': 0.9624,
               'h': 0.705,
               'N_ur': 3.046,
               'Omega_b': 0.0449,
               'Omega_cdm': 0.227,
               'Omega_k': 0.,
               'P_k_max_1/Mpc': 10.,
               'non linear': 'halofit',
               'z_max_pk': 100.
               }
      super(UnivHandEtAl13, self).__init__(params=params)


##################################################################################

# uses the OmM and h that Mariana used for her reconstructed velocities
class UnivMariana(Universe):

   def __init__(self):
      params = {
               'output': 'dTk vTk mPk',#'lCl tCl pCl mPk',
#                  'l_max_scalars': 2000,
#                  'lensing': 'yes',
               'A_s': 2.3e-9,
               'n_s': 0.9624,
               'h': 0.7,
               'N_ur': 3.046,
               'Omega_b': 0.0458571,
               'Omega_cdm': 0.2441429,
               'Omega_k': 0.,
               'P_k_max_1/Mpc': 10.,
               'non linear': 'halofit',
               'z_max_pk': 100.
               }
      super(UnivMariana, self).__init__(params=params)



##################################################################################

class UnivPlanck15(Universe):

   def __init__(self):
      params = {
               'output': 'dTk vTk mPk',#'lCl tCl pCl mPk',
#                  'l_max_scalars': 2000,
#                  'lensing': 'yes',
               'A_s': 2.3e-9,
               'n_s': 0.9624,
               'h': 0.6712,
               'N_ur': 3.046,
               'Omega_b': 0.0493,
               'Omega_cdm': 0.267,
               'Omega_k': 0.,
               'P_k_max_1/Mpc': 10.,
               'non linear': 'halofit',
               'z_max_pk': 100.
               }
      super(UnivPlanck15, self).__init__(params=params)


##################################################################################

class UnivNuWCurv(Universe):

   def __init__(self):

      # neutrino masses
      self.Mnu = 0.06 # eV, minimum possible masses
      self.normalHierarchy = True
      self.nuMasses = self.computeNuMasses(self.Mnu, normal=self.normalHierarchy)

      params = {
               # Cosmological parameters
               'Omega_cdm': 0.267,
               'Omega_b': 0.0493,
               'A_s': 2.3e-9,
               'n_s': 0.9624,
               'tau_reio': 0.06,
               'h': 0.6712,
               #
               # Massive neutrinos
               'N_ncdm': 3,
               'm_ncdm': str(self.nuMasses[0])+','+str(self.nuMasses[1])+','+str(self.nuMasses[2]),
               'deg_ncdm': '1, 1, 1',
               #
               # w0 and wa
               'Omega_Lambda': 0.,
               'w0_fld': -1.,
               'wa_fld': 0.,
               #
               # Curvature
               'Omega_k': 0.,
               #
               # parameters
               'output': 'mPk dTk vTk',#'lCl tCl pCl mPk',
               'P_k_max_1/Mpc': 10.,
               'non linear': 'halofit',
               'z_max_pk': 100.
               }
      super(UnivNuWCurv, self).__init__(params=params)

