from headers import *

##################################################################################
##################################################################################

class Catalog(object):

   def __init__(self, U, MassConversion, name="test", nameLong=None, pathInCatalog="", save=False):

      self.U = U
      self.MassConversion = MassConversion
      self.name = name
      if nameLong is None:
         self.nameLong = self.name
      else:
         self.nameLong = nameLong
      self.pathInCatalog = pathInCatalog
      
      # Output path
      self.pathOut = "./output/catalog/"+self.name
      if not os.path.exists(self.pathOut):
         os.makedirs(self.pathOut)
      # catalog path
      self.pathOutCatalog = self.pathOut + "/catalog.txt"
      
      # Figures path
      self.pathFig = "./figures/catalog/"+self.name
      if not os.path.exists(self.pathFig):
         os.makedirs(self.pathFig)
      
      
      if save:
         self.readInputCatalog()
         self.addHaloMass()
         self.writeCatalog()
      
      self.loadCatalog()
   

   ##################################################################################
   ##################################################################################

   def copy(self, name="test", nameLong=None):
      """Copy a catalog class, with the option of changing the name.
      """
      # First copy the output catalog
      # new catalog path
      newPathOut = "./output/catalog/"+name
      if not os.path.exists(newPathOut):
         os.makedirs(newPathOut)
      newPathOutCatalog = newPathOut + "/catalog.txt"
      # copy the output catalog
      copyfile(self.pathOutCatalog, newPathOutCatalog)
      
      # Then copy the catalog properties
      newCat = Catalog(self.U, self.MassConversion, name=name, nameLong=nameLong, pathInCatalog=self.pathInCatalog, save=False)
      return newCat


   ##################################################################################
   ##################################################################################

   def readInputCatalog(self):
      print "- read input catalog from "+self.pathInCatalog
      data = np.genfromtxt(self.pathInCatalog)
      self.nObj = len(data[:,0])
      #
      # sky coordinates and redshift
      self.RA = data[:,0] # [deg]
      self.DEC = data[:,1]   # [deg]
      self.Z = data[:,2]
      #
      # observed cartesian coordinates
      self.coordX = data[:,3]   # [Mpc/h]
      self.coordY = data[:,4]   # [Mpc/h]
      self.coordZ = data[:,5]   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      self.dX = data[:,6]   # [Mpc/h]
      self.dY = data[:,7]   # [Mpc/h]
      self.dZ = data[:,8]   # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      self.dXKaiser = data[:,9]   # [Mpc/h] from cartesian catalog difference
      self.dYKaiser = data[:,10]   # [Mpc/h]
      self.dZKaiser = data[:,11]   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      self.vX = data[:,12]   #[km/s]
      self.vY = data[:,13]   #[km/s]
      self.vZ = data[:,14]   #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      self.vR = data[:,15]  # [km/s]   from spherical catalo
      self.vTheta = data[:,16]   # [km/s]
      self.vPhi = data[:,17]  # [km/s]
      #
      # Stellar masses
      self.Mstellar = data[:,18]   # [M_sun], from Maraston et al


   ##################################################################################

   def addHaloMass(self):
      """Generate halo masses in M_sun,
      from stellar masses in M_sun.
      """
      print "- add halo masses"
      # flag: 1 if object has mass
      self.hasM = np.zeros(self.nObj)
      
      self.Mvir = np.zeros(self.nObj)
      for iObj in range(self.nObj):
         mStellar = self.Mstellar[iObj]
         if (mStellar>1.e3) and not np.isnan(mStellar):
            self.hasM[iObj] = True
            self.Mvir[iObj] = self.MassConversion.fmStarTomVir(mStellar)

#   def addOptical


   ##################################################################################

   def writeCatalog(self):
      print "- write full catalog to "+self.pathOutCatalog
      data = np.zeros((self.nObj,21))
      #
      # sky coordinates and redshift
      data[:,0] = self.RA # [deg]
      data[:,1] = self.DEC   # [deg]
      data[:,2] = self.Z
      #
      # observed cartesian coordinates
      data[:,3] = self.coordX   # [Mpc/h]
      data[:,4] = self.coordY   # [Mpc/h]
      data[:,5] = self.coordZ   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      data[:,6] = self.dX   # [Mpc/h]
      data[:,7] = self.dY   # [Mpc/h]
      data[:,8] = self.dZ   # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      data[:,9] = self.dXKaiser   # [Mpc/h] from cartesian catalog difference
      data[:,10] = self.dYKaiser   # [Mpc/h]
      data[:,11] = self.dZKaiser   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      data[:,12] = self.vX   #[km/s]
      data[:,13] = self.vY   #[km/s]
      data[:,14] = self.vZ   #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      data[:,15] = self.vR  # [km/s]   from spherical catalo
      data[:,16] = self.vTheta   # [km/s]
      data[:,17] = self.vPhi  # [km/s]
      #
      # Stellar mass
      data[:,18] = self.Mstellar   # [M_sun], from Maraston et al
      #
      # Halo mass
      data[:,19] = self.hasM  # flag=1 if mass is known
      data[:,20] = self.Mvir   # [M_sun]
      #
      np.savetxt(self.pathOutCatalog, data)


   def loadCatalog(self):
      print "- load full catalog from "+self.pathOutCatalog
      data = np.genfromtxt(self.pathOutCatalog)
      self.nObj = len(data[:,0])
      #
      # sky coordinates and redshift
      self.RA = data[:,0] # [deg]
      self.DEC = data[:,1]   # [deg]
      self.Z = data[:,2]
      #
      # observed cartesian coordinates
      self.coordX = data[:,3]   # [Mpc/h]
      self.coordY = data[:,4]   # [Mpc/h]
      self.coordZ = data[:,5]   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      self.dX = data[:,6]   # [Mpc/h]
      self.dY = data[:,7]   # [Mpc/h]
      self.dZ = data[:,8]   # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      self.dXKaiser = data[:,9]   # [Mpc/h] from cartesian catalog difference
      self.dYKaiser = data[:,10]   # [Mpc/h]
      self.dZKaiser = data[:,11]   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      self.vX = data[:,12]   #[km/s]
      self.vY = data[:,13]   #[km/s]
      self.vZ = data[:,14]   #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      self.vR = data[:,15]  # [km/s]   from spherical catalo
      self.vTheta = data[:,16]   # [km/s]
      self.vPhi = data[:,17]  # [km/s]
      #
      # Stellar masses
      self.Mstellar = data[:,18]   # [M_sun], from Maraston et al
      #
      # Halo mass
      self.hasM = data[:,19]
      self.Mvir = data[:,20]  # [M_sun]


   ##################################################################################
   ##################################################################################
   
   def addCatalog(self, newCat, save=False):
      """Combines the current catalog with a new catalog newCat.
      """
      # number of objects
      self.nObj += newCat.nObj
      #
      # sky coordinates and redshift
      self.RA = np.concatenate((self.RA, newCat.RA)) # [deg]
      self.DEC = np.concatenate((self.DEC, newCat.DEC))   # [deg]
      self.Z = np.concatenate((self.Z, newCat.Z))
      #
      # observed cartesian coordinates
      self.coordX = np.concatenate((self.coordX, newCat.coordX))   # [Mpc/h]
      self.coordY = np.concatenate((self.coordY, newCat.coordY))   # [Mpc/h]
      self.coordZ = np.concatenate((self.coordZ, newCat.coordZ))   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      self.dX = np.concatenate((self.dX, newCat.dX))  # [Mpc/h]
      self.dY = np.concatenate((self.dY, newCat.dY))   # [Mpc/h]
      self.dZ = np.concatenate((self.dZ, newCat.dZ))  # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      self.dXKaiser = np.concatenate((self.dXKaiser, newCat.dXKaiser))  # [Mpc/h] from cartesian catalog difference
      self.dYKaiser = np.concatenate((self.dYKaiser, newCat.dYKaiser))   # [Mpc/h]
      self.dZKaiser = np.concatenate((self.dZKaiser, newCat.dZKaiser))   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      self.vX = np.concatenate((self.vX, newCat.vX))   #[km/s]
      self.vY = np.concatenate((self.vY, newCat.vY))   #[km/s]
      self.vZ = np.concatenate((self.vZ, newCat.vZ))  #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      self.vR = np.concatenate((self.vR, newCat.vR))  # [km/s]   from spherical catalo
      self.vTheta = np.concatenate((self.vTheta, newCat.vTheta))   # [km/s]
      self.vPhi = np.concatenate((self.vPhi, newCat.vPhi))  # [km/s]
      #
      # Stellar masses
      self.Mstellar = np.concatenate((self.Mstellar, newCat.Mstellar))   # [M_sun], from Maraston et al
      #
      # Halo mass
      self.hasM = np.concatenate((self.hasM, newCat.hasM))
      self.Mvir = np.concatenate((self.Mvir, newCat.Mvir))  # [M_sun]

      # Write the full catalog to the output path, if needed
      if save:
         self.writeCatalog()



   ##################################################################################
   ##################################################################################
   
   def plotFootprint(self):
      """Overlay a scatter plot of the catalog positions on top of a healpix map,
      here the AdvACT hit count map.
      """
      fig=plt.figure(0)
      #
      # hit count map for AdvACT
      path = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/healpix_f150_daynight_all_div_mono.fits"
      hHitMap = hp.read_map(path)
      hp.mollview(np.log(np.abs(hHitMap)+1.e-5), fig=0, title="", coord=None, cbar=False, unit='')
      #
      # scatter plot of the catalog
      hp.projscatter(self.RA, self.DEC, alpha=0.01, lonlat=True, marker='.', c='r', rasterized=True)
      #
      fig.savefig(self.pathFig+"/footprint_"+self.name+".pdf", dpi=1200)
      fig.clf()


   ##################################################################################
   ##################################################################################
   
   def printProperties(self):
      print "Catalog: "+self.nameLong
      print "Number of objects = "+str(self.nObj)
      print "with mass: "+str(np.sum(self.hasM))+", i.e. fraction "+str(np.sum(self.hasM)/self.nObj)
      print "Z: mean = "+str(np.mean(self.Z))+", median = "+str(np.median(self.Z))
      m = self.Mstellar[self.hasM==1]
      print "M_star [M_sun]: mean = "+str(np.mean(m))+", median = "+str(np.median(m))
      m = self.Mvir[self.hasM==1]
      print "M_vir [M_sun]: mean = "+str(np.mean(m))+", median = "+str(np.median(m))
   


   ##################################################################################
   ##################################################################################

   def compareV1dRms(self):
      """expected std dev of velocities for LCDM:
      comparison between using median z or z-distribution
      """
      # interpolate RMS 1d velocity for speed
      f = lambda z: self.U.v3dRms(0., z, W3d_sth) / np.sqrt(3.)
      Z = np.linspace(0., 1., 201)
      V1dRms = np.array(map(f, Z))
      f = interp1d(Z, V1dRms, kind='linear', bounds_error=False, fill_value='extrapolate')
      
      print "Expected v1d_rms = "+str(np.mean(np.array(map(f, self.Z))))+" km/s"
      print "Expected v1d_rms(z=z_mean) = "+str(f(np.mean(self.Z)))+" km/s"
      print "RMS v_r, v_theta, v_phi = "+str(np.std(self.vR))+", "+str(np.std(self.vTheta))+", "+str(np.std(self.vPhi))+" km/s"



   ##################################################################################
   ##################################################################################
   
   
   def histogram(self, X, nBins=71, lim=(-1000., 1000.), sigma2Theory=None, name='x', nameLatex=r'$x$ [km/s]', semilogx=False, doGauss=False):
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
         histGaussFit *= self.nObj

      # Theory histogram
      if sigma2Theory is not None:
         av = 0.
         fPDF = lambda x: (2*np.pi*sigma2Theory)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2Theory))
         g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-3)[0]
         histTheory = np.array(map(g, range(nBins-1)))
         histTheory *= self.nObj

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
      ax.set_xlabel(nameLatex)
      ax.set_ylabel(r'number of objects')
      fig.savefig(self.pathFig+"/hist_"+name+".pdf")
      fig.clf()



   def plotHistograms(self):
      z0 = np.mean(self.Z)
      s2v1d = self.U.v3dRms(0., z0, W3d_sth)**2 / 3.
      
      # redshifts
      self.histogram(self.Z, nBins=71, lim=(0., 1.), name='z', nameLatex=r'$z$')
      
      # spherical velocities
      self.histogram(self.vR, nBins=71, lim=(-1000., 1000.), sigma2Theory=s2v1d, name='vr', nameLatex=r'$v_r$ [km/s]', doGauss=True)
      self.histogram(self.vTheta, nBins=71, lim=(-1000., 1000.), sigma2Theory=s2v1d, name='vtheta', nameLatex=r'$v_\theta$ [km/s]', doGauss=True)
      self.histogram(self.vPhi, nBins=71, lim=(-1000., 1000.), sigma2Theory=s2v1d, name='vphi', nameLatex=r'$v_\phi$ [km/s]', doGauss=True)
      
      # cartesian velocities
      self.histogram(self.vX, nBins=71, lim=(-1000., 1000.), sigma2Theory=s2v1d, name='vx', nameLatex=r'$v_x$ [km/s]', doGauss=True)
      self.histogram(self.vY, nBins=71, lim=(-1000., 1000.), sigma2Theory=s2v1d, name='vy', nameLatex=r'$v_y$ [km/s]', doGauss=True)
      self.histogram(self.vZ, nBins=71, lim=(-1000., 1000.), sigma2Theory=s2v1d, name='vz', nameLatex=r'$v_z$ [km/s]', doGauss=True)
      
      # stellar masses
      self.histogram(self.Mstellar, nBins=71, lim=(1.e9, 1.e12), name='mstellar', nameLatex=r'$M_\star$ [M$_\odot$]', semilogx=True)

      # virial masses
      self.histogram(self.Mvir, nBins=71, lim=(1.e12, 1.e15), name='mvir', nameLatex=r'$M_\text{vir}$ [M$_\odot$]', semilogx=True)
      
      # comoving virial radius
      # need masses in Msun/h
      Par = zip(self.Mvir*self.U.bg.h, self.Z)
      f = lambda par: self.U.frvir(par[0], par[1])   # in: Msun/h, out: Mpc/h
      Rvir = np.array(map(f, Par))  # in Mpc/h
      #Rvir /= self.U.bg.h  # Mpc
      self.histogram(Rvir/self.U.bg.h, nBins=71, lim=(0., 2.5), name='rvir', nameLatex=r'$R_\text{vir}$ [Mpc]')
      
      # virial angular radius
      Chi = np.array(map(self.U.bg.comoving_distance, self.Z)) # Mpc/h
      Thetavir = Rvir / Chi   # rad
      self.histogram(Thetavir * (180.*60./np.pi), nBins=71, lim=(0.5, 3.), name='thetavir', nameLatex=r'$\theta_\text{vir}$ [arcmin]')
      
      # expected kSZ?
      
      # expected tSZ?

      # displacements?


##################################################################################
##################################################################################
#
#class CMASS_S_Mariana(Catalog):
#
#   def __init__(self, U, MassConversion, save=False):
#      # galaxy or cluster catalog
#      self.name = "cmass_s_mariana"
#      self.nameLong = "CMASS S M"
#
#      # path to the input catalog
#      self.pathInCatalog = "../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt"
#
#      super(CMASS_S_Mariana, self).__init__(U, MassConversion, save=save)
#
#
###################################################################################
#
#class CMASS_N_Mariana(Catalog):
#
#   def __init__(self, U, MassConversion, save=False):
#      # galaxy or cluster catalog
#      self.name = "cmass_n_mariana"
#      self.nameLong = "CMASS N M"
#
#      # path to the input catalog
#      self.pathInCatalog = "../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt"
#
#      super(CMASS_N_Mariana, self).__init__(U, MassConversion, save=save)
#
#
###################################################################################
###################################################################################
#
#class CMASS_S_Kendrick(Catalog):
#
#   def __init__(self, U, MassConversion, save=False):
#      # galaxy or cluster catalog
#      self.name = "cmass_s_kendrick"
#      self.nameLong = "CMASS S K"
#
#      # path to the input catalog
#      self.pathInCatalog = "../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_S_kendrick.txt"
#
#      super(CMASS_S_Kendrick, self).__init__(U, MassConversion, save=save)
#
#
###################################################################################
#
#class CMASS_N_Kendrick(Catalog):
#
#   def __init__(self, U, MassConversion, save=False):
#      # galaxy or cluster catalog
#      self.name = "cmass_n_kendrick"
#      self.nameLong = "CMASS N K"
#
#      # path to the input catalog
#      self.pathInCatalog = "../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_N_kendrick.txt"
#
#      super(CMASS_N_Kendrick, self).__init__(U, MassConversion, save=save)
#
#
###################################################################################
#
#class LOWZ_S_Kendrick(Catalog):
#
#   def __init__(self, U, MassConversion, save=False):
#      # galaxy or cluster catalog
#      self.name = "lowz_s_kendrick"
#      self.nameLong = "LOWZ S K"
#
#      # path to the input catalog
#      self.pathInCatalog = "../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_S_kendrick.txt"
#
#      super(LOWZ_S_Kendrick, self).__init__(U, MassConversion, save=save)
#
#
###################################################################################
#
#class LOWZ_N_Kendrick(Catalog):
#
#   def __init__(self, U, MassConversion, save=False):
#      # galaxy or cluster catalog
#      self.name = "lowz_n_kendrick"
#      self.nameLong = "LOWZ N K"
#
#      # path to the input catalog
#      self.pathInCatalog = "../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_N_kendrick.txt"
#
#      super(LOWZ_N_Kendrick, self).__init__(U, MassConversion, save=save)
#
#
###################################################################################
###################################################################################
#
#
