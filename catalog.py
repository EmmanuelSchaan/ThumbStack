from headers import *

##################################################################################
##################################################################################

class Catalog(object):

   def __init__(self, U, MassConversion, save=False):

      self.U = U
      self.MassConversion = MassConversion

      # galaxy or cluster catalog name
      #self.name = "cmasssouth"
      #self.nameLong = "CMASS South"
      
      # path to the input catalog
      #self.pathInCatalog = "./input/?.txt"
      
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


   ##################################################################################

   def writeCatalog(self):
      print "- write full catalog to "+self.pathOutCatalog
      data = np.zeros((self.nObj,20))
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
      data[:,19] = self.Mvir   # [M_sun]
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
      self.Mvir = data[:,19]  # [M_sun]



   ##################################################################################
   ##################################################################################
   
   def ftestVel(self):
      """expected std dev of velocities for LCDM
      comparison between using median z or z-distribution
      """

      z0 = np.mean(self.Z)
      
      f = lambda z: self.U.v3dRms(0., z, W3d_sth)**2 / 3.
      Sigma2 = np.array(map(f, self.Z))
      sigma2 = np.mean(Sigma2)
   
      print sigma2
      print self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.


   ##################################################################################

   def plotVelSphHistogram(self):

      z0 = 0.57
      
      #################################################################
      # velR
      
      velR = np.copy(self.velR)
      mean = np.mean(velR)
      std = np.std(velR)
      print std

      # histogram parameters
      nbins = 71
      firstbin = -1000.
      lastbin = 1000.
      binwidth = (lastbin-firstbin)/(nbins-1)
      Bins = np.linspace(firstbin, lastbin, nbins)

      # histograms from data
      histVelR = np.histogram(velR, Bins)[0]
      histVelR = histVelR.astype(float)

      # histogram for a Gaussian with the variance from the data
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histGaussFit = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histGaussFit *= self.nObj

      # histogram for LCDM
      #f = lambda z: self.U.RMSVelocity(0., z, W3d_sth)**2 / 3.
      #Sigma2 = np.array(map(f, self.Z))
      #sigma2 = np.mean(Sigma2)
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.nObj

      # velR
      fig = plt.figure(0)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histVelR, binwidth, color='b', alpha=0.5, label=r'rec. $v_{r}$')
      ax.step(Bins[:-1], histLCDM, color='r', lw=3, where='post', label=r'$\Lambda CDM$')
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      #
      infos = str(len(velR))+" objects"+"\n"+"mean="+str(round(mean,1))+"km/s"+"\n"+r"$\sigma$="+str(round(std,1))+"km/s"+"\n"+r"expect $\sigma$="+str(round(np.sqrt(sigma2),1))+r"km/s for $\Lambda CDM$."
      ax.text(-900.,  10500., infos, fontsize=16,
              horizontalalignment='left',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      #ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$v_r$ [km/s]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/velocities/hist_vr.pdf")
      
      #################################################################
      # velTheta
      
      velTheta = np.copy(self.velTheta)
      mean = np.mean(velTheta)
      std = np.std(velTheta)
      print std

      # histogram parameters
      nbins = 71
      firstbin = -1000.
      lastbin = 1000.
      binwidth = (lastbin-firstbin)/(nbins-1)
      Bins = np.linspace(firstbin, lastbin, nbins)

      # histograms from data
      histVelTheta = np.histogram(velTheta, Bins)[0]
      histVelTheta = histVelTheta.astype(float)

      # histogram for a Gaussian with the variance from the data
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histGaussFit = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histGaussFit *= self.nObj

      # histogram for LCDM
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.nObj

      # velR
      fig = plt.figure(1)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histVelTheta, binwidth, color='b', alpha=0.5, label=r'rec. $v_{\theta}$')
      ax.step(Bins[:-1], histLCDM, color='r', lw=3, where='post', label=r'$\Lambda CDM$')
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      #
      infos = str(len(velTheta))+" objects"+"\n"+"mean="+str(round(mean,1))+"km/s"+"\n"+r"$\sigma$="+str(round(std,1))+"km/s"+"\n"+r"expect $\sigma$="+str(round(np.sqrt(sigma2),1))+r"km/s for $\Lambda CDM$."
      ax.text(-900.,  10500., infos, fontsize=16,
              horizontalalignment='left',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$v_\theta$ [km/s]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/velocities/hist_vtheta.pdf")
      
      #################################################################
      # velPhi
      
      velPhi = np.copy(self.velPhi)
      mean = np.mean(velPhi)
      std = np.std(velPhi)
      print std
      
      # histogram parameters
      nbins = 71
      firstbin = -1000.
      lastbin = 1000.
      binwidth = (lastbin-firstbin)/(nbins-1)
      Bins = np.linspace(firstbin, lastbin, nbins)
      
      # histograms from data
      histVelPhi = np.histogram(velPhi, Bins)[0]
      histVelPhi = histVelPhi.astype(float)
      
      # histogram for a Gaussian with the variance from the data
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histGaussFit = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histGaussFit *= self.nObj
      
      # histogram for LCDM
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.nObj
      
      # velR
      fig = plt.figure(2)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histVelPhi, binwidth, color='b', alpha=0.5, label=r'rec. $v_{\phi}$')
      ax.step(Bins[:-1], histLCDM, color='r', lw=3, where='post', label=r'$\Lambda CDM$')
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      #
      infos = str(len(velPhi))+" objects"+"\n"+"mean="+str(round(mean,1))+"km/s"+"\n"+r"$\sigma$="+str(round(std,1))+"km/s"+"\n"+r"expect $\sigma$="+str(round(np.sqrt(sigma2),1))+r"km/s for $\Lambda CDM$."
      ax.text(-900.,  10000., infos, fontsize=16,
              horizontalalignment='left',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$v_\phi$ [km/s]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/velocities/hist_vphi.pdf")
      
      plt.show()


   ##################################################################################

   def plotVelCartHistogram(self):
      
      z0 = 0.57
      
      #################################################################
      # velX
      
      velX = np.copy(self.velX)
      mean = np.mean(velX)
      std = np.std(velX)
      
      # histogram parameters
      nbins = 71
      firstbin = -1000.
      lastbin = 1000.
      binwidth = (lastbin-firstbin)/(nbins-1)
      Bins = np.linspace(firstbin, lastbin, nbins)
      
      # histograms from data
      histVel = np.histogram(velX, Bins)[0]
      histVel = histVel.astype(float)
      
      # histogram for a Gaussian with the variance from the data
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histGaussFit = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histGaussFit *= self.nObj
      
      # histogram for LCDM
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.nObj
      
      # velX
      fig = plt.figure(0)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histVel, binwidth, color='b', alpha=0.5, label=r'rec. $v_{x}$')
      ax.step(Bins[:-1], histLCDM, color='r', lw=3, where='post', label=r'$\Lambda CDM$')
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      #
      infos = str(len(velX))+" objects"+"\n"+"mean="+str(round(mean,1))+"km/s"+"\n"+r"$\sigma$="+str(round(std,1))+"km/s"+"\n"+r"expect $\sigma$="+str(round(np.sqrt(sigma2),1))+r"km/s for $\Lambda CDM$."
      ax.text(-900.,  10500., infos, fontsize=16,
              horizontalalignment='left',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$v_x$ [km/s]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/velocities/hist_vx.pdf")

      #################################################################
      # velY

      velY = np.copy(self.velY)
      mean = np.mean(velY)
      std = np.std(velY)

      # histogram parameters
      nbins = 71
      firstbin = -1000.
      lastbin = 1000.
      binwidth = (lastbin-firstbin)/(nbins-1)
      Bins = np.linspace(firstbin, lastbin, nbins)

      # histograms from data
      histVel = np.histogram(velY, Bins)[0]
      histVel = histVel.astype(float)

      # histogram for a Gaussian with the variance from the data
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histGaussFit = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histGaussFit *= self.nObj

      # histogram for LCDM
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.nObj

      # velY
      fig = plt.figure(1)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histVel, binwidth, color='b', alpha=0.5, label=r'rec. $v_y$')
      ax.step(Bins[:-1], histLCDM, color='r', lw=3, where='post', label=r'$\Lambda CDM$')
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      #
      infos = str(len(velY))+" objects"+"\n"+"mean="+str(round(mean,1))+"km/s"+"\n"+r"$\sigma$="+str(round(std,1))+"km/s"+"\n"+r"expect $\sigma$="+str(round(np.sqrt(sigma2),1))+r"km/s for $\Lambda CDM$."
      ax.text(-900.,  10000., infos, fontsize=16,
             horizontalalignment='left',
             verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$v_y$ [km/s]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/velocities/hist_vy.pdf")

      #################################################################
      # velZ

      velZ = np.copy(self.velZ)
      mean = np.mean(velZ)
      std = np.std(velZ)

      # histogram parameters
      nbins = 71
      firstbin = -1000.
      lastbin = 1000.
      binwidth = (lastbin-firstbin)/(nbins-1)
      Bins = np.linspace(firstbin, lastbin, nbins)

      # histograms from data
      histVel = np.histogram(velZ, Bins)[0]
      histVel = histVel.astype(float)

      # histogram for a Gaussian with the variance from the data
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histGaussFit = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histGaussFit *= self.nObj

      # histogram for LCDM
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.nObj

      # velZ
      fig = plt.figure(2)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histVel, binwidth, color='b', alpha=0.5, label=r'rec. $v_z$')
      ax.step(Bins[:-1], histLCDM, color='r', lw=3, where='post', label=r'$\Lambda CDM$')
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      #
      infos = str(len(velZ))+" objects"+"\n"+"mean="+str(round(mean,1))+"km/s"+"\n"+r"$\sigma$="+str(round(std,1))+"km/s"+"\n"+r"expect $\sigma$="+str(round(np.sqrt(sigma2),1))+r"km/s for $\Lambda CDM$."
      ax.text(-900.,  10000., infos, fontsize=16,
            horizontalalignment='left',
            verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$v_z$ [km/s]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/velocities/hist_vz.pdf")

      plt.show()


   ##################################################################################

   def plotMassHistogram(self):

      ################################################################################
      # stellar masses
      
      Mstellar = np.delete(self.Mstellar, np.where(self.Mstellar==0)[0])
      mean = np.mean(Mstellar)
      std = np.std(Mstellar)

      # histogram parameters
      nbins = 71
      firstbin = np.min(Mstellar)
      lastbin = np.max(Mstellar)
      Bins = np.logspace(np.log10(firstbin), np.log10(lastbin), nbins, 10.)
      Binwidth = Bins[1:]-Bins[:-1]
      
      # histograms from data
      histMstellar = np.histogram(Mstellar, Bins)[0]
      histMstellar = histMstellar.astype(float)

      # Mstellar
      fig = plt.figure(0)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histMstellar, Binwidth, color='b', alpha=0.7, label=r'CMASS South galaxies')
      #
      infos = str(len(Mstellar))+" objects"+"\n"+"mean="+format(mean, '.2e')+r"$M_\odot$"+"\n"+r"$\sigma$="+format(std, '.2e')+r"$M_\odot$"
      ax.text(firstbin,  10500., infos, fontsize=16,
              horizontalalignment='left',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xscale('log')
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$M_*$ [$M_\odot$]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/mass_conversion/hist_mStellar.pdf")
      
      
      ################################################################################
      # virial masses
      
      Mvir = np.delete(self.Mvir, np.where((self.Mvir==0)+np.isnan(self.Mvir))[0])
      mean = np.mean(Mvir)
      std = np.std(Mvir)

      # histogram parameters
      nbins = 71
      firstbin = np.min(Mvir)
      lastbin = np.max(Mvir)
      Bins = np.logspace(np.log10(firstbin), np.log10(lastbin), nbins, 10.)
      Binwidth = Bins[1:]-Bins[:-1]
      
      # histograms from data
      histMvir = np.histogram(Mvir, Bins)[0]
      histMvir = histMvir.astype(float)
      
      # Mvir
      fig = plt.figure(1)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histMvir, Binwidth, color='b', alpha=0.7, label=r'CMASS South galaxies')
      #
      infos = str(len(Mvir))+" objects"+"\n"+"mean="+format(mean, '.2e')+r"$M_\odot$"+"\n"+r"$\sigma$="+format(std, '.2e')+r"$M_\odot$"
      ax.text(firstbin,  5000., infos, fontsize=16,
              horizontalalignment='left',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xscale('log')
      #ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$M_\text{vir}$ [$M_\odot$]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/mass_conversion/hist_mVir_Kravtsov14.pdf")
      #fig.savefig("./figures/"+self.name+"/mass_conversion/hist_mVir_PC11.pdf")
      
      ################################################################################
#      # Mvir * dn/dMvir
#      fig = plt.figure(11)
#      ax = fig.add_subplot(111)
#      #
#      ax.bar(Bins[:-1], Bins[:-1]*histMvir, Binwidth, color='b', alpha=0.7, label=r'CMASS South galaxies')
#      #
#      infos = str(len(Mvir))+" objects"+"\n"+"mean="+format(mean, '.2e')+r"$M_\odot$"+"\n"+r"$\sigma$="+format(std, '.2e')+r"$M_\odot$"
#      ax.text(firstbin,  5000., infos, fontsize=16,
#              horizontalalignment='left',
#              verticalalignment='center')
#      #
#      ax.legend(loc=1)
#      ax.set_xscale('log')
#      #ax.set_yscale('log', nonposy='clip')
#      ax.set_xlim((firstbin, lastbin))
#      ax.set_xlabel(r'$M_\text{vir}$ [$M_\odot$]')
#      ax.set_ylabel(r'$M_\text{vir} * dn/dM_\text{vir}$')


      ################################################################################
      # virial radius

      f = lambda z: ( 3./(4*np.pi*self.U.rho_z(z)*self.U.DeltaVir(z)) )**(1./3.)
      Z = np.delete(self.Z, np.where((self.Mvir==0)+np.isnan(self.Mvir))[0])
      rVir = Mvir**(1./3.) * np.array(map(f, Z))   # rVir
      mean = np.mean(rVir)
      std = np.std(rVir)

      # histogram parameters
      nbins = 71
      firstbin = np.min(rVir)
      lastbin = np.max(rVir)
      #      Bins = np.logspace(np.log10(firstbin), np.log10(lastbin), nbins, 10.)
      Bins = np.linspace(firstbin, lastbin, nbins)
      Binwidth = Bins[1:]-Bins[:-1]

      # histograms from data
      histRVir = np.histogram(rVir, Bins)[0]
      histRVir = histRVir.astype(float)

      # rVir
      fig = plt.figure(2)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histRVir, Binwidth, color='b', alpha=0.7, label=r'CMASS South galaxies')
      #
      infos = str(len(Mvir))+" objects"+"\n"+"mean="+format(mean, '.1f')+r"Mpc/h"+"\n"+r"$\sigma$="+format(std, '.1f')+r"Mpc/h"
      ax.text(lastbin,  5000., infos, fontsize=16,
              horizontalalignment='right',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      #ax.set_xscale('log')
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$R_\text{vir}$ [Mpc/$h$]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/mass_conversion/hist_rVir_Kravtsov14.pdf")
      #fig.savefig("./figures/"+self.name+"/mass_conversion/hist_rVir_PC11.pdf")

 
      ################################################################################
      # angular virial radius
      
      f = lambda z: ( 3./(4*np.pi*self.U.rho_z(z)*self.U.DeltaVir(z)) )**(1./3.) / self.U.ComovDist(1./(1.+z), 1.)
      Z = np.delete(self.Z, np.where((self.Mvir==0)+np.isnan(self.Mvir))[0])
      thetaVir = Mvir**(1./3.) * np.array(map(f, Z))   # thetaVir in radians
      thetaVir *= (180./np.pi)*60.  # thetaVir in arcmin
      mean = np.mean(thetaVir)
      std = np.std(thetaVir)
      
      # histogram parameters
      nbins = 71
      firstbin = np.min(thetaVir)
      lastbin = np.max(thetaVir)
#      Bins = np.logspace(np.log10(firstbin), np.log10(lastbin), nbins, 10.)
      Bins = np.linspace(firstbin, lastbin, nbins)
      Binwidth = Bins[1:]-Bins[:-1]
      
      # histograms from data
      histThetaVir = np.histogram(thetaVir, Bins)[0]
      histThetaVir = histThetaVir.astype(float)
      
      # thetaVir
      fig = plt.figure(3)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histThetaVir, Binwidth, color='b', alpha=0.7, label=r'CMASS South galaxies')
      #
      infos = str(len(Mvir))+" objects"+"\n"+"mean="+format(mean, '.1f')+r"arcmin"+"\n"+r"$\sigma$="+format(std, '.1f')+r"arcmin"
      ax.text(lastbin,  5000., infos, fontsize=16,
              horizontalalignment='right',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      #ax.set_xscale('log')
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$\theta_\text{vir}$ [arcmin]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/mass_conversion/hist_thetaVir_Kravtsov14.pdf")
      #fig.savefig("./figures/"+self.name+"/mass_conversion/hist_thetaVir_PC11.pdf")

      plt.show()



   def plotKSZHistogram(self):
      
      z0 = 0.57
      
      #################################################################
      # velR * Mvir
      
      X = self.Mvir * self.velR
      mean = np.mean(X)
      std = np.std(X)
      print std
      
      # histogram parameters
      nbins = 71
      firstbin = min(X)
      lastbin = max(X)
      binwidth = (lastbin-firstbin)/(nbins-1)
      Bins = np.linspace(firstbin, lastbin, nbins)
      
      # histograms from data
      hist = np.histogram(X, Bins)[0]
      hist = hist.astype(float)
      
      # histogram for a Gaussian with the variance from the data
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histGaussFit = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histGaussFit *= self.nObj
      
      # plot
      fig = plt.figure(0)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], hist, binwidth, color='b', alpha=0.5, label=r'$M_\text{vir} v_{r}$')
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      #
      infos = str(len(X))+" objects"+"\n"+"mean="+str(round(mean,1))+"\n"+r"$\sigma$="+str(round(std,1))
      ax.text(-900.,  10500., infos, fontsize=16,
              horizontalalignment='left',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_yscale('log', nonposy='clip')
      ax.set_xlim((firstbin, lastbin))
      ax.set_ylim((1.e-1, 1.e5))
      ax.set_xlabel(r'$M_\text{vir} v_r$')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/actdeep5656/hist_mvir_vr.pdf")

      plt.show()


   ##################################################################################

   def plotDispSphHistogram(self):
      
      z0 = 0.57
      
      #################################################################
      # dispR
      
      disp = np.copy(self.dispR)
      # undo the correction for redshift-space correction
      disp /= 1. + self.U.growthLogDerivativeF(z0)
      mean = np.mean(disp)
      std = np.std(disp)
      
      # histogram parameters
      nbins = 71
      firstbin = -10.
      lastbin = 10.
      binwidth = (lastbin-firstbin)/(nbins-1)
      Bins = np.linspace(firstbin, lastbin, nbins)
      
      # histograms from data
      hist = np.histogram(disp, Bins)[0]
      hist = hist.astype(float)
      
      # histogram for a Gaussian with the variance from the data
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histGaussFit = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histGaussFit *= self.nObj
      
      # histogram for LCDM
      sigma2 = self.U.RMSDisplacement(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.nObj
      
      # dispR
      fig = plt.figure(0)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], hist, binwidth, color='b', alpha=0.5, label=r'rec. $\psi_{r}$')
      ax.step(Bins[:-1], histLCDM, color='r', lw=3, where='post', label=r'$\Lambda CDM$')
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      #
      infos = str(len(disp))+" objects"+"\n"+"mean="+str(round(mean,1))+"Mpc/h"+"\n"+r"$\sigma$="+str(round(std,1))+"Mpc/h"+"\n"+r"expect $\sigma$="+str(round(np.sqrt(sigma2),1))+r"Mpc/h for $\Lambda CDM$."
      ax.text(-10.,  8000., infos, fontsize=16,
              horizontalalignment='left',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$\psi_r$ [Mpc/h]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/velocities/hist_dispr.pdf")


      #################################################################
      # dispTheta

      disp = np.copy(self.dispTheta)
      mean = np.mean(disp)
      std = np.std(disp)

      # histogram parameters
      nbins = 71
      firstbin = -10.
      lastbin = 10.
      binwidth = (lastbin-firstbin)/(nbins-1)
      Bins = np.linspace(firstbin, lastbin, nbins)

      # histograms from data
      hist = np.histogram(disp, Bins)[0]
      hist = hist.astype(float)

      # histogram for a Gaussian with the variance from the data
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histGaussFit = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histGaussFit *= self.nObj

      # histogram for LCDM
      sigma2 = self.U.RMSDisplacement(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.nObj

      # dispTheta
      fig = plt.figure(1)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], hist, binwidth, color='b', alpha=0.5, label=r'rec. $\psi_{\theta}$')
      ax.step(Bins[:-1], histLCDM, color='r', lw=3, where='post', label=r'$\Lambda CDM$')
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      #
      infos = str(len(disp))+" objects"+"\n"+"mean="+str(round(mean,1))+"Mpc/h"+"\n"+r"$\sigma$="+str(round(std,1))+"Mpc/h"+"\n"+r"expect $\sigma$="+str(round(np.sqrt(sigma2),1))+r"Mpc/h for $\Lambda CDM$."
      ax.text(-10.,  8000., infos, fontsize=16,
          horizontalalignment='left',
          verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$\psi_\theta$ [Mpc/h]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/velocities/hist_disptheta.pdf")



      #################################################################
      # dispPhi

      disp = np.copy(self.dispPhi)
      mean = np.mean(disp)
      std = np.std(disp)

      # histogram parameters
      nbins = 71
      firstbin = -10.
      lastbin = 10.
      binwidth = (lastbin-firstbin)/(nbins-1)
      Bins = np.linspace(firstbin, lastbin, nbins)

      # histograms from data
      hist = np.histogram(disp, Bins)[0]
      hist = hist.astype(float)

      # histogram for a Gaussian with the variance from the data
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histGaussFit = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histGaussFit *= self.nObj

      # histogram for LCDM
      sigma2 = self.U.RMSDisplacement(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.nObj

      # dispTheta
      fig = plt.figure(2)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], hist, binwidth, color='b', alpha=0.5, label=r'rec. $\psi_{\phi}$')
      ax.step(Bins[:-1], histLCDM, color='r', lw=3, where='post', label=r'$\Lambda CDM$')
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      #
      infos = str(len(disp))+" objects"+"\n"+"mean="+str(round(mean,1))+"Mpc/h"+"\n"+r"$\sigma$="+str(round(std,1))+"Mpc/h"+"\n"+r"expect $\sigma$="+str(round(np.sqrt(sigma2),1))+r"Mpc/h for $\Lambda CDM$."
      ax.text(-10.,  8000., infos, fontsize=16,
            horizontalalignment='left',
            verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$\psi_\phi$ [Mpc/h]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.name+"/velocities/hist_dispphi.pdf")
      
      
      plt.show()



##################################################################################
##################################################################################

class CMASS_S_Mariana(Catalog):
   
   def __init__(self, U, MassConversion, save=False):
      # galaxy or cluster catalog
      self.name = "cmass_s_mariana"
      self.nameLong = "CMASS S M"
      
      # path to the input catalog
      self.pathInCatalog = "../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_S_mariana.txt"
   
      super(CMASS_S_Mariana, self).__init__(U, MassConversion, save=save)


##################################################################################

class CMASS_N_Mariana(Catalog):
   
   def __init__(self, U, MassConversion, save=False):
      # galaxy or cluster catalog
      self.name = "cmass_n_mariana"
      self.nameLong = "CMASS N M"
      
      # path to the input catalog
      self.pathInCatalog = "../../data/CMASS_DR12_mariana_20160200/output/cmass_dr12_N_mariana.txt"
   
      super(CMASS_N_Mariana, self).__init__(U, MassConversion, save=save)


##################################################################################
##################################################################################

class CMASS_S_Kendrick(Catalog):
   
   def __init__(self, U, MassConversion, save=False):
      # galaxy or cluster catalog
      self.name = "cmass_s_kendrick"
      self.nameLong = "CMASS S K"
      
      # path to the input catalog
      self.pathInCatalog = "../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_S_kendrick.txt"
   
      super(CMASS_S_Kendrick, self).__init__(U, MassConversion, save=save)


##################################################################################

class CMASS_N_Kendrick(Catalog):
   
   def __init__(self, U, MassConversion, save=False):
      # galaxy or cluster catalog
      self.name = "cmass_n_kendrick"
      self.nameLong = "CMASS N K"
      
      # path to the input catalog
      self.pathInCatalog = "../../data/BOSS_DR10_kendrick_20150407/output/cmass_dr10_N_kendrick.txt"
   
      super(CMASS_N_Kendrick, self).__init__(U, MassConversion, save=save)


##################################################################################

class LOWZ_S_Kendrick(Catalog):
   
   def __init__(self, U, MassConversion, save=False):
      # galaxy or cluster catalog
      self.name = "lowz_s_kendrick"
      self.nameLong = "LOWZ S K"
      
      # path to the input catalog
      self.pathInCatalog = "../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_S_kendrick.txt"
   
      super(LOWZ_S_Kendrick, self).__init__(U, MassConversion, save=save)


##################################################################################

class LOWZ_N_Kendrick(Catalog):
   
   def __init__(self, U, MassConversion, save=False):
      # galaxy or cluster catalog
      self.name = "lowz_n_kendrick"
      self.nameLong = "LOWZ N K"
      
      # path to the input catalog
      self.pathInCatalog = "../../data/BOSS_DR10_kendrick_20150407/output/lowz_dr10_N_kendrick.txt"
   
      super(LOWZ_N_Kendrick, self).__init__(U, MassConversion, save=save)


##################################################################################
##################################################################################


