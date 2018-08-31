from headers import *

##################################################################################


class CheckCMASSSouth(object):
   
   def __init__(self, U, save=False):
      
      # galaxy or cluster catalog
      self.nameCatalog = "cmasssouth"
      self.nameCatalogLong = "CMASS South"
      self.pathCatalog = "./input/cmasssouth_dr12.txt"
      
      # path to the catalog with halo masses
      self.pathOutputCatalog = "./output/"+self.nameCatalog+"/objects.txt"
      
      self.U = U
      
      # load the conversion from stellar to halo mass
      self.floadMassConversion()
      self.floadCatalog(save=save)
   


   ##################################################################################
   
   def faddHaloMass(self):
      data = np.genfromtxt(self.pathCatalog)
      Mstellar = data[:,12]
      Nobj = len(Mstellar)
      
      # keep all objects with a known stellar mass
      Flag = np.ones(Nobj)
      Mvir = np.zeros(Nobj)
      for iobj in range(Nobj):
         mStellar = Mstellar[iobj]
         if (mStellar<1.) or np.isnan(mStellar):
            Flag[iobj] = 0.
         else:
            Mvir[iobj] = self.fMstellarToMvir(mStellar)
               
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
 
      newData = np.zeros((nObjRemaining, 14))
      newData[:,:13] = data[iRemaining,:13]
      newData[:,13] = Mvir[iRemaining]
      print "- save catalog with halo masses to "+self.pathOutputCatalog
      np.savetxt(self.pathOutputCatalog, newData)
   
   
   def floadCatalog(self, save=False):
      
      # add halo mass to catalog if not already done
      if (not os.path.isfile(self.pathOutputCatalog)) or save:
         self.faddHaloMass()

      print "- read catalog with halo masses from "+self.pathOutputCatalog
      data = np.genfromtxt(self.pathOutputCatalog)
      self.coordX = data[:,0]
      self.coordY = data[:,1]
      self.coordZ = data[:,2]
      self.RA = data[:,3]
      self.DEC = data[:,4]
      self.Z = data[:,5]
      self.velX = data[:,6]
      self.velY = data[:,7]
      self.velZ = data[:,8]
      self.velR = data[:,9]
      self.velTheta = data[:,10]
      self.velPhi = data[:,11]
      self.Mstellar = data[:,12]
      self.Mvir = data[:,13]
      self.Nobj = len(self.coordX)
      
      
      # read back-displacements from reconstruction, in spherical coord.
      # just to plot their histograms and compare it to the ones of velocities
      dispData = np.genfromtxt("../../data/CMASS_DR12/displacements/dr_dtheta_dphi_DR12_CMASS_S.txt")
      # get the forward displacements
      # watch out that dispR has the redshift-space correction in it
      self.dispR = -dispData[:,0]
      self.dispTheta = -dispData[:,1]
      self.dispPhi = -dispData[:,2]


   ###############################################################################
   
   def floadMassConversion(self):
      
      # from Planck Intermediate 11, fig3
      # should be M200 and not Mvir...
      data = np.genfromtxt("./input/digitizing_PCI11_fig3/PCI11_fig3.txt")
      flogMstellarToLogMvir = UnivariateSpline(data[:,0], data[:,1], k=3,s=0)
      m_min = 10.**np.min(data[:,0])
      m_max = 10.**np.max(data[:,0])
      self.fMstellarToMvir = lambda mStellar: 10.**flogMstellarToLogMvir(np.log10(mStellar)) * (mStellar>m_min and mStellar<m_max)
      '''
      # from Kravtsov et al 2014, fig17
      # should be M200 and not Mvir
      # watch for the different order of M200 and Mstellar
      data = np.genfromtxt("./input/digitizing_kravtsov2014_fig17/kravtsov2014_fig17.txt")
      ffMstellarToMvir = UnivariateSpline(data[:,1], data[:,0], k=3,s=0)
      m_min = np.min(data[:,1])
      m_max = np.max(data[:,1])
      self.fMstellarToMvir = lambda mStellar: ffMstellarToMvir(mStellar) * (mStellar>m_min and mStellar<m_max)
      '''

   ##################################################################################
   
   # expected std dev of velocities for LCDM
   # comparison between using median z or z-distribution
   def ftestVel(self):
      z0 = 0.57
      
      f = lambda z: self.U.RMSVelocity(0., z, W3d_sth)**2 / 3.
      Sigma2 = np.array(map(f, self.Z))
      sigma2 = np.mean(Sigma2)
   
      print sigma2
      print self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.

   
   
   
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
      histGaussFit *= self.Nobj

      # histogram for LCDM
      #f = lambda z: self.U.RMSVelocity(0., z, W3d_sth)**2 / 3.
      #Sigma2 = np.array(map(f, self.Z))
      #sigma2 = np.mean(Sigma2)
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.Nobj

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
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$v_r$ [km/s]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.nameCatalog+"/velocities/hist_vr.pdf")
      
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
      histGaussFit *= self.Nobj

      # histogram for LCDM
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.Nobj

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
      #fig.savefig("./figures/"+self.nameCatalog+"/velocities/hist_vtheta.pdf")
      
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
      histGaussFit *= self.Nobj
      
      # histogram for LCDM
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.Nobj
      
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
      #fig.savefig("./figures/"+self.nameCatalog+"/velocities/hist_vphi.pdf")
      
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
      histGaussFit *= self.Nobj
      
      # histogram for LCDM
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.Nobj
      
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
      #fig.savefig("./figures/"+self.nameCatalog+"/velocities/hist_vx.pdf")

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
      histGaussFit *= self.Nobj

      # histogram for LCDM
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.Nobj

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
      #fig.savefig("./figures/"+self.nameCatalog+"/velocities/hist_vy.pdf")

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
      histGaussFit *= self.Nobj

      # histogram for LCDM
      sigma2 = self.U.RMSVelocity(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.Nobj

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
      #fig.savefig("./figures/"+self.nameCatalog+"/velocities/hist_vz.pdf")

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
      #fig.savefig("./figures/"+self.nameCatalog+"/mass_conversion/hist_mStellar.pdf")
      
      
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
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$M_\text{vir}$ [$M_\odot$]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.nameCatalog+"/mass_conversion/hist_mVir_Kravtsov14.pdf")
      #fig.savefig("./figures/"+self.nameCatalog+"/mass_conversion/hist_mVir_PC11.pdf")
 
 
      ################################################################################
      # virial radii
      
      
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
      Bins = np.logspace(np.log10(firstbin), np.log10(lastbin), nbins, 10.)
      Binwidth = Bins[1:]-Bins[:-1]
      
      # histograms from data
      histThetaVir = np.histogram(thetaVir, Bins)[0]
      histThetaVir = histThetaVir.astype(float)
      
      # thetaVir
      fig = plt.figure(2)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histThetaVir, Binwidth, color='b', alpha=0.7, label=r'CMASS South galaxies')
      #
      infos = str(len(Mvir))+" objects"+"\n"+"mean="+format(mean, '.1f')+r"arcmin"+"\n"+r"$\sigma$="+format(std, '.1f')+r"arcmin"
      ax.text(firstbin,  5000., infos, fontsize=16,
              horizontalalignment='left',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xscale('log')
      ax.set_xlim((firstbin, lastbin))
      ax.set_xlabel(r'$\theta_\text{vir}$ [arcmin]')
      ax.set_ylabel(r'number of objects')
      #fig.savefig("./figures/"+self.nameCatalog+"/mass_conversion/hist_thetaVir_Kravtsov14.pdf")
      #fig.savefig("./figures/"+self.nameCatalog+"/mass_conversion/hist_thetaVir_PC11.pdf")

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
      histGaussFit *= self.Nobj
      
      # histogram for LCDM
      sigma2 = self.U.RMSDisplacement(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.Nobj
      
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
      #fig.savefig("./figures/"+self.nameCatalog+"/velocities/hist_dispr.pdf")


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
      histGaussFit *= self.Nobj

      # histogram for LCDM
      sigma2 = self.U.RMSDisplacement(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.Nobj

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
      #fig.savefig("./figures/"+self.nameCatalog+"/velocities/hist_disptheta.pdf")



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
      histGaussFit *= self.Nobj

      # histogram for LCDM
      sigma2 = self.U.RMSDisplacement(0., z0, W3d_sth)**2 / 3.
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histLCDM = np.array(map( g, np.linspace(0, nbins-2, nbins-1) ))
      histLCDM *= self.Nobj

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
      #fig.savefig("./figures/"+self.nameCatalog+"/velocities/hist_dispphi.pdf")


      plt.show()