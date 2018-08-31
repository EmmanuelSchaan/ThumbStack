from headers import *

##################################################################################


class Catalog(object):
   
   def __init__(self, U, MassConversion):
 
       self.U = U
       self.MassConversion = MassConversion

      # galaxy or cluster catalog name
      #self.name = "cmasssouth"
      #self.nameLong = "CMASS South"
      
      # path to the complete catalog with halo masses
      #self.path = "./output/"+self.name+"/objects.txt"


   ##################################################################################
   
   # expected std dev of velocities for LCDM
   # comparison between using median z or z-distribution
   def ftestVel(self):
      z0 = np.mean(self.Z)
      
      f = lambda z: self.U.RMSVelocity(0., z, W3d_sth)**2 / 3.
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
      histGaussFit *= self.Nobj
      
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
##################################################################################

class CMASSSouth(Catalog):
   
   def __init__(self, U, MassConversion, save=False):
      
      super(CMASSSouth, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "cmasssouth"
      self.nameLong = "CMASS South"
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"

      # add halo masses to the catalog if necessary
      self.pathInput = "./input/cmasssouth_dr12.txt"
      self.floadCatalog(save=save)
   
   
   ##################################################################################
   
   def faddHaloMass(self):
      data = np.genfromtxt(self.pathInput)
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
            Mvir[iobj] = self.MassConversion.fMstellarToMvir(mStellar)
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      
      newData = np.zeros((nObjRemaining, 14))
      newData[:,:13] = data[iRemaining,:13]
      newData[:,13] = Mvir[iRemaining]
      print "- save catalog with halo masses to "+self.path
      np.savetxt(self.path, newData)
   
   
   def floadCatalog(self, save=False):
      
      # add halo mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddHaloMass()
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
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
      # just for test purposes, to plot their histograms and compare it to the ones of velocities
      dispData = np.genfromtxt("../../data/CMASS_DR12/displacements/dr_dtheta_dphi_DR12_CMASS_S.txt")
      # get the forward displacements
      # watch out that dispR has the redshift-space correction in it
      self.dispR = -dispData[:,0]
      self.dispTheta = -dispData[:,1]
      self.dispPhi = -dispData[:,2]


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
      #fig.savefig("./figures/"+self.name+"/velocities/hist_dispphi.pdf")
      
                              
      plt.show()


##################################################################################
##################################################################################

#class Random(Catalog):
#   
#   def __init__(self, U, MassConversion, name="random", save=False):
#      
#      super(CMASSSouth, self).__init__(U, MassConversion)
#      
#      # galaxy or cluster catalog
#      self.name = name
#      self.nameLong = self.name
#      
#      # path to the catalog with halo masses
#      self.path = "./output/"+self.name+"/objects.txt"
#      
#      # add halo masses to the catalog if necessary
#      self.pathInput = "./input/cmasssouth_dr12.txt"
#      self.floadCatalog(save=save)
#
#
#   ##################################################################################
#
#   def floadCatalog(self, save=False):
#      nTot = 1.e4
#      catalogRA = np.random.uniform(ra_min, ra_max, size=nTot)
#      catalogDEC = np.random.uniform(dec_min, dec_max, size=nTot)
#      std_vr = self.U.RMSVelocity(0., 0.57, W3d_sth)/np.sqrt(3.)
#      catalogVelR = np.random.normal(loc=0.0, scale=std_vr, size=nTot)
#      self.Z = 0.57*np.ones(nTot)
#      
#      data = np.genfromtxt(self.path)
#      self.coordX = np.zeros(nTot)
#      self.coordY = np.zeros(nTot)
#      self.coordZ = np.zeros(nTot)
#
#      self.velX = np.zeros(nTot)
#      self.velY = np.zeros(nTot)
#      self.velZ = np.zeros(nTot)
#      self.velTheta = np.zeros(nTot)
#      self.velPhi = np.zeros(nTot)
#      self.Mstellar = 2.e11*np.ones(nTot)
#      self.Mvir = self.MassConversion.fMstellarToMvir(2.e11)*np.ones(nTot)
#      self.Nobj = len(self.coordX)

##################################################################################
##################################################################################

class CMASSSouthSmoo(Catalog):
   
   def __init__(self, U, MassConversion, save=False, smoo=10, meanMass=False, rand=False):
      
      super(CMASSSouthSmoo, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "cmasssouth_smoo"+str(smoo)
      if rand==True:
         self.name += "randradec"
      self.nameLong = "CMASS South smoo"+str(smoo)
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      # add halo masses to the catalog if necessary
      self.pathInput = "./input/cmasssouth_dr12_smoo"+str(smoo)+".txt"
      self.floadCatalog(save=save, meanMass=meanMass, rand=rand)
   
   
   ##################################################################################
   
   def faddHaloMass(self, rand=False):
      data = np.genfromtxt(self.pathInput)
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
            Mvir[iobj] = self.MassConversion.fMstellarToMvir(mStellar)
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      
      newData = np.zeros((nObjRemaining, 14))
      newData[:,:13] = data[iRemaining,:13]
      newData[:,13] = Mvir[iRemaining]
      if rand==True:
         print "adding random offset to RA and DEC"
         #np.random.shuffle(newData[:,3])  # shuffle RA
         #np.random.shuffle(newData[:,4])  # shuffle DEC
         Theta = np.random.uniform(low=0., high=2.*np.pi, size=len(newData[:,0]))
         offset = 15./60.  # degrees
         newData[:,3] += offset * np.cos(Theta)
         newData[:,4] += offset * np.sin(Theta)
      print "- save catalog with halo masses to "+self.path
      np.savetxt(self.path, newData)
   
   
   def floadCatalog(self, save=False, meanMass=False, rand=False):
      
      # add halo mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddHaloMass(rand=rand)
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
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

      if meanMass:
         self.Mvir[:] = 2.e13/self.U.h # 2.e13 Msun/h from Miyatake et al 2013


##################################################################################
##################################################################################

class MockCMASSSouthACTDeep5656(Catalog):
   
   def __init__(self, U, MassConversion, save=False, shuffledVrec=False, meanMass=False):
      
      super(MockCMASSSouthACTDeep5656, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "mockcmasssouth_actdeep5656"
      if shuffledVrec:
         self.name += "_shuffledvrec"
      if meanMass:
         self.name += "_meanmass"
      self.nameLong = self.name
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      # add halo masses to the catalog if necessary
      self.pathInput = "../../data/mocks_actdeep5656_Simo/20march2015/Deep56_coords_mock_map_actbeam_Ncluster22332_z0.57_vrms0.00103_c5.0.txt"
      self.floadCatalog(save=save, shuffledVrec=shuffledVrec, meanMass=meanMass)
   
   
   ##################################################################################
   
   def faddStellarMass(self, shuffledVrec=False, meanMass=True):
      data = np.genfromtxt(self.pathInput)
      #RA = data[:,0]
      #DEC = data[:,1]
      #velR = data[:,2] * 3.e5 # speed in km/s
      mVir = data[:,3]
      Nobj = len(mVir)
      
      # compute corresponding stellar masses
      Flag = np.ones(Nobj)
      mStellar = np.zeros(Nobj)
      for iobj in range(Nobj):
         mvir = mVir[iobj]
         mstellar = self.MassConversion.fMvirToMstellar(mvir)
         if (mstellar<1.) or np.isnan(mstellar):
            Flag[iobj] = 0.
         else:
            mStellar[iobj] = mstellar
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      newData = np.zeros((nObjRemaining, 5))
      newData[:,:4] = data[iRemaining,:4]
      newData[:,4] = mStellar[iRemaining]
      
      if shuffledVrec:
         # shuffle velR = data[:,9]
         np.random.shuffle(newData[:,9])
      
      if meanMass:
         # replace Mstellar and Mvir by their means
         newData[:,3] = np.mean(newData[:,3]) * np.ones(nObjRemaining)
         newData[:,4] = np.mean(newData[:,4]) * np.ones(nObjRemaining)

      
      print "- save catalog with halo masses to "+self.path
      np.savetxt(self.path, newData)
   
   
   def floadCatalog(self, save=False, shuffledVrec=False, meanMass=True):

      # add halo mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddStellarMass(shuffledVrec=shuffledVrec, meanMass=meanMass)
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
      self.Nobj = len(data[:,0])
      self.coordX = np.zeros(self.Nobj)
      self.coordY = np.zeros(self.Nobj)
      self.coordZ = np.zeros(self.Nobj)
      self.RA = data[:,0]
      self.DEC = data[:,1]
      self.Z = 0.57 * np.ones(self.Nobj)
      self.velX = np.zeros(self.Nobj)
      self.velY = np.zeros(self.Nobj)
      self.velZ = np.zeros(self.Nobj)
      self.velR = data[:,2]
      self.velTheta = np.zeros(self.Nobj)
      self.velPhi = np.zeros(self.Nobj)
      self.Mstellar = data[:,4]
      self.Mvir = data[:,3]


##################################################################################
##################################################################################

class BigMockCMASSSouthACTDeep5656(Catalog):
   
   def __init__(self, U, MassConversion, save=False, shuffledVrec=False, meanMass=False):
      
      super(BigMockCMASSSouthACTDeep5656, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "bigmockcmasssouth_actdeep5656"
      if shuffledVrec:
         self.name += "_shuffledvrec"
      if meanMass:
         self.name += "_meanmass"
      self.nameLong = self.name
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      # add halo masses to the catalog if necessary
      self.pathInput = "../../data/mocks_actdeep5656_Simo/20march2015/Deep56_coords_mock_map_actbeam_Ncluster48014_z0.57_vrms0.00103_c5.0.txt"
      self.floadCatalog(save=save, shuffledVrec=shuffledVrec, meanMass=meanMass)
   
   
   ##################################################################################
   
   def faddStellarMass(self, shuffledVrec=False, meanMass=True):
      data = np.genfromtxt(self.pathInput)
      #RA = data[:,0]
      #DEC = data[:,1]
      #velR = data[:,2] * 3.e5 # speed in km/s
      mVir = data[:,3]
      Nobj = len(mVir)
      
      # compute corresponding stellar masses
      Flag = np.ones(Nobj)
      mStellar = np.zeros(Nobj)
      for iobj in range(Nobj):
         mvir = mVir[iobj]
         mstellar = self.MassConversion.fMvirToMstellar(mvir)
         if (mstellar<1.) or np.isnan(mstellar):
            Flag[iobj] = 0.
         else:
            mStellar[iobj] = mstellar
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      newData = np.zeros((nObjRemaining, 5))
      newData[:,:4] = data[iRemaining,:4]
      newData[:,4] = mStellar[iRemaining]
      
      if shuffledVrec:
         # shuffle velR = data[:,9]
         np.random.shuffle(newData[:,9])
      
      if meanMass:
         # replace Mstellar and Mvir by their means
         newData[:,3] = np.mean(newData[:,3]) * np.ones(nObjRemaining)
         newData[:,4] = np.mean(newData[:,4]) * np.ones(nObjRemaining)
      
      
      print "- save catalog with halo masses to "+self.path
      np.savetxt(self.path, newData)
   
   
   def floadCatalog(self, save=False, shuffledVrec=False, meanMass=True):
      
      # add halo mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddStellarMass(shuffledVrec=shuffledVrec, meanMass=meanMass)
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
      self.Nobj = len(data[:,0])
      self.coordX = np.zeros(self.Nobj)
      self.coordY = np.zeros(self.Nobj)
      self.coordZ = np.zeros(self.Nobj)
      self.RA = data[:,0]
      self.DEC = data[:,1]
      self.Z = 0.57 * np.ones(self.Nobj)
      self.velX = np.zeros(self.Nobj)
      self.velY = np.zeros(self.Nobj)
      self.velZ = np.zeros(self.Nobj)
      self.velR = data[:,2]
      self.velTheta = np.zeros(self.Nobj)
      self.velPhi = np.zeros(self.Nobj)
      self.Mstellar = data[:,4]
      self.Mvir = data[:,3]



##################################################################################
##################################################################################

#class CMASSSouthDR10(Catalog):
#   
#   def __init__(self, U, MassConversion, meanMass=False, save=False):
#      
#      super(CMASSSouthDR10, self).__init__(U, MassConversion)
#      
#      # galaxy or cluster catalog
#      self.name = "cmasssouth_dr10"
#      self.nameLong = "CMASS South DR10"
#      
#      # path to the catalog
#      self.path = "./output/"+self.name+"/objects.txt"
#      
#      # add halo masses to the catalog if necessary
#      self.pathInput = "./input/cmasssouth_dr10.txt"
#      
#      # path to the catalog
#      self.path = self.pathInput
#      
#      self.floadCatalog(meanMass=meanMass, save=save)
#   
#   
#   ##################################################################################
#
#   def floadCatalog(self, meanMass=False, save=False):
#      
#      # add halo mass to catalog if not already done
#      #if (not os.path.isfile(self.path)) or save:
#      #   self.faddHaloMass()
#      
#      print "- read catalog with halo masses from "+self.path
#      data = np.genfromtxt(self.path)
#      Nobj = len(data[:,0])
#      self.coordX = np.zeros(Nobj)
#      self.coordY = np.zeros(Nobj)
#      self.coordZ = np.zeros(Nobj)
#      self.RA = data[:,0]
#      self.DEC = data[:,1]
#      self.Z = 0.57 * np.ones(Nobj)
#      self.velX = np.zeros(Nobj)
#      self.velY = np.zeros(Nobj)
#      self.velZ = np.zeros(Nobj)
#      self.velR = np.zeros(Nobj)
#      self.velTheta = np.zeros(Nobj)
#      self.velPhi = np.zeros(Nobj)
#      self.Mstellar = np.zeros(Nobj)
#      self.Mvir = np.zeros(Nobj)
#      self.Nobj = Nobj
#
#      if meanMass:
#         self.Mvir[:] = 2.e13/self.U.h # 2.e13 Msun/h from Miyatake et al 2013


##################################################################################
##################################################################################

class CMASSSouthSmooMeanMassFull(Catalog):
   
   def __init__(self, U, MassConversion, save=False, smoo=10):
      
      super(CMASSSouthSmooMeanMassFull, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "cmasssouth_smoo"+str(smoo)+"_meanmass_full"
      self.nameLong = "CMASS South smoo"+str(smoo)
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      # add halo masses to the catalog if necessary
      self.pathInput = "./input/cmasssouth_dr12_smoo"+str(smoo)+".txt"
      self.floadCatalog(save=save)
   
   
   ##################################################################################
   
   def faddHaloMassAndShuffleVrec(self):
      data = np.genfromtxt(self.pathInput)
      Mstellar = data[:,12]
      Nobj = len(Mstellar)
      '''
      # shuffle velR = data[:,9]
      np.random.shuffle(data[:,9])
      
      # keep all objects with a known stellar mass
      Flag = np.ones(Nobj)
      Mvir = np.zeros(Nobj)
      for iobj in range(Nobj):
         mStellar = Mstellar[iobj]
         if (mStellar<1.) or np.isnan(mStellar):
            Flag[iobj] = 0.
         else:
            Mvir[iobj] = self.MassConversion.fMstellarToMvir(mStellar)
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      '''
      newData = np.zeros((Nobj, 14))
      newData[:,:12] = data[:,:12]
      newData[:,12] = 2.420e11*np.ones(Nobj)
      newData[:,13] = 2.104e13*np.ones(Nobj)
      print "- save catalog with halo masses to "+self.path
      np.savetxt(self.path, newData)
   
   
   def floadCatalog(self, save=False):
      
      # add halo mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddHaloMassAndShuffleVrec()
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
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



##################################################################################
##################################################################################

class CMASSSouthDR10Kendrick(Catalog):
   
   def __init__(self, U, MassConversion, meanMass=False, save=False):
      
      super(CMASSSouthDR10Kendrick, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "cmasssouth_dr10_kendrick"
      self.nameLong = "CMASS South DR10"
      
      # input catalog (without halo masses)
      self.pathInput = "./input/cmasssouth_dr10_kendrick.txt"
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      self.floadCatalog(meanMass=False, save=save)

   
   
   ##################################################################################
   
   
   def faddHaloMass(self):
      data = np.genfromtxt(self.pathInput)
      Mstellar = data[:,4]
      Nobj = len(Mstellar)
      
      # keep all objects with a known stellar mass
      Flag = np.ones(Nobj)
      Mvir = np.zeros(Nobj)
      for iobj in range(Nobj):
         mStellar = Mstellar[iobj]
         if (mStellar<1.) or np.isnan(mStellar):
            Flag[iobj] = 0.
         else:
            Mvir[iobj] = self.MassConversion.fMstellarToMvir(mStellar)
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      
      newData = np.zeros((nObjRemaining, 6))
      newData[:,:5] = data[iRemaining,:5]
      newData[:,5] = Mvir[iRemaining]
      print "- save catalog with halo masses to "+self.path
      #print "initially nobj = ", Nobj, ". Remaining", nObjRemaining
      np.savetxt(self.path, newData)
   
   
   
   def floadCatalog(self, meanMass=False, save=False):
      
      # add halo mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddHaloMass()
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
      Nobj = len(data[:,0])
      self.coordX = np.zeros(Nobj)
      self.coordY = np.zeros(Nobj)
      self.coordZ = np.zeros(Nobj)
      self.RA = data[:,0]
      self.DEC = data[:,1]
      self.Z = data[:,2]
      self.velX = np.zeros(Nobj)
      self.velY = np.zeros(Nobj)
      self.velZ = np.zeros(Nobj)
      self.velR = data[:,3]
      self.velTheta = np.zeros(Nobj)
      self.velPhi = np.zeros(Nobj)
      self.Mstellar = data[:,4]
      self.Mvir = data[:,5]
      self.Nobj = Nobj

      if meanMass:
         self.Mvir[:] = 2.e13/self.U.h # 2.e13 Msun/h from Miyatake et al 2013



##################################################################################
##################################################################################

class LOWZSouthMariana(Catalog):
   
   def __init__(self, U, MassConversion, meanMass=False, save=False):
      
      super(LOWZSouthMariana, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "lowzsouth_dr12"
      self.nameLong = "LOWZ South DR12"
      
      # input catalog (without halo masses)
      #self.pathInput = "./input/lowzsouth_dr10_kendrick.txt"
      self.pathInput = "../prepare_cmass_mariana/output/lowzsouth_dr12.txt"
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      self.floadCatalog(meanMass=False, save=save)
   
   
   
   ##################################################################################
   
   
   def faddHaloMass(self, rand=False, shift=False, raShift=0, decShift=0, parity=0):
      data = np.genfromtxt(self.pathInput)
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
            Mvir[iobj] = self.MassConversion.fMstellarToMvir(mStellar)
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      
      newData = np.zeros((nObjRemaining, 14))
      newData[:,:13] = data[iRemaining,:13]
      newData[:,13] = Mvir[iRemaining]
      
      if rand:
         # rough bounds of ACT Deep5656
         ra_min = -15.
         ra_max = 46.
         dec_min = -10.
         dec_max = 7.
         # generate random RA and DEC
         newData[:,3] = np.random.uniform(ra_min, ra_max, size=nObjRemaining)  # ra
         newData[:,4] = np.random.uniform(dec_min, dec_max, size=nObjRemaining)   # dec
      if shift:
         # shift RA and DEC
         newData[:,3] += raShift
         newData[:,4] += decShift
      if parity==1:
         newData[:,3] = 38.5 - newData[:,0]
      if parity==2:
         newData[:,4] = -1.5 - newData[:,1]
      if parity==3:
         newData[:,3] = 38.5 - newData[:,0]
         newData[:,4] = -1.5 - newData[:,1]
      
      print "- save catalog with halo masses to "+self.path
      #print "initially nobj = ", Nobj, ". Remaining", nObjRemaining
      np.savetxt(self.path, newData)


   def floadCatalog(self, meanMass=False, rand=False, shift=False, raShift=0, decShift=0, parity=0, save=False):
      
      # add halo mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddHaloMass(rand=rand, shift=False, raShift=raShift, decShift=decShift, parity=parity)
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
      Nobj = len(data[:,0])
      self.coordX = data[:,0]
      self.coordY = data[:,1]
      self.coordZ = data[:,2]
      self.RA = data[:,3]
      self.DEC = data[:,4]
      self.Z = data[:,5]
      self.velX = np.zeros(Nobj)
      self.velY = np.zeros(Nobj)
      self.velZ = np.zeros(Nobj)
      #
      self.velR = data[:,9]
      self.velTheta = data[:,10]
      self.velPhi = data[:,11]
      self.Mstellar = data[:,12] #self.MassConversion.fMvirToMstellar(2.e13/self.U.h) * np.ones(Nobj)
      self.Mvir = data[:,13]  #2.e13/self.U.h * np.ones(Nobj) # 2.e13 Msun/h from Miyatake et al 2013
      self.Nobj = len(self.RA)

      if meanMass:
         self.Mvir[:] = 2.e13/self.U.h # 2.e13 Msun/h from Miyatake et al 2013



##################################################################################
##################################################################################

class LOWZSouthDR10Kendrick(Catalog):
   
   def __init__(self, U, MassConversion, meanMass=False, save=False):
      
      super(LOWZSouthDR10Kendrick, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "lowzsouth_dr10_kendrick"
      self.nameLong = "LOWZ South DR10"
      
      # input catalog (without halo masses)
      #self.pathInput = "./input/lowzsouth_dr10_kendrick.txt"
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      self.floadCatalog(meanMass=False, save=save)
   
   
   
   ##################################################################################

   def floadCatalog(self, meanMass=False, save=False):
      
      print "- read DR10 catalog"

      path = "../../data/CMASS_DR10/vrec_kendrick/v1/galaxy_DR10v8_LOWZ_South_vrad_using_randoms.fits"
      hdu_number = 0 # 0 for the primary hdu
      catalog = pyfits.getdata(path, hdu_number)
      #collateIndx = catalog[:]['collateindx']
      
      Nobj = len(catalog[:]['z'])
      self.coordX = np.zeros(Nobj)
      self.coordY = np.zeros(Nobj)
      self.coordZ = np.zeros(Nobj)
      self.RA = catalog[:]['ra_deg']
      self.DEC = catalog[:]['dec_deg']
      self.Z = catalog[:]['z']
      self.velX = np.zeros(Nobj)
      self.velY = np.zeros(Nobj)
      self.velZ = np.zeros(Nobj)
      #
      self.velR = catalog[:]['vr'] * 3.e5
      #
      self.velTheta = np.zeros(Nobj)
      self.velPhi = np.zeros(Nobj)
      self.Mstellar = self.MassConversion.fMvirToMstellar(2.e13/self.U.h) * np.ones(Nobj)
      self.Mvir = 2.e13/self.U.h * np.ones(Nobj) # 2.e13 Msun/h from Miyatake et al 2013
      self.Nobj = len(self.RA)























##################################################################################
##################################################################################

class CMASSSouthMarianaKendrick(Catalog):
   
   def __init__(self, U, MassConversion, meanMass=False, rand=False, randId=0, shift=False, raShift=0, decShift=0, parity=0, save=False):
      
      super(CMASSSouthMarianaKendrick, self).__init__(U, MassConversion)
      
      # input catalog (without halo masses)
      self.nameLong = "CMASS South"
      self.pathInput = "./input/cmasssouth_mariana_kendrick.txt"
      
      # galaxy or cluster catalog
      self.name = "cmasssouth_mariana_kendrick"
      if rand:
         self.name = "cmasssouth_mariana_kendrick/random"+str(randId)
      if shift:
         self.name = "cmasssouth_mariana_kendrick/shift_ra"+str(raShift)+"_dec"+str(decShift)+"_parity"+str(parity)
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      
      self.floadCatalog(meanMass=meanMass, rand=rand, shift=shift, raShift=raShift, decShift=decShift, parity=parity, save=save)
   
   
   
   ##################################################################################
   
   
   def faddHaloMass(self, rand=False, shift=False, raShift=0, decShift=0, parity=0):
      data = np.genfromtxt(self.pathInput)
      Mstellar = data[:,7]
      Nobj = len(Mstellar)
      
      # keep all objects with a known stellar mass
      Flag = np.ones(Nobj)
      Mvir = np.zeros(Nobj)
      for iobj in range(Nobj):
         mStellar = Mstellar[iobj]
         if (mStellar<1.) or np.isnan(mStellar):
            Flag[iobj] = 0.
         else:
            Mvir[iobj] = self.MassConversion.fMstellarToMvir(mStellar)
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      #print "Nobj w ok masses=", nObjRemaining
      
      newData = np.zeros((nObjRemaining, 9))
      newData[:,:8] = data[iRemaining,:8]
      newData[:,8] = Mvir[iRemaining]
      
      if rand:
         # rough bounds of ACT Deep5656
         ra_min = -15.
         ra_max = 46.
         dec_min = -10.
         dec_max = 7.
         # generate random RA and DEC
         newData[:,0] = np.random.uniform(ra_min, ra_max, size=nObjRemaining)  # ra
         newData[:,1] = np.random.uniform(dec_min, dec_max, size=nObjRemaining)   # dec
      if shift:
         # shift RA and DEC
         newData[:,0] += raShift
         newData[:,1] += decShift
      if parity==1:
         newData[:,0] = 38.5 - newData[:,0]
      if parity==2:
         newData[:,1] = -1.5 - newData[:,1]
      if parity==3:
         newData[:,0] = 38.5 - newData[:,0]
         newData[:,1] = -1.5 - newData[:,1]
      
      print "- save catalog with halo masses to "+self.path
      #print "initially nobj = ", Nobj, ". Remaining", nObjRemaining
      np.savetxt(self.path, newData)
   
   
   
   def floadCatalog(self, meanMass=False, rand=False, shift=False, raShift=0, decShift=0, parity=0, save=False):
      
      # add halo mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddHaloMass(rand=rand, shift=shift, raShift=raShift, decShift=decShift, parity=parity)
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
      Nobj = len(data[:,0])
      self.coordX = np.zeros(Nobj)
      self.coordY = np.zeros(Nobj)
      self.coordZ = np.zeros(Nobj)
      self.RA = data[:,0]
      self.DEC = data[:,1]
      self.Z = data[:,2]
      self.velX = np.zeros(Nobj)
      self.velY = np.zeros(Nobj)
      self.velZ = np.zeros(Nobj)
      #
      self.velR = data[:,3]   # Mariana's velocities
      #self.velR = data[:,6]   # Kendrick's velocities
      self.velTheta = data[:,4]  # Mariana
      self.velPhi = data[:,5]  # Mariana
      #
      self.Mstellar = data[:,7]
      self.Mvir = data[:,8]
      self.Nobj = Nobj
      
      # for nulls
      #np.random.shuffle(self.velR)
      #self.velR = np.copy(self.velTheta)

      if meanMass:
         self.Mvir[:] = 2.e13/self.U.h # 2.e13 Msun/h from Miyatake et al 2013



##################################################################################
##################################################################################

class MockCMASSSouthDR10Kendrick(Catalog):
   
   def __init__(self, U, MassConversion, id=0, offCentering=0, shuffledV=False, save=False):
      
      super(MockCMASSSouthDR10Kendrick, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "mockcmasssouth_dr10_kendrick/mock_"+str(id)
      self.nameLong = "mock CMASS South DR10 "+str(id)
      if offCentering<>0:
         self.name += "_off"+"{:.1f}".format(offCentering)
         self.nameLong += " off "+ "{:.1f}".format(offCentering) +" arcmin"
      if shuffledV:
         self.name += "_shuffledv"
      
      # input catalog (without halo masses)
      self.pathInput = "../../data/mocks_kendrick_actdeep5656/Deep56_coords_CMASS_South_Kendrick_MOCKNUM"+str(id)+"_c5.0.dat"
      if shuffledV:
         self.pathInput = "../../data/mocks_kendrick_actdeep5656/randomv_Deep56_coords_CMASS_South_Kendrick_MOCKNUM"+str(id)+"_c5.0.dat"
      
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      self.floadCatalog(save=save, offCentering=offCentering)
   
   
   
   ##################################################################################
   
   def faddStellarMass(self):
      data = np.genfromtxt(self.pathInput)
      mVir = data[:,5]
      Nobj = len(mVir)
      
      # compute corresponding stellar masses
      Flag = np.ones(Nobj)
      mStellar = np.zeros(Nobj)
      for iobj in range(Nobj):
         mvir = mVir[iobj]
         mstellar = self.MassConversion.fMvirToMstellar(mvir)
         if (mstellar<1.) or np.isnan(mstellar):
            Flag[iobj] = 0.
         else:
            mStellar[iobj] = mstellar
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      newData = np.zeros((nObjRemaining, 7))
      newData[:,:6] = data[iRemaining,:6]
      newData[:,6] = mStellar[iRemaining]
      
      print "- save catalog with halo masses to "+self.path
      np.savetxt(self.path, newData)

   
   
   
   def floadCatalog(self, save=False, offCentering=0):
      
      # add stellar mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddStellarMass()
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
      Nobj = len(data[:,0])
      self.coordX = np.zeros(Nobj)
      self.coordY = np.zeros(Nobj)
      self.coordZ = np.zeros(Nobj)
      
      self.RA = data[:,0]
      self.DEC = data[:,1]
      if offCentering<>0:
         self.RA += np.random.normal(0., offCentering/60., Nobj)
         self.DEC += np.random.normal(0., offCentering/60., Nobj)

      self.Z = data[:,4]   # 4 is z_obs, 3 is z_true
      self.velX = np.zeros(Nobj)
      self.velY = np.zeros(Nobj)
      self.velZ = np.zeros(Nobj)

      #self.velR = data[:,2] * 3.e5  # reconstructed velocity
      self.velR = 3.e5 * (data[:,4] - data[:,3])/(1. + data[:,3])  # vtrue = c*(z_obs-z_true)/(1+z_true)

      self.velTheta = np.zeros(Nobj)
      self.velPhi = np.zeros(Nobj)
      self.Mstellar = data[:,6]
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!! Extra factor of h !!!!!!!!!!!!!!!!!!!!!!!
      self.Mvir = data[:,5] / 0.7
      
      # add scatter from Planck XI
      '''
      print "- add scatter in mass"
      self.MassMultNoise = 10.**np.random.normal(0., np.log10(2.), Nobj)
      self.MassMultNoise /= np.mean(self.MassMultNoise)
      self.Mvir *= self.MassMultNoise
      '''
      
      # replace Mvir by its mean
      #self.Mvir = np.mean(self.Mvir) * np.ones(Nobj)
      
      self.Nobj = Nobj



##################################################################################
##################################################################################

class MockC10Clusters(Catalog):
   
   def __init__(self, U, MassConversion, save=False):
      
      super(MockC10Clusters, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "mockcmasssouth_dr10_kendrick/mock_10clusters"
      self.nameLong = "mock 10 clusters"
      
      # input catalog (without halo masses)
      #self.pathInput = "../../data/mocks_simo_10clusters/10clusters_normalized_coords.txt"
      self.pathInput = "../../data/mocks_simo_10clusters/randomv_10clusters_normalized_coords.txt"
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      self.floadCatalog(save=save)
   
   
   
   ##################################################################################
   
   def faddStellarMass(self):
      data = np.genfromtxt(self.pathInput)
      mVir = data[:,5]
      Nobj = len(mVir)
      
      # compute corresponding stellar masses
      Flag = np.ones(Nobj)
      mStellar = np.zeros(Nobj)
      for iobj in range(Nobj):
         mvir = mVir[iobj]
         mstellar = self.MassConversion.fMvirToMstellar(mvir)
         if (mstellar<1.) or np.isnan(mstellar):
            Flag[iobj] = 0.
         else:
            mStellar[iobj] = mstellar
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      newData = np.zeros((nObjRemaining, 7))
      newData[:,:6] = data[iRemaining,:6]
      newData[:,6] = mStellar[iRemaining]
      
      print "- save catalog with halo masses to "+self.path
      np.savetxt(self.path, newData)
   
   #(ra, dec, v_rec/c, ztrue, zobs, Mvir)
   
   
   def floadCatalog(self, save=False):
      
      # add stellar mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddStellarMass()
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
      Nobj = len(data[:,0])
      self.coordX = np.zeros(Nobj)
      self.coordY = np.zeros(Nobj)
      self.coordZ = np.zeros(Nobj)
      
      self.RA = data[:,0]
      self.DEC = data[:,1]
      
      self.Z = data[:,4]   # 4 is z_obs, 3 is z_true
      self.velX = np.zeros(Nobj)
      self.velY = np.zeros(Nobj)
      self.velZ = np.zeros(Nobj)
      
      #self.velR = data[:,2] * 3.e5  # reconstructed velocity
      self.velR = 3.e5 * (data[:,4] - data[:,3])/(1. + data[:,3])  # vtrue = c*(z_obs-z_true)/(1+z_true)
      
      self.velTheta = np.zeros(Nobj)
      self.velPhi = np.zeros(Nobj)
      self.Mstellar = data[:,6]
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!! Extra factor of h !!!!!!!!!!!!!!!!!!!!!!!
      #self.Mvir = data[:,5] / 0.7

      self.Nobj = Nobj




##################################################################################
##################################################################################

class CMASSNorthDR10Kendrick(Catalog):
   
   def __init__(self, U, MassConversion, save=False):
      
      super(CMASSNorthDR10Kendrick, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "cmassnorth_dr10_kendrick"
      self.nameLong = "CMASS North DR10"
      
      # input catalog (without halo masses)
      self.pathInput = "./input/cmassnorth_dr10_kendrick.txt"
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      self.floadCatalog(save=save)
   
   
   
   ##################################################################################
   
   
   def faddHaloMass(self):
      data = np.genfromtxt(self.pathInput)
      Mstellar = data[:,4]
      Nobj = len(Mstellar)
      
      # keep all objects with a known stellar mass
      Flag = np.ones(Nobj)
      Mvir = np.zeros(Nobj)
      for iobj in range(Nobj):
         mStellar = Mstellar[iobj]
         if (mStellar<1.) or np.isnan(mStellar):
            Flag[iobj] = 0.
         else:
            Mvir[iobj] = self.MassConversion.fMstellarToMvir(mStellar)
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      
      newData = np.zeros((nObjRemaining, 6))
      newData[:,:5] = data[iRemaining,:5]
      newData[:,5] = Mvir[iRemaining]
      print "- save catalog with halo masses to "+self.path
      #print "initially nobj = ", Nobj, ". Remaining", nObjRemaining
      np.savetxt(self.path, newData)
   
   
   
   def floadCatalog(self, save=False):
      
      # add halo mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddHaloMass()
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
      Nobj = len(data[:,0])
      self.coordX = np.zeros(Nobj)
      self.coordY = np.zeros(Nobj)
      self.coordZ = np.zeros(Nobj)
      self.RA = data[:,0]
      self.DEC = data[:,1]
      self.Z = data[:,2]
      self.velX = np.zeros(Nobj)
      self.velY = np.zeros(Nobj)
      self.velZ = np.zeros(Nobj)
      self.velR = data[:,3]
      self.velTheta = np.zeros(Nobj)
      self.velPhi = np.zeros(Nobj)
      self.Mstellar = data[:,4]
      self.Mvir = data[:,5]
      self.Nobj = Nobj


##################################################################################
##################################################################################

class CMASSNorthSmoo(Catalog):
   
   def __init__(self, U, MassConversion, save=False, smoo=10):
      
      super(CMASSNorthSmoo, self).__init__(U, MassConversion)
      
      # galaxy or cluster catalog
      self.name = "cmassnorth_smoo"+str(smoo)
      self.nameLong = "CMASS North smoo"+str(smoo)
      
      # path to the catalog with halo masses
      self.path = "./output/"+self.name+"/objects.txt"
      
      # add halo masses to the catalog if necessary
      self.pathInput = "./input/cmassnorth_dr12_smoo"+str(smoo)+".txt"
      self.floadCatalog(save=save)
   
   
   ##################################################################################
   
   def faddHaloMass(self):
      data = np.genfromtxt(self.pathInput)
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
            Mvir[iobj] = self.MassConversion.fMstellarToMvir(mStellar)
      
      iRemaining = np.where(Flag==1)[0]
      nObjRemaining = len(iRemaining)
      
      newData = np.zeros((nObjRemaining, 14))
      newData[:,:13] = data[iRemaining,:13]
      newData[:,13] = Mvir[iRemaining]
      print "- save catalog with halo masses to "+self.path
      np.savetxt(self.path, newData)
   
   
   def floadCatalog(self, save=False):
      
      # add halo mass to catalog if not already done
      if (not os.path.isfile(self.path)) or save:
         self.faddHaloMass()
      
      print "- read catalog with halo masses from "+self.path
      data = np.genfromtxt(self.path)
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







