from headers import *

##################################################################################


class MassConversion(object):
   
   def __init__(self):
      
      self.floadPC11()
      self.floadKravtsov2014()
      
      # choose Kravtsov et al 2014 as the default
      self.fMstellarToMvir = self.fMstellarToMvirKravtsov2014
      self.fMvirToMstellar = self.fMvirToMstellarKravtsov2014
      
#      # choose Planck Intermediate XI as the default
#      self.fMstellarToMvir = self.fMstellarToMvirPC11
#      self.fMvirToMstellar = self.fMvirToMstellarPC11

   ##################################################################################

   def floadPC11(self):
      # from Planck Intermediate 11, fig3
      # should be M200 and not Mvir...
      # watch that data is log10 of mass and not mass
      data = np.genfromtxt("./input/mass_conversion/digitizing_PCI11_fig3/PCI11_fig3.txt")
      
      ffMstellarToMvir = UnivariateSpline(data[:,0], data[:,1], k=3,s=0)
      mStellarMin = np.min(data[:,0])
      mStellarMax = np.max(data[:,0])
      self.fMstellarToMvirPC11 = lambda mStellar: ffMstellarToMvir(mStellar) * (mStellar>mStellarMin and mStellar<mStellarMax)

      ffMvirToMstellar = UnivariateSpline(data[:,1], data[:,0], k=3,s=0)
      mVirMin = np.min(data[:,1])
      mVirMax = np.max(data[:,1])
      self.fMvirToMstellarPC11 = lambda mVir: ffMvirToMstellar(mVir) * (mVir>mVirMin and mVir<mVirMax)


   def floadKravtsov2014(self):
      # from Kravtsov et al 2014, fig17
      # should be M200 and not Mvir
      data = np.genfromtxt("./input/mass_conversion/digitizing_kravtsov2014_fig17/kravtsov2014_fig17.txt")

      ffMstellarToMvir = UnivariateSpline(data[:,0], data[:,1], k=3,s=0)
      mStellarMin = np.min(data[:,0])
      mStellarMax = np.max(data[:,0])
      self.fMstellarToMvirKravtsov2014 = lambda mStellar: ffMstellarToMvir(mStellar) * (mStellar>mStellarMin and mStellar<mStellarMax)

      ffMvirToMstellar = UnivariateSpline(data[:,1], data[:,0], k=3,s=0)
      mVirMin = np.min(data[:,1])
      mVirMax = np.max(data[:,1])
      self.fMvirToMstellarKravtsov2014 = lambda mVir: ffMvirToMstellar(mVir) * (mVir>mVirMin and mVir<mVirMax)


   ##################################################################################


   def ftestInterpPC11(self):
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # raw array
      data = np.genfromtxt("./input//mass_conversion/digitizing_PCI11_fig3/PCI11_fig3.txt")
      ax.loglog(data[:,0], data[:,1], 'ko', label=r'Planck Intermediate 11, fig3')
      #
      # interpolation mStellar -> mVir
      Mstellar = np.logspace(np.log10(1.e10), np.log10(1.e12), 1001, 10.)
      MvirPC11 = np.array(map(self.fMstellarToMvirPC11, Mstellar))
      ax.loglog(Mstellar, MvirPC11, 'b', lw=2, label=r'interpolation $m_{*} \rightarrow m_\text{vir}$')
      #
      # interpolation mVir -> mStellar
      Mvir = np.logspace(np.log10(1.e10), np.log10(1.e15), 1001, 10.)
      MstellarPC11 = np.array(map(self.fMvirToMstellarPC11, Mvir))
      ax.loglog(MstellarPC11, Mvir, 'r--', lw=2, label=r'interpolation $m_\text{vir} \rightarrow m_{*}$')
      #
      ax.legend(loc=2, numpoints=1)
      ax.grid()
      ax.set_xlabel(r'$M_*$ [M$_\odot$]')
      ax.set_ylabel(r'$M_{200}$ [M$_\odot$]')
      #fig.savefig("./figures/mass_conversion/interpPC11.pdf")

      plt.show()


   def ftestInterpKravtsov2014(self):
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      data = np.genfromtxt("./input//mass_conversion/digitizing_kravtsov2014_fig17/kravtsov2014_fig17.txt")
      ax.loglog(data[:,0], data[:,1], 'ko', label=r'Kravtsov et al 2014, fig17')
      #
      # interpolation mStellar -> mVir
      Mstellar = np.logspace(np.log10(1.e7), np.log10(1.e13), 1001, 10.)
      MvirKravtsov2014 = np.array(map(self.fMstellarToMvirKravtsov2014, Mstellar))
      ax.loglog(Mstellar, MvirKravtsov2014, 'b', lw=2, label=r'interpolation $m_{*} \rightarrow m_\text{vir}$')
      #
      # interpolation mVir -> mStellar
      Mvir = np.logspace(np.log10(1.e10), np.log10(1.e16), 1001, 10.)
      MstellarKravtsov2014 = np.array(map(self.fMvirToMstellarKravtsov2014, Mvir))
      ax.loglog(MstellarKravtsov2014, Mvir, 'r--', lw=2, label=r'interpolation $m_\text{vir} \rightarrow m_{*}$')
      #
      ax.legend(loc=2, numpoints=1)
      ax.grid()
      ax.set_xlabel(r'$M_*$ [M$_\odot$]')
      ax.set_ylabel(r'$M_{200}$ [M$_\odot$]')
      #fig.savefig("./figures/mass_conversion/interpKravtsov2014.pdf")
      
      plt.show()



   def fcomparePC11Kravtsov2014(self):
      Mstellar = np.logspace(np.log10(1.e10), np.log10(1.e12), 1001, 10.)
      MvirPC11 = np.array(map(self.fMstellarToMvirPC11, Mstellar))
      MvirKravtsov2014 = np.array(map(self.fMstellarToMvirKravtsov2014, Mstellar))


      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.loglog(Mstellar, MvirPC11, 'b', lw=2, label=r'Planck Intermediate 11, fig3')
      ax.loglog(Mstellar, MvirKravtsov2014, 'r', lw=2, label=r'Kravtsov et al 2014, fig17')
      #
      ax.legend(loc=2)
      ax.grid()
      ax.set_xlabel(r'$M_*$ [M$_\odot$]')
      ax.set_ylabel(r'$M_{200}$ [M$_\odot$]')
      #fig.savefig("./figures/mass_conversion/comparePC11Kravtsov2014.pdf")


      fig=plt.figure(1)
      ax=fig.add_subplot(111)
      #
      ax.semilogx(Mstellar, MvirPC11/MvirKravtsov2014, 'b', lw=2, label=r'PC11 / Kravtsov2014')
      #
      ax.legend(loc=2)
      ax.grid()
      ax.set_xlabel(r'$M_*$ [M$_\odot$]')
      ax.set_ylabel(r'ratio of $M_{200}$ for PC11 and Kravtsov2014')
      #fig.savefig("./figures/mass_conversion/ratioPC11Kravtsov2014.pdf")

      plt.show()



   def plotStellarToVirRatio(self):
      Mstellar = np.logspace(np.log10(1.e9), np.log10(1.e12), 1001, 10.)
      MvirPC11 = np.array(map(self.fMstellarToMvirPC11, Mstellar))
      MvirKravtsov2014 = np.array(map(self.fMstellarToMvirKravtsov2014, Mstellar))
      
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.semilogx(MvirPC11, Mstellar/MvirPC11,  'b', lw=2, label=r'Planck Intermediate 11, fig3')
      ax.semilogx(MvirKravtsov2014, Mstellar/MvirKravtsov2014, 'r', lw=2, label=r'Kravtsov et al 2014, fig17')
      #
      ax.legend(loc=1)
      ax.grid()
      ax.set_ylabel(r'$M_*/M_{200}$')
      ax.set_xlabel(r'$M_{200}$ [M$_\odot$]')
      #fig.savefig("./figures/mass_conversion/stellar_to_vir_ratio.pdf")
      
      plt.show()



