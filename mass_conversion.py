from headers import *

##################################################################################
##################################################################################

class MassConversionKravtsov14(object):
   """Conversions between stellar mass of the central galaxy
   and halo mass of the parent halo.
   From Kravtsov+14.
   """

   def __init__(self):
      
      # Params for M_star from M_vir,
      # accounting for scatter,
      # from Table 3 in Kravtsov+14.
      self.log10M1 = 11.39
      self.log10e = -1.685
      self.alpha = -1.740   # minus sign missing in Table 3
      self.delta = 4.335
      self.gamma = 0.531
      
      '''
      # Params for M_star from M_200m,
      # accounting for scatter,
      # from Table 3 in Kravtsov+14.
      self.log10M1 = 11.41
      self.log10e = -1.720
      self.alpha = -1.727   # minus sign missing in Table 3
      self.delta = 4.305
      self.gamma = 0.544
      '''

      # interpolate M_star = f(M_vir)
      self.mVir = np.logspace(np.log10(1.e9), np.log10(1.e18), 501, 10.) # [M_sun]
      self.mStar = np.array(map(self.fmStar, self.mVir))  # [M_sun]
      self.fmStarTomVir = interp1d(self.mStar, self.mVir, kind='cubic', bounds_error=True, fill_value=0.)
      self.fmVirTomStar = interp1d(self.mVir, self.mStar, kind='cubic', bounds_error=True, fill_value=0.)


   def f(self, x):
      """Fitting function, Kravtsov+14, eq A4
      """
      result = np.log10(1.+np.exp(x))**self.gamma
      result *= self.delta
      result /= 1. + np.exp(10.**(-x))
      result += -np.log10(10.**(self.alpha*x) + 1.)
      return result

   def fmStar(self, mVir):
      """Computes stellar mass [M_sun]
      from halo mass Mvir [Msun].
      """
      result = self.log10e + self.log10M1
      result += self.f(np.log10(mVir) - self.log10M1)
      result -= self.f(0.)
      result = 10.**result
      return result

   ##################################################################################

   def plot(self):
      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      ax.loglog(self.mVir, self.mStar, '.')
      #
      ax.set_xlabel(r'$M_\text{vir}$ [$M_\odot$]')
      ax.set_ylabel(r'$M_\star$ [$M_\odot$]')

      plt.show()






