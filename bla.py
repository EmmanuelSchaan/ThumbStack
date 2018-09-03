# M_star from M_vir, accounting for scatter
self.log10M1 = 11.43
self.log10e = -1.66
self.alpha = 1.750
self.delta = 4.290
self.gamma = 0.595


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
   result += self.(np.log10(mVir) - self.log10M1)
   result -= self.f(0.)
   result = 10.**result
   return result
