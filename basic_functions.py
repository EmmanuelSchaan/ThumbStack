from headers import *

##################################################################################
#  Mathematical functions


def W3d_sth(x):
   """Fourier transform of a 3d spherical top hat window function.
   Use x = k*R as input,
   where R is the tophat radius and k the wave vector.
   Input and output are dimensionless.
   """
   if x < 1.e-3:  # for small x, replace by expansion, for numerical stability
      f = 1. - 0.1* x**2 + 0.00357143* x**4
   else:
      f = (3./(x**3)) * ( np.sin(x) - x * np.cos(x) )
   return f


def dW3d_sth(x):
   """Derivative of the FT of the top hat.
   Input and output are dimensionless.
   """
   f = 3. * (3. * x * np.cos(x) - 3. * np.sin(x) + (x**2) * np.sin(x)) / (x**4)
   return f


def W2d_cth(x):
   """FT of a 2d circular top hat window function.
   Input and output are dimensionless.
   """
   return 2.*special.jn(1, x) / x

def W1d_th(x):
   """FT of a 1d tophat
   normalized to unity at k=0 (ie real space integral is 1)
   Input and output are dimensionless.
   """
   return sinc(x/2.)
   
def Si(x):
   return special.sici(x)[0]

def Ci(x):
   return special.sici(x)[1]

def sinc(x):
   return special.sph_jn(0, x)[0][0]

def j0(x):
   """relevant for isotropic Fourier transform in 2d
   """
   return special.jn(0, x)


def i0(x):
   """Modified Bessel function of the first kind
   """
   return special.iv(0, x)


##################################################################################
# formatting numbers

def intExpForm(input):
   """
   clean scientific notation for file names
   removes trailing decimal point if not needed
   """
   a = '%e' % np.float(input)
   # mantissa: remove trailing zeros
   # then remove dot if no decimal digits
   mantissa = a.split('e')[0].rstrip('0').rstrip('.')
   # exponent: remove + sign if there, and leading zeros
   exponent = np.int(a.split('e')[1])
   exponent = np.str(exponent)
   if exponent=='0':
      return mantissa
   else:
      return mantissa + 'e' + exponent


def floatExpForm(input):
   """same as intExpForm, except always leaves the decimal point
   """
   a = '%e' % np.float(input)
   # mantissa: remove trailing zeros
   # then remove dot if no decimal digits
   mantissa = a.split('e')[0].rstrip('0')
   # exponent: remove + sign if there, and leading zeros
   exponent = np.int(a.split('e')[1])
   exponent = np.str(exponent)
   if exponent=='0':
      return mantissa
   else:
      return mantissa + 'e' + exponent


##################################################################################


def myHistogram(X, nBins=71, lim=(-1000., 1000.), sigma2Theory=None, path='./test.pdf', nameLatex=r'$x$ [km/s]', semilogx=False, semilogy=False, doGauss=False):
   """Generic histogram plotter.
   Flattens the input array X first thing.
   """
   # Flatten the array in case
   X = X.flatten()
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
      histGaussFit *= len(X)

   # Theory histogram
   if sigma2Theory is not None:
      av = 0.
      fPDF = lambda x: (2*np.pi*sigma2Theory)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2Theory))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-3)[0]
      histTheory = np.array(map(g, range(nBins-1)))
      histTheory *= len(X)

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
   if semilogy:
      ax.set_yscale('log', nonposy='clip')
   ax.set_ylim((0.5*np.min(histX[histX>0]), 2.*np.max(histX)))
   ax.set_xlabel(nameLatex)
   ax.set_ylabel(r'number of objects')
   fig.savefig(path, bbox_inches='tight')
   fig.clf()
#      plt.show()

