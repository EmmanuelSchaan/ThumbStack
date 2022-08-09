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


def myHistogram(X, nBins=71, lim=None, S2Theory=[], path=None, plot=False, nameLatex=r'$x$', semilogx=False, semilogy=False, doGauss=False):
   """Generic histogram plotter.
   Flattens the input array X first thing.
   """
   # Flatten the array in case
   X = X.flatten()
   # value limits for the histogram
   if lim is None:
      lim = (np.min(X), np.max(X))
   # Bin edges
   if semilogx:
      Bins = np.logspace(np.log10(lim[0]), np.log10(lim[1]), nBins, 10.)
   else:
      Bins = np.linspace(lim[0], lim[1], nBins)
   binwidth = Bins[1:] - Bins[:-1]

   # Data histogram
   histX = np.histogram(X, Bins)[0]
   histX = histX.astype(float)

   # Plot
   fig = plt.figure(0)
   ax = fig.add_subplot(111)
   #
   # Histogram from data
   ax.bar(Bins[:-1], histX, binwidth, color='b', alpha=0.5, label=r'Data')
   #
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
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      #
      ax.axvline(1.*std, c='g', ls='--')
      ax.axvline(3.*std, c='g', ls='--')
      ax.axvline(5.*std, c='g', ls='--')
      ax.axvline(-1.*std, c='g', ls='--')
      ax.axvline(-3.*std, c='g', ls='--')
      ax.axvline(-5.*std, c='g', ls='--')
   #
   # Theory histogram
   for s2Theory in S2Theory:
      av = 0.
      fPDF = lambda x: (2*np.pi*s2Theory)**(-1./2.) * np.exp(-(x-av)**2 / (2*s2Theory))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-3)[0]
      histTheory = np.array(map(g, range(nBins-1)))
      histTheory *= len(X)
      ax.step(Bins[:-1], histTheory, color='r', lw=3, where='post')#, label=r'Theory')
   #
   ax.legend(loc=1)
   ax.set_xlim((lim[0], lim[1]))
   if semilogx:
      ax.set_xscale('log', nonposx='clip')
   if semilogy:
      ax.set_yscale('log', nonposy='clip')
   #
#   if doGauss:
#      yMin = min(histGaussFit.min(), 0.5*np.min(histX[histX>0]))
#      yMax = max(histGaussFit.max(), 2.*np.max(histX))
#   else:
#      yMin = 0.5*np.min(histX[histX>0])
#      yMax = 2.*np.max(histX)
   yMin = 0.5*np.min(histX[histX>0])
   yMax = 2.*np.max(histX)
   ax.set_ylim((yMin, yMax))
   ax.set_xlabel(nameLatex)
   ax.set_ylabel(r'number of objects')
   if path is not None:
      fig.savefig(path, bbox_inches='tight')
   if plot:
      plt.show()
   else:
      fig.clf()



def my2dHistogram(X, Y, nBins=(71, 71), limx=None, limy=None, limc=None, fTheory=[], path=None, plot=False, nameLatexX=r'$x$', nameLatexY=r'$y$', logx=False, logy=False, logColor=False, cmap=plt.cm.jet):
   """Generic 2d histogram plotter.
   """
   # limits for bins and colors
   if limx is None:
      limx = (np.min(X), np.max(X))
   if limy is None:
      limy = (np.min(Y), np.max(Y))
   if limc is None:
      limc = (None, None)
   # x-bin edges
   if logx:
      BinsX = np.logspace(np.log10(limx[0]), np.log10(limx[1]), nBins[0], 10.)
   else:
      BinsX = np.linspace(limx[0], limx[1], nBins[0])
   # y-bin edges
   if logy:
      BinsY = np.logspace(np.log10(limy[0]), np.log10(limy[1]), nBins[1], 10.)
   else:
      BinsY = np.linspace(limy[0], limy[1], nBins[1])

   # Plot
   fig = plt.figure(0)
   ax = fig.add_subplot(111)
   #
   # 2d histogram
   H, xEdges, yEdges = np.histogram2d(X, Y, bins=[BinsX, BinsY], range=[[limx[0], limx[1]],[limy[0], limy[1]]])
   H = H.T
   x, y = np.meshgrid(xEdges, yEdges)
   if logColor:
      im = ax.pcolormesh(x, y, H, cmap=cmap, linewidth=0, rasterized=True, norm=LogNorm(limc[0], limc[1]))
   else:
      im = ax.pcolormesh(x, y, H, cmap=cmap, linewidth=0, rasterized=True, vmin=limc[0], vmax=limc[1])
   #
   # theory curves, if any
   for f in fTheory:
      Ytheory = np.array(map(f, X))
      ax.plot(X, Ytheory, 'y')
   #
   plt.colorbar(im)
   if logx:
      ax.set_xscale('log', nonposx='clip')
   if logy:
      ax.set_yscale('log', nonposy='clip')
   ax.set_xlim((limx[0], limx[1]))
   ax.set_ylim((limy[0], limy[1]))
   ax.set_xlabel(nameLatexX)
   ax.set_ylabel(nameLatexY)
   if path is not None:
      fig.savefig(path, bbox_inches='tight')
   if plot:
      plt.show()
   else:
      fig.clf()


def splitBins(x, nBins):
   '''Give the nBins+1 edges of the nBins bins,
   such that there is an equal number of elements of x
   in each bin.
   '''
   nBins = np.int(nBins)
   # sort the x values
   xSorted = np.sort(x)
   # corresponding probability
   proba = 1. * np.arange(len(x)) / (len(x) - 1)
   # interpolate the CDF
   # scipy is smart: even if twice the same x appears,
   # it will know to give the interpolating function the right jump
   cdf = interp1d(xSorted, proba, kind='linear', bounds_error=False, fill_value=(0., 1.))

   # fill an array with the z-bounds of the bins
   binEdges = np.zeros(nBins+1)
   binEdges[0] = np.min(x)
   binEdges[-1] = np.max(x)
   for iBin in range(nBins-1):
      # force equal number of objects per bins
      f = lambda xMax: cdf(xMax) - (iBin+1.)/nBins
      binEdges[iBin+1] = optimize.brentq(f , binEdges[0], binEdges[-1])
   return binEdges


def computeSnr(d, theory, cov, dof=None):
   """Compute null rejection, SNR (=detection significance)
   """
   if dof is None:
      dof = len(d)

   # Compute chi^2_null
   chi2Null = d.dot( np.linalg.inv(cov).dot(d) )
   # goodness of fit for null hypothesis
   print("number of dof:"+str(dof))
   print("null chi2Null="+str(chi2Null))
   print("SNR = null sqrt(chi2Null)="+str(np.sqrt(chi2Null)))
#   pteNull = 1.- stats.chi2.cdf(chi2Null, dof)
#   print("null pte="+str(pteNull))
   pteNull = stats.distributions.chi2.sf(chi2Null, dof)
   print("null pte="+str(pteNull))
   # pte as a function of sigma, for a Gaussian random variable
   fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.)) - pteNull
   sigmaNull = optimize.brentq(fsigmaToPTE , 0., 1.e5)
   print("null pte significance="+str(sigmaNull)+"sigmas")

   # find the best fit amplitude for the theory curve
   def fdchi2(p):
      a = p[0]
      result = (d-a*theory).dot( np.linalg.inv(cov).dot(d-a*theory) )
      result -= chi2Null
      return result
   # Minimize the chi squared
   p0 = 1.
   res = optimize.minimize(fdchi2, p0)
   abest = res.x[0]
   #sbest= res.x[1]
   print("best-fit amplitude="+str(abest))
   print("number of dof:"+str(dof - 1))

   # goodness of fit for best fit
   chi2Best = fdchi2([abest])+chi2Null
   print("best-fit chi2="+str(chi2Best))
   pteBest = 1.- stats.chi2.cdf(chi2Best, dof-1.)
   print("best-fit pte="+str(pteBest))
   # pte as a function of sigma, for a Gaussian random variable
   fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.)) - pteBest
   sigma = optimize.brentq(fsigmaToPTE , 0., 1.e3)
   print("best-fit pte significance="+str(sigma)+"sigmas")

   # favour of best fit over null
   print("best-fit sqrt(delta chi2)="+str(np.sqrt(abs(fdchi2([abest]))))+"sigmas")
   fsigmaToPTE = lambda sigma: special.erfc(sigma/np.sqrt(2.))
   pte = fsigmaToPTE( np.sqrt(abs(fdchi2([abest]))) )
   print("pte (if Gaussian)="+str(pte))

