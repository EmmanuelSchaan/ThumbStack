from headers import *


def my2dHistogram(X, Y, nBins=(71, 71), limx=None, limy=None, limc=None, fTheory=[], path='./test.pdf', nameLatexX=r'$x$', nameLatexY=r'$y$', logx=False, logy=False, logColor=False, cmap=plt.cm.jet):
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
   H, xEdges, yEdges = np.histogram2d(X, Y, bins=[BinsX, BinsY], range=[[limx[0], limx[1]],[limy[0], limy[1]]], density=False)
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
   fig.savefig(path, bbox_inches='tight')
   fig.clf()
#   plt.show()







X = np.logspace(np.log10(1.), np.log10(1.e4), 1001, 10)


Y = 1.e3 + X


#X += np.random.normal(loc=0., scale=0.01*X, size=len(X))
#Y += np.random.normal(loc=0., scale=0.01*Y, size=len(Y))


f = lambda x: 1.e3 + x

my2dHistogram(X, Y, nBins=(71, 71), limx=None, limy=None, limc=None, fTheory=[f], path='./test.pdf', nameLatexX=r'$x$', nameLatexY=r'$y$', logx=False, logy=False, logColor=True, cmap=plt.cm.Greys)

