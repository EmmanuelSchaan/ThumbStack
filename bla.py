   def plotHistogram(self, X, nBins=71, lim=(-1000., 1000.), sigma2Theory=None, name='x', nameLatex=r'$x$ [m/s]'):
      # histogram parameters
      binwidth = (lim[0]-lim[1])/(nBins-1)
      Bins = np.linspace(lim[0], lim[1], nBins)

      # Data histogram
      histX = np.histogram(X, Bins)[0]
      histX = histX.astype(float)

      # histogram for a Gaussian with the variance from the data
      mean = np.mean(X)
      std = np.std(X)
      sigma2 = std**2
      av = mean
      fPDF = lambda x: (2*np.pi*sigma2)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2))
      g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
      histGaussFit = np.array(map( g, np.linspace(0, nBins-2, nBins-1) ))
      histGaussFit *= self.nObj

      # Theory histogram
      if sigma2Theory is not None:
         av = 0.
         fPDF = lambda x: (2*np.pi*sigma2Theory)**(-1./2.) * np.exp(-(x-av)**2 / (2*sigma2Theory))
         g = lambda i: integrate.quad(fPDF, Bins[i], Bins[i+1], epsabs=0, epsrel=1.e-5)[0]
         histTheory = np.array(map( g, np.linspace(0, nBins-2, nBins-1) ))
         histTheory *= self.nObj

      # Plot
      fig = plt.figure(0)
      ax = fig.add_subplot(111)
      #
      ax.bar(Bins[:-1], histX, binwidth, color='b', alpha=0.5, label=r'Data')
      ax.step(Bins[:-1], histGaussFit, color='g', lw=3, where='post', label=r'Gaussian')
      if sigma2Theory is not None:
         ax.step(Bins[:-1], histTheory, color='r', lw=3, where='post', label=r'Theory')
      #
      infos = str(len(X))+" objects"+"\n"+"mean="+str(round(mean,1))+"\n"+r"$\sigma$="+str(round(std,1))
      ax.text(-900.,  10500., infos, fontsize=16,
              horizontalalignment='left',
              verticalalignment='center')
      #
      ax.legend(loc=1)
      ax.set_xlim((lim[0], lim[1]))
      ax.set_xlabel(nameLatex)
      ax.set_ylabel(r'number of objects')
      fig.savefig(self.pathFig+"/hist_"+name+".pdf")
      fig.clf()
