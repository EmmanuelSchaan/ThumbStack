   def plot(self, data=None, save=False, path=None, cmap='viridis'):
      if data is None:
         data = self.data.copy()
      sigma = np.std(data.flatten())
      vmin = np.min(data.flatten())
      vmax = np.max(data.flatten())

      fig=plt.figure(0)
      ax=fig.add_subplot(111)
      #
      # pcolor wants x and y to be edges of cell,
      # ie one more element, and offset by half a cell
      x = self.dX * (np.arange(self.nX+1) - 0.5)
      y = self.dY * (np.arange(self.nY+1) - 0.5)
      x,y = np.meshgrid(x, y, indexing='ij')
      #
      cp=ax.pcolormesh(x*180./np.pi, y*180./np.pi, data, linewidth=0, rasterized=True)
      #
      # choose color map: jet, summer, winter, Reds, gist_gray, YlOrRd, bwr, seismic
      cp.set_cmap(cmap)
      #cp.set_clim(0.,255.)
      #cp.set_clim(-3.*sigma, 3.*sigma)
      cp.set_clim(vmin, vmax)
      fig.colorbar(cp)
      #
      plt.axis('scaled')
      ax.set_xlim(np.min(x)*180./np.pi, np.max(x)*180./np.pi)
      ax.set_ylim(np.min(y)*180./np.pi, np.max(y)*180./np.pi)
      ax.set_xlabel('$x$ [deg]')
      ax.set_ylabel('$y$ [deg]')
      #
      if save==True:
         if path is None:
            path = "./figures/lens_simulator/"+self.name+".pdf"
         print "saving plot to "+path
         fig.savefig(path, bbox_inches='tight')
         fig.clf()
      else:
         plt.show()
