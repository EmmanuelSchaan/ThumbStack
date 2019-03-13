           # perform the filtering
            filtMap[iRAp], filtMask[iRAp], filtNoiseStdDev[iRAp], diskArea[iRAp] = self.diskRingFilter(opos, stampMap, stampMask, stampHit, r0, r1, test=test)

      return filtMap, filtMask, filtNoiseStdDev, diskArea



   def saveFiltering(self, nProc=1):
      
      # initialize arrays
      self.filtMap = np.zeros((self.Catalog.nObj, self.nRAp))
      self.filtMask = np.zeros((self.Catalog.nObj, self.nRAp))
      self.filtNoiseStdDev = np.zeros((self.Catalog.nObj, self.nRAp))
      self.diskArea = np.zeros((self.Catalog.nObj, self.nRAp))

      # loop over all objects in catalog
#      result = np.array(map(self.analyzeObject, range(self.Catalog.nObj)))
      tStart = time()
      with sharedmem.MapReduce(np=nProc) as pool:
         result = np.array(pool.map(self.analyzeObject, range(self.Catalog.nObj)))
      tStop = time()
      print "took", (tStop-tStart)/60., "min"
      
      # unpack and save to file
      self.filtMap = result[:,0,:].copy()
      np.savetxt(self.pathOut+"/filtmap.txt", self.filtMap)
      #
      self.filtMask = result[:,1,:].copy()
      np.savetxt(self.pathOut+"/filtmask.txt", self.filtMask)
      #
      self.filtNoiseStdDev = result[:,2,:].copy()
      np.savetxt(self.pathOut+"/filtnoisestddev.txt", self.filtNoiseStdDev)
      #
      self.diskArea = result[:,3,:].copy()
      np.savetxt(self.pathOut+"/diskarea.txt", self.diskArea)
