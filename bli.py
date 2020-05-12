
def intersectCatalog(self, newCat, save=False, vDiff=False, nProc=3):
   '''Take the intersection of the two catalogs.
   Keep the galaxy properties from the first catalog (self).
   If vDiff is True, use the difference of the two velocities, 
   for a null test.
   '''

   # find the intersection
   hasMatch = np.zeros(self.nObj, dtype=bool)

   def matchObj(iObj):
      if iObj%10000==0:
         print "matching object", iObj
      ra = self.RA[iObj]
      dec = self.DEC[iObj]
      z = self.Z[iObj]

      diff = (newCat.RA - ra)**2 / (1.e-3)**2   # accuracy of Mariana's RA (3.6 arcsec)
      diff += (newCat.DEC - dec)**2 / (1.e-4)**2   # accuracy of Mariana's DEC (0.36 arcsec)
      diff += (newCat.Z - z) **2 / (1.e-4)**2   # accuracy of Mariana's redshifts 
      diff = np.sqrt(diff)

      minDiff = np.min(diff)
      #print "min diff", minDiff
      if (minDiff<1.):
         IMatch = np.where(diff==minDiff)[0]
         if len(IMatch) > 1:
            print "Problem: got", len(IMatch), "matches"
         hasMatch[iObj] = True
         iMatch = IMatch[0]
         #print iObj, minDiff
      else:
         iMatch = -1
      return iMatch



   with sharedmem.MapReduce(np=nProc) as pool:
      #IMatch = np.array(pool.map(matchObj, range(self.nObj)))
      IMatch = np.array(pool.map(matchObj, range(500)))


   I0Match = np.where(IMatch<>-1.)[0]
   print "First catalog has", self.nObj, "objects"
   print "Second catalog has", newCat.nObj, "objects"
   print "Intersection has", len(I0Match), "objects"


   self.nObj = len(I0Match)
   #
   # sky coordinates and redshift
   self.RA = self.RA[I0Match]
   self.DEC = self.DEC[I0Match]
   self.Z = self.Z[I0Match]
   #
   # observed cartesian coordinates
   self.coordX = self.coordX[I0Match]
   self.coordY = self.coordY[I0Match]
   self.coordZ = self.coordZ[I0Match]
   #
   # displacement from difference,
   # not including the Kaiser displacement,
   # from differences of the observed and reconstructed fields
   self.dX = self.dX[I0Match]
   self.dY = self.dY[I0Match]
   self.dZ = self.dZ[I0Match]
   #
   # Kaiser-only displacement
   # originally from differences of the observed and reconstructed fields
   self.dXKaiser = self.dXKaiser[I0Match]
   self.dYKaiser = self.dYKaiser[I0Match]
   self.dZKaiser = self.dZKaiser[I0Match]
   #
   # velocity in cartesian coordinates
   self.vX = self.vX[I0Match]
   self.vY = self.vY[I0Match]
   self.vZ = self.vZ[I0Match]
   #
   # velocity in spherical coordinates,
   # from catalog of spherical displacements
   self.vR = self.vR[I0Match] - (vDiff==True) * newCat.vR[IMatch[I0Match]]
   self.vTheta = self.vTheta[I0Match] - (vDiff==True) * newCat.vTheta[IMatch[I0Match]]
   self.vPhi = self.vPhi[I0Match] - (vDiff==True) * newCat.vPhi[IMatch[I0Match]]
   #
   # Stellar masses
   self.Mstellar = self.Mstellar[I0Match]
   #
   # Halo mass
   self.hasM = self.hasM[I0Match]
   self.Mvir = self.Mvir[I0Match]
   #
   # Integrated optical depth [dimless]: int d^2theta n_e^2d sigma_T = (total nb of electrons) * sigma_T / (a chi)^2
   self.integratedTau = self.integratedTau[I0Match]
   #
   # Integrated kSZ signal [muK * sr]: int d^2theta n_e sigma_T (-v/c) Tcmb
   self.integratedKSZ = self.integratedKSZ[I0Match]
   #
   # Integrated Y signal [sr]: int d^2theta n_e sigma_T (kB Te / me c^2)
   # needs to be multiplied by Tcmb * f(nu) to get muK
   self.integratedY = self.integratedY[I0Match]


   # Write the full catalog to the output path, if needed
   if save:
      self.writeCatalog()
