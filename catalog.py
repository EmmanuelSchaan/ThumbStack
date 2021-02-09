from headers import *

##################################################################################
##################################################################################

class Catalog(object):

   def __init__(self, U, MassConversion, name="test", nameLong=None, pathInCatalog="", rV=1.,  save=False, nObj=None):
      '''nObj: used to keep the first nObj objects of the catalog, useful for quick debugging
      '''

      self.U = U
      self.MassConversion = MassConversion
      self.name = name
      self.rV = rV   # velocity real-space correlation coefficient
      if nameLong is None:
         self.nameLong = self.name
      else:
         self.nameLong = nameLong
      self.pathInCatalog = pathInCatalog
      
      # Output path
      self.pathOut = "./output/catalog/"+self.name
      if not os.path.exists(self.pathOut):
         os.makedirs(self.pathOut)
      # catalog path
      self.pathOutCatalog = self.pathOut + "/catalog.txt"
      # path for vtk file (to visualize with VisIt)
      self.pathOutVtk = self.pathOut + "/catalog.vtk"
      
      # Figures path
      self.pathFig = "./figures/catalog/"+self.name
      if not os.path.exists(self.pathFig):
         os.makedirs(self.pathFig)
      
      
      if save:
         self.readInputCatalog()
         self.addHaloMass()
         self.addIntegratedTau()
         self.addIntegratedKSZ()
         self.addIntegratedY()
         self.writeCatalog()

      self.loadCatalog(nObj=nObj)
   

   ##################################################################################
   ##################################################################################

   def copy(self, name="test", nameLong=None):
      """Copy a catalog class, with the option of changing the name.
      """
      # First copy the output catalog
      # new catalog path
      newPathOut = "./output/catalog/"+name
      if not os.path.exists(newPathOut):
         os.makedirs(newPathOut)
      newPathOutCatalog = newPathOut + "/catalog.txt"
      # copy the output catalog
      copyfile(self.pathOutCatalog, newPathOutCatalog)
      
      # Then copy the catalog properties
      newCat = Catalog(self.U, self.MassConversion, name=name, nameLong=nameLong, pathInCatalog=self.pathInCatalog, save=False)
      return newCat


   def extractCatalog(self, I, name="test", nameLong=None):
      """create and return a new catalog object,
      keeping only the objects with indices in I
      """
      # new catalog path
      newPathOut = "./output/catalog/"+name
      if not os.path.exists(newPathOut):
         os.makedirs(newPathOut)
      newPathOutCatalog = newPathOut + "/catalog.txt"
      # read the current catalog
      data = np.genfromtxt(self.pathOutCatalog)
      # keep only objects with indices in I
      data = data[I,:]
      # save it to new catalog path
      np.savetxt(newPathOutCatalog, data)
      # Then copy the catalog properties
      newCat = Catalog(self.U, self.MassConversion, name=name, nameLong=nameLong, pathInCatalog=self.pathInCatalog, save=False)
      return newCat


   ##################################################################################
   ##################################################################################

   def readInputCatalog(self):
      print "- read input catalog from "+self.pathInCatalog
      data = np.genfromtxt(self.pathInCatalog)
      self.nObj = len(data[:,0])
      #
      # sky coordinates and redshift
      self.RA = data[:,0] # [deg]
      self.DEC = data[:,1]   # [deg]
      self.Z = data[:,2]
      #
      # observed cartesian coordinates
      self.coordX = data[:,3]   # [Mpc/h]
      self.coordY = data[:,4]   # [Mpc/h]
      self.coordZ = data[:,5]   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      self.dX = data[:,6]   # [Mpc/h]
      self.dY = data[:,7]   # [Mpc/h]
      self.dZ = data[:,8]   # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      self.dXKaiser = data[:,9]   # [Mpc/h] from cartesian catalog difference
      self.dYKaiser = data[:,10]   # [Mpc/h]
      self.dZKaiser = data[:,11]   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      self.vX = data[:,12]   #[km/s]
      self.vY = data[:,13]   #[km/s]
      self.vZ = data[:,14]   #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      self.vR = data[:,15]  # [km/s]   from spherical catalog, >0 away from us
      self.vTheta = data[:,16]   # [km/s]
      self.vPhi = data[:,17]  # [km/s]
      #
      # Stellar masses
      self.Mstellar = data[:,18]   # [M_sun], from Maraston et al


   ##################################################################################

   def addHaloMass(self):
      """Generate halo masses in M_sun,
      from stellar masses in M_sun.
      """
      print "- add halo masses"
      # flag: 1 if object has mass
      self.hasM = np.zeros(self.nObj)
      
      self.Mvir = np.zeros(self.nObj)
      for iObj in range(self.nObj):
         mStellar = self.Mstellar[iObj]
         if (mStellar>1.e3) and not np.isnan(mStellar):
            self.hasM[iObj] = True
            self.Mvir[iObj] = self.MassConversion.fmStarTomVir(mStellar)

      # for object without a mass, use the mean mass from the others
      if np.sum(self.hasM)>0:
         meanMstellar = np.mean(self.Mstellar[self.hasM.astype('bool')])
         self.Mstellar[~self.hasM.astype('bool')] = meanMstellar
         #
         meanMvir = np.mean(self.Mvir[self.hasM.astype('bool')])
         self.Mvir[~self.hasM.astype('bool')] = meanMvir
      # if no object has a mass, make a random guess, rather than keeping 0
      else:
         self.Mstellar = 2.6e11   # random guess
         self.Mvir = self.MassConversion.fmStarTomVir(2.6e11)



   def addIntegratedTau(self):
      """integrated optical depth to Thompson scattering: int d^2theta n_e^2d sigma_T
      = (total nb of electrons) * sigma_T / (a chi)^2
      [sr]
      """
      print "- add integrated tau"
      # convert from total mass to baryon mass
      # assuming cosmological abundance of baryons, and all baryons are in gas
      self.integratedTau = self.Mvir * self.U.bg.Omega0_b/self.U.bg.Omega0_m
      
      # convert from total baryon mass to electron total number
      me = 9.10938291e-31  # electron mass (kg)
      mH = 1.67262178e-27  # proton mass (kg)
      mHe = 4.*mH # helium nucleus mass (kg)
      xH = 0.76   # primordial hydrogen fraction by mass
      nH_ne = 2.*xH/(xH+1.)
      nHe_ne = (1.-xH)/(2.*(1.+xH))
      msun = 1.989e30   # solar mass (kg)
      factor = (me + nH_ne*mH + nHe_ne*mHe) * (1./msun)   # total mass per electron in (Msun)
      self.integratedTau /= factor
      
      # multiply by Thomson cross section (physical)
      mpc = 3.08567758e16*1.e6   # 1Mpc in m
      sigma_T = 6.6524e-29 # Thomson cross section in m^2
      self.integratedTau *= sigma_T / (mpc / self.U.bg.h)**2
      
      # divide by (a chi)^2
      self.integratedTau /= (self.U.bg.comoving_distance(self.Z) / (1.+self.Z))**2

   
   def addIntegratedKSZ(self):
      """Integrated kSZ signal: int d^2theta n_e sigma_T (-v/c) Tcmb
      in [muK * sr]
      """
      print "- add integrated kSZ"
      self.integratedKSZ = - self.integratedTau * (self.vR/3.e5) * 2.726e6

   
   def addIntegratedY(self, nu=150.e9):
      """Integrated tSZ signal: int d^2theta n_e sigma_T (k_B T_e / m_e c^2)
      in [sr].
      To get dT in muK*sr, multiply by Tcmb * f(nu).
      Simple power-law fit to Greco et al 2014, fig4.
      """
      print "- add integrated y"
      
      # in arcmin^2
      yCcyltilda = (self.Mstellar/1.e11)**3.2 * 1.e-6

      # in arcmin^2
      yCcyl = yCcyltilda * (self.U.hubble(self.Z) / self.U.hubble(0.))**(2./3.)
      yCcyl /= (self.U.bg.comoving_distance(self.Z) / (500.*self.U.bg.h))**2
      # in sr
      yCcyl *= (np.pi/180./60.)**2
      
      self.integratedY = yCcyl


   ##################################################################################

   def writeCatalog(self):
      print "- write full catalog to "+self.pathOutCatalog
      data = np.zeros((self.nObj,24))
      #
      # sky coordinates and redshift
      data[:,0] = self.RA # [deg]
      data[:,1] = self.DEC   # [deg]
      data[:,2] = self.Z
      #
      # observed cartesian coordinates
      data[:,3] = self.coordX   # [Mpc/h]
      data[:,4] = self.coordY   # [Mpc/h]
      data[:,5] = self.coordZ   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      data[:,6] = self.dX   # [Mpc/h]
      data[:,7] = self.dY   # [Mpc/h]
      data[:,8] = self.dZ   # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      data[:,9] = self.dXKaiser   # [Mpc/h] from cartesian catalog difference
      data[:,10] = self.dYKaiser   # [Mpc/h]
      data[:,11] = self.dZKaiser   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      data[:,12] = self.vX   #[km/s]
      data[:,13] = self.vY   #[km/s]
      data[:,14] = self.vZ   #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      data[:,15] = self.vR  # [km/s]   from spherical catalog, >0 away from us
      data[:,16] = self.vTheta   # [km/s]
      data[:,17] = self.vPhi  # [km/s]
      #
      # Stellar mass
      data[:,18] = self.Mstellar   # [M_sun], from Maraston et al
      #
      # Halo mass
      data[:,19] = self.hasM  # flag=1 if mass is known
      data[:,20] = self.Mvir   # [M_sun]
      #
      # Integrated optical depth [sr]: int d^2theta n_e^2d sigma_T = (total nb of electrons) * sigma_T / (a chi)^2
      data[:,21] = self.integratedTau   # [sr]
      #
      # Integrated kSZ signal [muK * sr]: int d^2theta n_e^2d sigma_T (-v/c) Tcmb
      data[:, 22] = self.integratedKSZ # [muK * sr]
      #
      # Integrated Y signal [sr]: int d^2theta n_e^2d sigma_T (kB Te / me c^2)
      # needs to be multiplied by Tcmb * f(nu) to get muK
      data[:, 23] = self.integratedY # [sr]
      #
      np.savetxt(self.pathOutCatalog, data)


   def loadCatalog(self, nObj=None):
      print "- load full catalog from "+self.pathOutCatalog
      data = np.genfromtxt(self.pathOutCatalog)
      self.nObj = len(data[:nObj,0])
      #
      # sky coordinates and redshift
      self.RA = data[:nObj,0] # [deg]
      self.DEC = data[:nObj,1]   # [deg]
      self.Z = data[:nObj,2]
      #
      # observed cartesian coordinates
      self.coordX = data[:nObj,3]   # [Mpc/h]
      self.coordY = data[:nObj,4]   # [Mpc/h]
      self.coordZ = data[:nObj,5]   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      self.dX = data[:nObj,6]   # [Mpc/h]
      self.dY = data[:nObj,7]   # [Mpc/h]
      self.dZ = data[:nObj,8]   # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      self.dXKaiser = data[:nObj,9]   # [Mpc/h] from cartesian catalog difference
      self.dYKaiser = data[:nObj,10]   # [Mpc/h]
      self.dZKaiser = data[:nObj,11]   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      self.vX = data[:nObj,12]   #[km/s]
      self.vY = data[:nObj,13]   #[km/s]
      self.vZ = data[:nObj,14]   #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      self.vR = data[:nObj,15]  # [km/s]   from spherical catalog, >0 away from us
      self.vTheta = data[:nObj,16]   # [km/s]
      self.vPhi = data[:nObj,17]  # [km/s]
      #
      # Stellar masses
      self.Mstellar = data[:nObj,18]   # [M_sun], from Maraston et al
      #
      # Halo mass
      self.hasM = data[:nObj,19]
      self.Mvir = data[:nObj,20]  # [M_sun]
      #
      # Integrated optical depth [sr]: int d^2theta n_e^2d sigma_T = (total nb of electrons) * sigma_T / (a chi)^2
      self.integratedTau = data[:nObj,21]   # [sr]
      #
      # Integrated kSZ signal [muK * sr]: int d^2theta n_e sigma_T (-v/c) Tcmb
      self.integratedKSZ = data[:nObj, 22] # [muK * sr]
      #
      # Integrated Y signal [sr]: int d^2theta n_e sigma_T (kB Te / me c^2)
      # needs to be multiplied by Tcmb * f(nu) to get muK
      self.integratedY = data[:nObj, 23] # [sr]


   ##################################################################################
   ##################################################################################


   def writeVtk(self):
      print "- save catalog to vtk file:"
      print self.pathOutVtk
      # create vtk file
      f = open(self.pathOutVtk,'w')
      f.write('# vtk DataFile Version 3.0\n')
      f.write('Catalog data for visualization\n')
      f.write('ASCII\n')
      f.write('\n')

      # Cartesian coord
      f.write('DATASET POLYDATA\n')
      f.write('POINTS '+str(self.nObj)+' DOUBLE\n')
      for iobj in range(self.nObj):
         f.write(format(self.coordX[iobj], '.10f')+' '+format(self.coordY[iobj], '.10f')+' '+format(self.coordZ[iobj], '.10f')+'\n')
      f.write('\n')
      f.write('POINT_DATA '+str(self.nObj)+'\n')

      # Velocity: cartesian
      f.write('VECTORS vel DOUBLE\n')
      for iobj in range(self.nObj):
         f.write(format(self.vX[iobj], '.10e')+' '+format(self.vY[iobj], '.10e')+' '+format(self.vZ[iobj], '.10e')+'\n')
      f.write('\n')

      # Displacement: cartesian
      f.write('VECTORS disp DOUBLE\n')
      for iobj in range(self.nObj):
         f.write(format(self.dX[iobj], '.10e')+' '+format(self.dY[iobj], '.10e')+' '+format(self.dZ[iobj], '.10e')+'\n')
      f.write('\n')

      # Kaiser displacement:cartesian
      f.write('VECTORS dispKaiser DOUBLE\n')
      for iobj in range(self.nObj):
         f.write(format(self.dXKaiser[iobj], '.10e')+' '+format(self.dYKaiser[iobj], '.10e')+' '+format(self.dZKaiser[iobj], '.10e')+'\n')
      f.write('\n')

      # Velocity: spherical
      f.write('VECTORS velSph DOUBLE\n')
      for iobj in range(self.nObj):
         f.write(format(self.vR[iobj], '.10e')+' '+format(self.vTheta[iobj], '.10e')+' '+format(self.vPhi[iobj], '.10e')+'\n')
      f.write('\n')

      # RA
      f.write('SCALARS RA DOUBLE\n')
      f.write('LOOKUP_TABLE default\n')
      for iobj in range(self.nObj):
         f.write(format(self.RA[iobj], '.10e')+'\n')
      f.write('\n')

      # DEC
      f.write('SCALARS DEC DOUBLE\n')
      f.write('LOOKUP_TABLE default\n')
      for iobj in range(self.nObj):
         f.write(format(self.DEC[iobj], '.10e')+'\n')
      f.write('\n')

      # Redshift
      f.write('SCALARS Redshift DOUBLE\n')
      f.write('LOOKUP_TABLE default\n')
      for iobj in range(self.nObj):
         f.write(format(self.Z[iobj], '.10e')+'\n')
      f.write('\n')

      # Stellar mass
      f.write('SCALARS mStellar DOUBLE\n')
      f.write('LOOKUP_TABLE default\n')
      for iobj in range(self.nObj):
         f.write(format(self.Mstellar[iobj], '.10e')+'\n')
      f.write('\n')

      ## output the virial mass, in Msun/h
      #f.write('SCALARS mVir DOUBLE\n')
      #f.write('LOOKUP_TABLE default\n')
      #for iobj in range(self.nObj):
      #   f.write(format(self.Mvir[iobj], '.10e')+'\n')
      #f.write('\n')

      ## output the expected kSZ, in muK
      #f.write('SCALARS expKSZ DOUBLE\n')
      #f.write('LOOKUP_TABLE default\n')
      #for iobj in range(self.nObj):
      #   f.write(format(self.ExpectedKSZ[iobj], '.10e')+'\n')
      #f.write('\n')

      ## output the expected tSZ, in muK
      #f.write('SCALARS expTSZ DOUBLE\n')
      #f.write('LOOKUP_TABLE default\n')
      #for iobj in range(self.nObj):
      #   f.write(format(self.ExpectedTSZ[iobj], '.10e')+'\n')
      #f.write('\n')

      f.close()




   ##################################################################################
   ##################################################################################
   
   def addCatalog(self, newCat, save=False):
      """Combines the current catalog with a new catalog newCat.
      """
      # number of objects
      self.nObj += newCat.nObj
      #
      # sky coordinates and redshift
      self.RA = np.concatenate((self.RA, newCat.RA)) # [deg]
      self.DEC = np.concatenate((self.DEC, newCat.DEC))   # [deg]
      self.Z = np.concatenate((self.Z, newCat.Z))
      #
      # observed cartesian coordinates
      self.coordX = np.concatenate((self.coordX, newCat.coordX))   # [Mpc/h]
      self.coordY = np.concatenate((self.coordY, newCat.coordY))   # [Mpc/h]
      self.coordZ = np.concatenate((self.coordZ, newCat.coordZ))   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      self.dX = np.concatenate((self.dX, newCat.dX))  # [Mpc/h]
      self.dY = np.concatenate((self.dY, newCat.dY))   # [Mpc/h]
      self.dZ = np.concatenate((self.dZ, newCat.dZ))  # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      self.dXKaiser = np.concatenate((self.dXKaiser, newCat.dXKaiser))  # [Mpc/h] from cartesian catalog difference
      self.dYKaiser = np.concatenate((self.dYKaiser, newCat.dYKaiser))   # [Mpc/h]
      self.dZKaiser = np.concatenate((self.dZKaiser, newCat.dZKaiser))   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      self.vX = np.concatenate((self.vX, newCat.vX))   #[km/s]
      self.vY = np.concatenate((self.vY, newCat.vY))   #[km/s]
      self.vZ = np.concatenate((self.vZ, newCat.vZ))  #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      self.vR = np.concatenate((self.vR, newCat.vR))  # [km/s]   from spherical catalo
      self.vTheta = np.concatenate((self.vTheta, newCat.vTheta))   # [km/s]
      self.vPhi = np.concatenate((self.vPhi, newCat.vPhi))  # [km/s]
      #
      # Stellar masses
      self.Mstellar = np.concatenate((self.Mstellar, newCat.Mstellar))   # [M_sun], from Maraston et al
      #
      # Halo mass
      self.hasM = np.concatenate((self.hasM, newCat.hasM))
      self.Mvir = np.concatenate((self.Mvir, newCat.Mvir))  # [M_sun]
      #
      # Integrated optical depth [dimless]: int d^2theta n_e^2d sigma_T = (total nb of electrons) * sigma_T / (a chi)^2
      self.integratedTau = np.concatenate((self.integratedTau, newCat.integratedTau))   # [dimless]
      #
      # Integrated kSZ signal [muK * sr]: int d^2theta n_e sigma_T v/c Tcmb
      self.integratedKSZ = np.concatenate((self.integratedKSZ, newCat.integratedKSZ)) # [muK * sr]
      #
      # Integrated Y signal [sr]: int d^2theta n_e sigma_T (kB Te / me c^2)
      # needs to be multiplied by Tcmb * f(nu) to get muK
      self.integratedY = np.concatenate((self.integratedY, newCat.integratedY)) # [sr]

      # Write the full catalog to the output path, if needed
      if save:
         self.writeCatalog()



   ##################################################################################
   ##################################################################################


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
         IMatch = np.array(pool.map(matchObj, range(self.nObj)))
         #IMatch = np.array(pool.map(matchObj, range(500)))


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


   ##################################################################################
   ##################################################################################
   
   def plotFootprint(self, hMap=None):
      """Overlay a scatter plot of the catalog positions on top of a healpix map,
      here the AdvACT hit count map.
      """

      if hMap is None:
         #pathMap = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2020_02_28/" + "act_planck_s08_s18_cmb_f150_daynight_map.fits"
         pathMask = "./output/cmb_map/pact20200228_r2/" + "mask_full_foot_gal60_ps.fits"
         pathHit = "/global/cscratch1/sd/eschaan/project_ucsc/data/planck_act_coadd_2020_02_28_r2/" + "act_planck_s08_s18_cmb_f150_daynight_ivar.fits"
         #
         #hMap = enmap.read_map(pathMap)[0] * enmap.read_map(pathMask)
         hMap = enmap.read_map(pathHit)[0] 
         hMap *= enmap.read_map(pathMask)[0]
         hMap = np.log(np.abs(hMap)+1.e-5)
         hMap = enmap.to_healpix(hMap)

      fig=plt.figure(0)
      #
      # hit count map for AdvACT
      hp.mollview(hMap, fig=0, title="", coord=None, cbar=False, unit='')
      #
      # scatter plot of the catalog
      hp.projscatter(self.RA, self.DEC, alpha=0.01, lonlat=True, marker='.', c='r', rasterized=True)
      #
      fig.savefig(self.pathFig+"/footprint_"+self.name+".pdf", dpi=1200, bbox_inches='tight')
      fig.clf()


   ##################################################################################
   ##################################################################################
   
   def addCatalog(self, newCat, save=False):
      """Combines the current catalog with a new catalog newCat.
      """
      # number of objects
      self.nObj += newCat.nObj
      #
      # sky coordinates and redshift
      self.RA = np.concatenate((self.RA, newCat.RA)) # [deg]
      self.DEC = np.concatenate((self.DEC, newCat.DEC))   # [deg]
      self.Z = np.concatenate((self.Z, newCat.Z))
      #
      # observed cartesian coordinates
      self.coordX = np.concatenate((self.coordX, newCat.coordX))   # [Mpc/h]
      self.coordY = np.concatenate((self.coordY, newCat.coordY))   # [Mpc/h]
      self.coordZ = np.concatenate((self.coordZ, newCat.coordZ))   # [Mpc/h]
      #
      # displacement from difference,
      # not including the Kaiser displacement,
      # from differences of the observed and reconstructed fields
      self.dX = np.concatenate((self.dX, newCat.dX))  # [Mpc/h]
      self.dY = np.concatenate((self.dY, newCat.dY))   # [Mpc/h]
      self.dZ = np.concatenate((self.dZ, newCat.dZ))  # [Mpc/h]
      #
      # Kaiser-only displacement
      # originally from differences of the observed and reconstructed fields
      self.dXKaiser = np.concatenate((self.dXKaiser, newCat.dXKaiser))  # [Mpc/h] from cartesian catalog difference
      self.dYKaiser = np.concatenate((self.dYKaiser, newCat.dYKaiser))   # [Mpc/h]
      self.dZKaiser = np.concatenate((self.dZKaiser, newCat.dZKaiser))   # [Mpc/h]
      #
      # velocity in cartesian coordinates
      self.vX = np.concatenate((self.vX, newCat.vX))   #[km/s]
      self.vY = np.concatenate((self.vY, newCat.vY))   #[km/s]
      self.vZ = np.concatenate((self.vZ, newCat.vZ))  #[km/s]
      #
      # velocity in spherical coordinates,
      # from catalog of spherical displacements
      self.vR = np.concatenate((self.vR, newCat.vR))  # [km/s]   from spherical catalo
      self.vTheta = np.concatenate((self.vTheta, newCat.vTheta))   # [km/s]
      self.vPhi = np.concatenate((self.vPhi, newCat.vPhi))  # [km/s]
      #
      # Stellar masses
      self.Mstellar = np.concatenate((self.Mstellar, newCat.Mstellar))   # [M_sun], from Maraston et al
      #
      # Halo mass
      self.hasM = np.concatenate((self.hasM, newCat.hasM))
      self.Mvir = np.concatenate((self.Mvir, newCat.Mvir))  # [M_sun]
      #
      # Integrated optical depth [dimless]: int d^2theta n_e^2d sigma_T = (total nb of electrons) * sigma_T / (a chi)^2
      self.integratedTau = np.concatenate((self.integratedTau, newCat.integratedTau))   # [dimless]
      #
      # Integrated kSZ signal [muK * sr]: int d^2theta n_e sigma_T v/c Tcmb
      self.integratedKSZ = np.concatenate((self.integratedKSZ, newCat.integratedKSZ)) # [muK * sr]
      #
      # Integrated Y signal [sr]: int d^2theta n_e sigma_T (kB Te / me c^2)
      # needs to be multiplied by Tcmb * f(nu) to get muK
      self.integratedY = np.concatenate((self.integratedY, newCat.integratedY)) # [sr]

      # Write the full catalog to the output path, if needed
      if save:
         self.writeCatalog()



   ##################################################################################
   ##################################################################################
   
#   def plotFootprint(self):
#      """Overlay a scatter plot of the catalog positions on top of a healpix map,
#      here the AdvACT hit count map.
#      """
#      fig=plt.figure(0)
#      #
#      # hit count map for AdvACT
#      path = "/global/cscratch1/sd/eschaan/project_ksz_act_planck/data/planck_act_coadd_2018_08_10/healpix_f150_daynight_all_div_mono.fits"
#      hHitMap = hp.read_map(path)
#      hp.mollview(np.log(np.abs(hHitMap)+1.e-5), fig=0, title="", coord=None, cbar=False, unit='')
#      #
#      # scatter plot of the catalog
#      hp.projscatter(self.RA, self.DEC, alpha=0.01, lonlat=True, marker='.', c='r', rasterized=True)
#      #
#      fig.savefig(self.pathFig+"/footprint_"+self.name+".pdf", dpi=1200)
#      fig.clf()


   ##################################################################################
   ##################################################################################
   
   def printProperties(self):
      print "Catalog: "+self.nameLong
      print "Number of objects = "+str(self.nObj)
      print "with mass: "+str(np.sum(self.hasM))+", i.e. fraction "+str(np.sum(self.hasM)/self.nObj)
      print "Z: mean = "+str(np.mean(self.Z))+", median = "+str(np.median(self.Z))
      m = self.Mstellar[self.hasM==1]
      print "M_star [M_sun]: mean = "+str(np.mean(m))+", median = "+str(np.median(m))
      m = self.Mvir[self.hasM==1]
      print "M_vir [M_sun]: mean = "+str(np.mean(m))+", median = "+str(np.median(m))
   


   ##################################################################################
   ##################################################################################

   def compareV1dRms(self):
      """expected std dev of velocities for LCDM:
      comparison between using median z or z-distribution
      """
      # interpolate RMS 1d velocity for speed
      f = lambda z: self.U.v3dRms(0., z, W3d_sth) / np.sqrt(3.)
      Z = np.linspace(0., 1., 201)
      V1dRms = np.array(map(f, Z))
      f = interp1d(Z, V1dRms, kind='linear', bounds_error=False, fill_value='extrapolate')
      
      print "Expected v1d_rms = "+str(np.mean(np.array(map(f, self.Z))))+" km/s"
      print "Expected v1d_rms(z=z_mean) = "+str(f(np.mean(self.Z)))+" km/s"
      print "RMS v_r, v_theta, v_phi = "+str(np.std(self.vR))+", "+str(np.std(self.vTheta))+", "+str(np.std(self.vPhi))+" km/s"



   ##################################################################################
   ##################################################################################


   def plotHistograms(self):
      z0 = np.mean(self.Z)
      s2v1d = self.U.v3dRms(0., z0, W3d_sth)**2 / 3.
      
      # redshifts
      path = self.pathFig+"/hist_z.pdf"
      myHistogram(self.Z, nBins=71, lim=(0., 1.), path=path, nameLatex=r'$z$', semilogy=True)
      
      # spherical velocities
      path = self.pathFig+"/hist_vr.pdf"
      myHistogram(self.vR, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_r$ [km/s]', doGauss=True)
      path = self.pathFig+"/hist_vtheta.pdf"
      myHistogram(self.vTheta, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_\theta$ [km/s]', doGauss=True)
      path = self.pathFig+"/hist_vphi.pdf"
      myHistogram(self.vPhi, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_\phi$ [km/s]', doGauss=True)
      
      # cartesian velocities
      path = self.pathFig+"/hist_vx.pdf"
      myHistogram(self.vX, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_x$ [km/s]', doGauss=True)
      path = self.pathFig+"/hist_vy.pdf"
      myHistogram(self.vY, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_y$ [km/s]', doGauss=True)
      path = self.pathFig+"/hist_vz.pdf"
      myHistogram(self.vZ, nBins=71, lim=(-1000., 1000.), S2Theory=[s2v1d], path=path, nameLatex=r'$v_z$ [km/s]', doGauss=True)
      
      # stellar masses
      path = self.pathFig+"/hist_mstellar.pdf"
      myHistogram(self.Mstellar, nBins=71, path=path, nameLatex=r'$M_\star$ [M$_\odot$]', semilogx=True, semilogy=True)

      # virial masses
      path = self.pathFig+"/hist_mvir.pdf"
      myHistogram(self.Mvir, nBins=71, path=path, nameLatex=r'$M_\text{vir}$ [M$_\odot$]', semilogx=True, semilogy=True)
      
      # comoving virial radius
      # need masses in Msun/h
      Par = zip(self.Mvir*self.U.bg.h, self.Z)
      f = lambda par: self.U.frvir(par[0], par[1])   # in: Msun/h, out: Mpc/h
      Rvir = np.array(map(f, Par))  # in Mpc/h
      #Rvir /= self.U.bg.h  # Mpc
      path = self.pathFig+"/hist_rvir.pdf"
      myHistogram(Rvir/self.U.bg.h, nBins=71, path=path, nameLatex=r'$R_\text{vir}$ [Mpc]', semilogx=True, semilogy=True)
      
      # virial angular radius
      Chi = np.array(map(self.U.bg.comoving_distance, self.Z)) # [Mpc/h]
      Thetavir = Rvir / Chi   # [rad]
      path = self.pathFig+"/hist_thetavir.pdf"
      x = Thetavir * (180.*60./np.pi)  # [arcmin]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\theta_\text{vir}$ [arcmin]', semilogx=True, semilogy=True)
      
      # integrated tau [arcmin^2]
      path = self.pathFig+"/hist_integratedtau.pdf"
      x = self.integratedTau * (180.*60./np.pi)**2 # [arcmin^2]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; \tau$ [arcmin$^2$]', semilogx=True, semilogy=True)
      
      # mean tau within Rvir [dimless]
      path = self.pathFig+"/hist_meantauvir.pdf"
      x = self.integratedTau / (np.pi * Thetavir**2) # [dimless]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; \tau / \left( \pi \theta_\text{vir} \right)$ [dimless]', semilogx=True, semilogy=True)

      # expected kSZ [muK*arcmin^2]
      path = self.pathFig+"/hist_ksz.pdf"
      x = self.integratedKSZ * (180.*60./np.pi)**2 # [muK*arcmin^2]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; \delta T_\text{kSZ}$ [$\mu$K.arcmin$^2$]', doGauss=True, semilogy=True)

      # mean kSZ within Rvir [muK]
      path = self.pathFig+"/hist_meankszvir.pdf"
      x = self.integratedKSZ / (np.pi * Thetavir**2) # [muK]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; \delta T_\text{kSZ} / \left( \pi \theta_\text{vir} \right)$ [$\mu$K]', doGauss=True, semilogy=True)

      # expected Y [arcmin^2]
      path = self.pathFig+"/hist_y.pdf"
      x = self.integratedY * (180.*60./np.pi)**2 # [arcmin^2]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; y_\text{tSZ}$ [arcmin$^2$]', semilogx=True, semilogy=True)

      # mean Y within Rvir [dimless]
      path = self.pathFig+"/hist_meanyvir.pdf"
      x = self.integratedY / (np.pi * Thetavir**2) # [dimless]
      myHistogram(x, nBins=71, path=path, nameLatex=r'$\int d^2\theta \; y_\text{tSZ} / \left( \pi \theta_\text{vir} \right)$ [dimless]', semilogx=True, semilogy=True)


      # displacements?
      # displacements?


##################################################################################
##################################################################################


#   def generateMockMaps(self, carMap, sigma=None, depixwin=False, test=False):
#      print "- Generate mock maps"
#      # map of pixel areas
#      pixSizeMap = carMap.pixsizemap()
#      # input for the counts
#      srcsCount = np.zeros((self.nObj, 3))
#      srcsCount[:,0] = self.DEC.copy() * np.pi/180.   # [rad]
#      srcsCount[:,1] = self.RA.copy() * np.pi/180. # [rad]
#      srcsCount[:,2] = 1.
#      # input for the LOS velocity
#      srcsVel = srcsCount.copy()
#      srcsVel[:,2] = - self.vR / 3.e5 
#
#      import pointsrcs
#
#      # Dirac profiles
#      countDirac = pointsrcs.sim_srcs(carMap.shape, carMap.wcs, srcsCount, 1.e-5*np.pi/(180.*60.))
#      print(np.sum(np.abs(countDirac)))
#      velDirac = pointsrcs.sim_srcs(carMap.shape, carMap.wcs, srcsVel, 1.e-5*np.pi/(180.*60.))
#      # normalize to integrate to 1 over angles in [muK*arcmin^2]
#      countDirac /= pixSizeMap * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
#      velDirac /= pixSizeMap * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
#      # normalize the mock maps, such that:
#      # int dOmega count = 1 [muK*arcmin^2]
#      # int dOmega vel = -(v/c) * sigma(v/c) [muK*arcmin^2]
#      # the reason for the vel normalization is that the kSZ estimator correlates with v/c,
#      # then divides by sigma^2(v/c), then re-multiplies by sigma(v/c),
#      # so that the estimated kSZ has the right amplitude.
#      # This way, the estimated tSZ and kSZ should converge to 1 muK*arcmin^2,
#      # and will be easily comparable to the theory curve.
#      velDirac /= np.std(self.vR / 3.e5)
#      # save the maps
#      enmap.write_map(self.pathOut+"mock_count_dirac_car.fits", countDirac)
#      enmap.write_map(self.pathOut+"mock_vel_dirac_car.fits", velDirac)
#      
#      if sigma is not None:
#         countGauss = pointsrcs.sim_srcs(carMap.shape, carMap.wcs, srcsCount, sigma*np.pi/(180.*60.))
#         print(np.sum(np.abs(countGauss)))
#         velGauss = pointsrcs.sim_srcs(carMap.shape, carMap.wcs, srcsVel, sigma*np.pi/(180.*60.))
#         # normalize to integrate to 1 over angles in [muK*arcmin^2]
#         countGauss /= pixSizeMap * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
#         velGauss /= pixSizeMap * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
#         # The amplitude given to the Gaussian was the peak amplitude:
#         # convert from peak amplitude to Gaussian normalization
#         countGauss /= 2. * np.pi * (sigma * np.pi / (180. * 60.))**2
#         velGauss /= 2. * np.pi * (sigma * np.pi / (180. * 60.))**2
#         # normalize the mock maps, such that:
#         # int dOmega count = 1 [muK*arcmin^2]
#         # int dOmega vel = -(v/c) * sigma(v/c) [muK*arcmin^2]
#         # the reason for the vel normalization is that the kSZ estimator correlates with v/c,
#         # then divides by sigma^2(v/c), then re-multiplies by sigma(v/c),
#         # so that the estimated kSZ has the right amplitude.
#         # This way, the estimated tSZ and kSZ should converge to 1 muK*arcmin^2,
#         # and will be easily comparable to the theory curve.
#         velGauss /= np.std(self.vR / 3.e5)
#         # save the maps
#         enmap.write_map(self.pathOut+"mock_count_gauss_car.fits", countGauss)
#         enmap.write_map(self.pathOut+"mock_vel_gauss_car.fits", velGauss)




   def generateMockMaps(self, carMap, sigma=None, test=False):
      """Generate mock maps with 1 at the pixel location of each  object, 0 everywhere else.
      If sigma [arcmin] is specified, produces also Gaussian smoothed versions,
      normalized such that   int d^2theta profile = 1, where theta is in [rad].
      If depixwin==True, the Gaussian profile map is deconvolved with one power
      of the pixel window function. This messes up the individual profile, 
      but ensures that the stacked profile is correct, by correcting the fact that 
      the galaxy profiles were placed at the center of the nearest pixel, 
      as opposed to the exact position within the pixel.
      This operation is not done on the Dirac maps, since the slight miscentering
      has no observable impact there.
      Assumes that carMap has shape [nX, nY], ie it is not a T,Q,U map, just a T map.
      """
      print "- Generate mock maps"
      tStart = time()
      # create empty maps
      countDirac = carMap.copy()
      countDirac[:,:] = 0.
      velDirac = countDirac.copy()
      # get map of exact pixel sizes
      pixSizeMap = countDirac.pixsizemap()
      # get map of ra and dec, just to check
      posmap = countDirac.posmap()

      for iObj in range(self.nObj):
#      for iObj in range(10):
         if iObj%100000==0:
            print "    -", iObj
         # object coordinates [deg]
         ra = self.RA[iObj]
         dec = self.DEC[iObj]
         # coordinates in [rad]
         sourcecoord = np.array([dec, ra]) * np.pi/180.
         # find pixel indices (float) corresponding to ra, dec
         iY, iX = enmap.sky2pix(countDirac.shape, countDirac.wcs, sourcecoord, safe=True, corner=False)
         if test:
            print 'ra, dec =', ra, dec, iY, iX
            print countDirac.shape
         # Check that the object is within the map boundaries
         # before rounding the indices
         if iX>=0 and iX<=(countDirac.shape[1]-1) and iY>=0 and iY<=(countDirac.shape[0]-1):
             # nearest pixel
             # watch out for difference round VS np.round!
             jY = np.int(round(iY))
             jX = np.int(round(iX))
             if test:
               print("Object "+str(iObj)+" overlaps")
             # fill the pixel
             countDirac[jY, jX] = 1.
             velDirac[jY, jX] = - self.vR[iObj] / 3.e5   # v_r/c  [dimless]
             # check that I filled the right pixel
             if countDirac.at(sourcecoord, prefilter=False, mask_nan=False, order=0)<>1:
                print "Filled the wrong pixel for  object", iObj
                print "wanted ra, dec=", ra, dec # [deg]
                print "chosen closest ra, dec=", posmap[::-1, jY, jX] * 180./np.pi  # [deg]
                print "difference in arcmin=", (posmap[::-1, jY, jX] * 180./np.pi - np.array([ra, dec]))*60.  # residual in [arcmin]
                print "ra index=", iX, jX, np.int(np.round(iX)), countDirac.shape[1]
                print "dec index=", iY, jY, np.int(np.round(iY)), countDirac.shape[0]

             # normalize to integrate to 1 over angles in [muK*arcmin^2]
             countDirac[jY, jX] /= pixSizeMap[jY, jX] * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
             velDirac[jY, jX] /= pixSizeMap[jY, jX] * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
             
      # normalize the mock maps, such that:
      # int dOmega count = 1 [muK*arcmin^2]
      # int dOmega vel = -(v/c) / sigma(v/c) [muK*arcmin^2]
      # the reason for the vel normalization is that the kSZ estimator correlates with v/c,
      # then divides by sigma^2(v/c), then re-multiplies by sigma(v/c),
      # so that the estimated kSZ has the right amplitude.
      # This way, the estimated tSZ and kSZ should converge to 1 muK*arcmin^2,
      # and will be easily comparable to the theory curve.
#      countDirac /= countDirac.pixsize() * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
#      velDirac /= velDirac.pixsize() * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
      velDirac /= np.std(self.vR / 3.e5)

      # save the maps
      print "Saving maps to:"
      print self.pathOut+"mock_count_dirac_car.fits"
      print self.pathOut+"mock_vel_dirac_car.fits"
      enmap.write_map(self.pathOut+"mock_count_dirac_car.fits", countDirac)
      enmap.write_map(self.pathOut+"mock_vel_dirac_car.fits", velDirac)
      

      # For the Gaussian profiles, use pixell.pointsrcs.sim_srcs,
      # which places the object at the exact position, rather than the center
      # of the closest pixel.
      if sigma is not None:


         # input for the counts
         srcsCount = np.zeros((self.nObj, 3))
         srcsCount[:,0] = self.DEC.copy() * np.pi/180.   # [rad]
         srcsCount[:,1] = self.RA.copy() * np.pi/180. # [rad]
         srcsCount[:,2] = 1.
         # input for the LOS velocity
         srcsVel = srcsCount.copy()
         srcsVel[:,2] = - self.vR / 3.e5 

         import pointsrcs

         countGauss = pointsrcs.sim_srcs(carMap.shape, carMap.wcs, srcsCount, sigma*np.pi/(180.*60.))
         velGauss = pointsrcs.sim_srcs(carMap.shape, carMap.wcs, srcsVel, sigma*np.pi/(180.*60.))
#         # normalize to integrate to 1 over angles in [muK*arcmin^2]
#         countGauss /= pixSizeMap * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
#         velGauss /= pixSizeMap * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
         # The amplitude given to the Gaussian was the peak amplitude:
         # convert from peak amplitude to Gaussian normalization
         countGauss /= 2. * np.pi * sigma**2 # (sigma * np.pi / (180. * 60.))**2
         velGauss /= 2. * np.pi * sigma**2   #(sigma * np.pi / (180. * 60.))**2
         # normalize the mock maps, such that:
         # int dOmega count = 1 [muK*arcmin^2]
         # int dOmega vel = -(v/c) / sigma(v/c) [muK*arcmin^2]
         # the reason for the vel normalization is that the kSZ estimator correlates with v/c,
         # then divides by sigma^2(v/c), then re-multiplies by sigma(v/c),
         # so that the estimated kSZ has the right amplitude.
         # This way, the estimated tSZ and kSZ should converge to 1 muK*arcmin^2,
         # and will be easily comparable to the theory curve.
         velGauss /= np.std(self.vR / 3.e5)
         # save the maps
         print "Saving maps to:"
         print self.pathOut+"mock_count_gauss_car.fits"
         print self.pathOut+"mock_vel_gauss_car.fits"
         enmap.write_map(self.pathOut+"mock_count_gauss_car.fits", countGauss)
         enmap.write_map(self.pathOut+"mock_vel_gauss_car.fits", velGauss)

      tStop = time()
      print("Took "+str((tStop-tStart)/60.)+" min")







#   def generateMockMaps(self, carMap, sigma=None, depixwin=False, test=False):
#      """Generate mock maps with 1 at the pixel location of each  object, 0 everywhere else.
#      If sigma [arcmin] is specified, produces also Gaussian smoothed versions,
#      normalized such that   int d^2theta profile = 1, where theta is in [rad].
#      If depixwin==True, the Gaussian profile map is deconvolved with one power
#      of the pixel window function. This messes up the individual profile, 
#      but ensures that the stacked profile is correct, by correcting the fact that 
#      the galaxy profiles were placed at the center of the nearest pixel, 
#      as opposed to the exact position within the pixel.
#      This operation is not done on the Dirac maps, since the slight miscentering
#      has no observable impact there.
#      Assumes that carMAp has shape [nX, nY], ie it is not a T,Q,U map, just a T map.
#      """
#      print "- Generate mock maps"
#      tStart = time()
#      # create empty maps
#      countDirac = carMap.copy()
#      countDirac[:,:] = 0.
#      velDirac = countDirac.copy()
#
#      # get map of exact pixel sizes
#      pixSizeMap = countDirac.pixsizemap()
#
#      # get map of ra and dec, just to check
#      posmap = countDirac.posmap()
#
#      for iObj in range(self.nObj):
##      for iObj in range(10):
#         if iObj%100000==0:
#            print "    -", iObj
#         # object coordinates [deg]
#         ra = self.RA[iObj]
#         dec = self.DEC[iObj]
#         # coordinates in [rad]
#         sourcecoord = np.array([dec, ra]) * np.pi/180.
#         
#         # find pixel indices (float) corresponding to ra, dec
#         iY, iX = enmap.sky2pix(countDirac.shape, countDirac.wcs, sourcecoord, safe=True, corner=False)
#         if test:
#            print 'ra, dec =', ra, dec, iY, iX
#            print countDirac.shape
#
#         # Check that the object is within the map boundaries
#         # before rounding the indices
#         if iX>=0 and iX<=(countDirac.shape[1]-1) and iY>=0 and iY<=(countDirac.shape[0]-1):
#             # nearest pixel
#             # watch out for difference round VS np.round!
#             jY = np.int(round(iY))
#             jX = np.int(round(iX))
#             
#             if test:
#               print("Object "+str(iObj)+" overlaps")
#
#             # fill the pixel
#             countDirac[jY, jX] = 1.
#             velDirac[jY, jX] = - self.vR[iObj] / 3.e5   # v_r/c  [dimless]
#
#             # check that I filled the right pixel
#             if countDirac.at(sourcecoord, prefilter=False, mask_nan=False, order=0)<>1:
#                print "Filled the wrong pixel for  object", iObj
#                print "wanted ra, dec=", ra, dec # [deg]
#                print "chosen closest ra, dec=", posmap[::-1, jY, jX] * 180./np.pi  # [deg]
#                print "difference in arcmin=", (posmap[::-1, jY, jX] * 180./np.pi - np.array([ra, dec]))*60.  # residual in [arcmin]
#                print "ra index=", iX, jX, np.int(np.round(iX)), countDirac.shape[1]
#                print "dec index=", iY, jY, np.int(np.round(iY)), countDirac.shape[0]
#                
#
#             # normalize to integrate to 1 over angles in [muK*arcmin^2]
#             countDirac[jY, jX] /= pixSizeMap[jY, jX] * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
#             velDirac[jY, jX] /= pixSizeMap[jY, jX] * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
#             
#
#
#      # normalize the mock maps, such that:
#      # int dOmega count = 1 [muK*arcmin^2]
#      # int dOmega vel = -(v/c) * sigma(v/c) [muK*arcmin^2]
#      # the reason for the vel normalization is that the kSZ estimator correlates with v/c,
#      # then divides by sigma^2(v/c), then re-multiplies by sigma(v/c),
#      # so that the estimated kSZ has the right amplitude.
#      # This way, the estimated tSZ and kSZ should converge to 1 muK*arcmin^2,
#      # and will be easily comparable to the theory curve.
##      countDirac /= countDirac.pixsize() * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
##      velDirac /= velDirac.pixsize() * (180.*60./np.pi)**2 # divide by pixel area in arcmin^2 
#      velDirac /= np.std(self.vR / 3.e5)
#
##      # undo one power of the window function,
##      # to undo the fact that I placed the galaxies at the center of the nearest pixel,
##      # as opposed to their exact position within the pixel
##      if depixwin:
##         countDirac = enmap.apply_window(countDirac, pow=-1)
##         velDirac = enmap.apply_window(velDirac, pow=-1)
#
#      # save the maps
#      enmap.write_map(self.pathOut+"mock_count_dirac_car.fits", countDirac)
#      enmap.write_map(self.pathOut+"mock_vel_dirac_car.fits", velDirac)
#      
#      if sigma is not None:
#         # convolve maps with a Gaussian  profile of given sigma (not fwhm)
#         sigma *= np.pi/180./60.  # convert from arcmin to [rad]
#         countGauss = enmap.smooth_gauss(countDirac, sigma)
#         velGauss = enmap.smooth_gauss(velDirac, sigma)
#
#         # undo one power of the window function,
#         # to undo the fact that I placed the galaxies at the center of the nearest pixel,
#         # as opposed to their exact position within the pixel
#         # individual profiles will look weird, but the stack should then be correct
#         if depixwin:
#            countGauss = enmap.apply_window(countGauss, pow=-1)
#            velGauss = enmap.apply_window(velGauss, pow=-1)
#         
#         enmap.write_map(self.pathOut+"mock_count_gauss_car.fits", countGauss)
#         enmap.write_map(self.pathOut+"mock_vel_gauss_car.fits", velGauss)
#
#      tStop = time()
#      print("Took "+str((tStart-tStop)/60.)+" min")
