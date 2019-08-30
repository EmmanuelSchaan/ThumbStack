





# create a map full of zeros
shape, wcs = enmap.geometry(np.array([[0.,-np.pi/2.],[2.*np.pi,np.pi/2.]]), res=0.25 * np.pi/180./60., proj='car')
carMap = enmap.zeros(shape, wcs)

# coordinates of point to set to 1 [radians]
ra = -np.pi/2. +  0.0003
dec = 0.

# find the pixel indices
iY, iX = enmap.sky2pix(carMap.shape, carMap.wcs, [dec, ra], safe=True, corner=False)
print iY, iX

#  make them integers, and in the proper range
iY = np.int(round(iY) % carMap.shape[0])
iX = np.int(round(iX) % carMap.shape[1])

# set the map value to 1
carMap[iY, iX] = 1.

# check that I filled in the right pixel
sourcecoord = np.array([dec, ra])   # convert from degrees to radians
# use nearest neighbor interpolation
print carMap.at(sourcecoord, prefilter=False, mask_nan=False, order=0)



   def generateMockMap(self, carMap):
      """Not  yet  working!!!
      """
      
      # put 1 in the closest pixel to each object
      countDirac = carMap.copy()
      countDirac[:,:] = 0.
      
#      for iObj in range(self.nObj):
      for iObj in range(1):
         # object coordinates
         ra = self.RA[iObj]
         dec = self.DEC[iObj]
         
         # find pixel indices (float) corresponding to ra, dec
         iY, iX = enmap.sky2pix(carMap.shape, carMap.wcs, [dec, ra], safe=True, corner=False)
#         iY = round(iY) % carMap.shape[0]
#         iX = round(iX) % carMap.shape[1]

         # nearest pixel
         iY = np.int(round(iY) % carMap.shape[0])
         iX = np.int(round(iX) % carMap.shape[1])
         print iY, iX

#         # nearest pixel
#         iY = np.int(round(iY))
#         iX = np.int(round(iX))

         countDirac[iY, iX] = 1.
         if self.sky2map(ra, dec, countDirac)==1.:
            print "it worked"
         else:
            print "Problem!!!, filled in the wrong pixel"
#      enmap.write_map(self.pathOut + "/mock_count_dirac_cea.fits")


      

      
#      opos = carMap.posmap()
#      # output map position [{dec,ra},ny,nx]
#      dec = opos[0,:,:]
#      ra = opos[1,:,:]
#
#
#
#
##         raResRad = (ra[0,:] - self.RA[iObj]*(np.pi/180.)) % (2.*np.pi)
##         decResRad = (dec[:,0] - self.DEC[iObj]*(np.pi/180.)) % (2.*np.pi)
##         iX =  np.argmin(raResRad**2)
##         iY =  np.argmin(decResRad**2)
#
#
#         raResRad = (ra - self.RA[iObj]*(np.pi/180.)) % (2.*np.pi)
#         decResRad = (dec - self.DEC[iObj]*(np.pi/180.)) % (2.*np.pi)
#         dist2 = raResRad**2 + decResRad**2
#         iY, iX = np.unravel_index(np.argmin(dist2, axis=None), dist2.shape)
#
#
#
#         countDirac[iY, iX] = 1.
#         if self.sky2map(self.RA[iObj], self.DEC[iObj], countDirac)<>1.:
#            print "Problem!!!, filled in the wrong pixel"
#      enmap.write_map(self.pathOut + "/mock_count_dirac_cea.fits")


#      self.pathOut + "/mock_count_gauss_sigma1.5arcmin_cea.fits"
#      self.pathOut + "/mock_v_gauss_sigma1.5arcmin_cea.fits"
#      pass








   def sky2map(self, ra, dec, map):
      '''Gives the map value at coordinates (ra, dec).
      ra, dec in degrees.
      Uses nearest neighbor, no interpolation.
      Will return 0 if the coordinates requested are outside the map
      '''
      # interpolate the map to the given sky coordinates
      sourcecoord = np.array([dec, ra]) * utils.degree   # convert from degrees to radians
      # use nearest neighbor interpolation
      return map.at(sourcecoord, prefilter=False, mask_nan=False, order=0)
