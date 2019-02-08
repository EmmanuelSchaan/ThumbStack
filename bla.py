# check the hit map as a function of dec, at fixed ra = 0
def sky2map(dec, ra=45.):
   yPix,xPix = enmap.sky2pix(bigmap.shape, bigmap.wcs, np.array([dec,ra])*utils.degree)
   # round to closest pixel
   yPix = np.int(np.round(yPix)) % bigmap.shape[0]
   xPix = np.int(np.round(xPix)) % bigmap.shape[1]
   # evaluate the map at the desired pixel
   return bigmap[yPix, xPix]


# check the hit map as a function of dec, at fixed ra = 0
def sky2map(dec, ra=45.):
   # interpolate the map to the given sky coordinates
   sourcecoord = np.array([dec, ra])*utils.degree
   return bigmap.at(sourcecoord, prefilter=False, mask_nan=False)
