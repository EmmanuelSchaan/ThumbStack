



dec = 20 ; ra = 180
coords = np.deg2rad(np.array((dec,ra)))
ypix,xpix = enmap.sky2pix(shape,wcs,coords)

bigmap[ypix, xpix]#???? Is this the right order?


ra = 20. # deg
dec = 3. # deg
bigmap.at(np.array([ra, dec])*utils.degree)


