import numpy as np

def moveaxis(a, o, n):
	if o < 0: o = o+a.ndim
	if n < 0: n = n+a.ndim
	if n <= o: return np.rollaxis(a, o, n)
	else: return np.rollaxis(a, o, n+1)

def rotmatrix(ang, raxis, axis=0):
	"""Construct a 3d rotation matrix representing a rotation of
	ang degrees around the specified rotation axis raxis, which can be "x", "y", "z"
	or 0, 1, 2. If ang is a scalar, the result will be [3,3]. Otherwise,
	it will be ang.shape + (3,3)."""
	ang  = np.asarray(ang)
	raxis = raxis.lower()
	c, s = np.cos(ang), np.sin(ang)
	R = np.zeros(ang.shape + (3,3))
	if   raxis == 0 or raxis == "x": R[...,0,0]=1;R[...,1,1]= c;R[...,1,2]=-s;R[...,2,1]= s;R[...,2,2]=c
	elif raxis == 1 or raxis == "y": R[...,0,0]=c;R[...,0,2]= s;R[...,1,1]= 1;R[...,2,0]=-s;R[...,2,2]=c
	elif raxis == 2 or raxis == "z": R[...,0,0]=c;R[...,0,1]=-s;R[...,1,0]= s;R[...,1,1]= c;R[...,2,2]=1
	else: raise ValueError("Rotation axis %s not recognized" % raxis)
	return moveaxis(R, 0, axis)

def ang2rect(angs, zenith=True, axis=0):
	"""Convert a set of angles [{phi,theta},...] to cartesian
	coordinates [{x,y,z},...]. If zenith is True (the default),
	the theta angle will be taken to go from 0 to pi, and measure
	the angle from the z axis. If zenith is False, then theta
	goes from -pi/2 to pi/2, and measures the angle up from the xy plane."""
	phi, theta = moveaxis(angs, axis, 0)
	ct, st, cp, sp = np.cos(theta), np.sin(theta), np.cos(phi), np.sin(phi)
	if zenith: res = np.array([st*cp,st*sp,ct])
	else:      res = np.array([ct*cp,ct*sp,st])
	return moveaxis(res, 0, axis)

def rect2ang(rect, zenith=True, axis=0):
	"""The inverse of ang2rect."""
	x,y,z = moveaxis(rect, axis, 0)
	r     = (x**2+y**2)**0.5
	phi   = np.arctan2(y,x)
	if zenith: theta = np.arctan2(r,z)
	else:      theta = np.arctan2(z,r)
	return moveaxis(np.array([phi,theta]), 0, axis)

def euler_mat(euler_angles, kind="zyz"):
	"""Defines the rotation matrix M for a ABC euler rotation,
	such that M = A(alpha)B(beta)C(gamma), where euler_angles =
	[alpha,beta,gamma]. The default kind is ABC=ZYZ."""
	alpha, beta, gamma = euler_angles
	R1 = rotmatrix(gamma, kind[2])
	R2 = rotmatrix(beta,  kind[1])
	R3 = rotmatrix(alpha, kind[0])
	return np.einsum("...ij,...jk->...ik",np.einsum("...ij,...jk->...ik",R3,R2),R1)

def euler_rot(euler_angles, coords, kind="zyz"):
	coords = np.asarray(coords)
	co     = coords.reshape(2,-1)
	M      = euler_mat(euler_angles, kind)
	rect   = ang2rect(co, False)
	rect   = np.einsum("...ij,j...->i...",M,rect)
	co     = rect2ang(rect, False)
	return co.reshape(coords.shape)

def recenter(angs, center):
	"""recenter(angs, [from_ra, from_dec]) or recenter(angs, [from_ra, from_dec, to_ra, to_dec]).
	In the first form, performs a coordinate rotation such that a point at (form_ra, from_to)
	ends up at the north pole. In the second form, that point is instead put at (to_ra, to_dec)."""
	# Performs the rotation E(0,-theta,-phi). Originally did
	# E(phi,-theta,-phi), but that is wrong (at least for our
	# purposes), as it does not preserve the relative orientation
	# between the boresight and the sun. For example, if the boresight
	# is at the same elevation as the sun but 10 degrees higher in az,
	# then it shouldn't matter what az actually is, but with the previous
	# method it would.
	#
	# Now supports specifying where to recenter by specifying center as
	# lon_from,lat_from,lon_to,lat_to
	if len(center) == 4: ra0, dec0, ra1, dec1 = center
	elif len(center) == 2: ra0, dec0, ra1, dec1 = center[0], center[1], 0, np.pi/2
	return euler_rot([ra1,dec0-dec1,-ra0], angs, kind="zyz")

def decenter(angs, center):
	"""Inverse operation of recenter."""
	if len(center) == 4: ra0, dec0, ra1, dec1 = center
	elif len(center) == 2: ra0, dec0, ra1, dec1 = center[0], center[1], 0, np.pi/2
	return euler_rot([ra0,dec1-dec0,-ra1],  angs, kind="zyz")
