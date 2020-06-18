import numpy as np
import matplotlib.pyplot as plt





pathInDir150 = "./output/thumbstack/cmass_kendrick_pactf150daynight20200228maskgal60r2/"
pathInDir90 = "./output/thumbstack/cmass_kendrick_pactf90daynight20200228maskgal60r2/"


###################################################################################
# kSZ

path = pathInDir150 + "cov_diskring_ksz_varweight_joint_cmass_kendrick_pactf150daynight20200228maskgal60r2_cmass_kendrick_pactf90daynight20200228maskgal60r2_bootstrap.txt"
covKszJoint = np.genfromtxt(path)
path = pathInDir150 + "cov_diskring_ksz_varweight_bootstrap.txt"
covKsz150 = np.genfromtxt(path)
path = pathInDir90 + "cov_diskring_ksz_varweight_bootstrap.txt"
covKsz90 = np.genfromtxt(path)

s150 = np.sqrt(np.diag(covKsz150))
s90 = np.sqrt(np.diag(covKsz90))

s150FromJoint = np.sqrt(np.diag(covKszJoint))[:len(s150)]
s90FromJoint = np.sqrt(np.diag(covKszJoint))[len(s150):]


fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.axhline(0.)
ax.plot(s150FromJoint/s150 - 1., label=r'150 from joint / 150 - 1')
ax.plot(s90FromJoint/s90 - 1., label=r'90 from joint / 90 - 1')
#
ax.legend(fontsize='x-small')
ax.set_title(r'kSZ')

plt.show()



###################################################################################
# tSZ

path = pathInDir150 + "cov_diskring_tsz_varweight_joint_cmass_kendrick_pactf150daynight20200228maskgal60r2_cmass_kendrick_pactf90daynight20200228maskgal60r2_bootstrap.txt"
covTszJoint = np.genfromtxt(path)
path = pathInDir150 + "cov_diskring_tsz_varweight_bootstrap.txt"
covTsz150 = np.genfromtxt(path)
path = pathInDir90 + "cov_diskring_tsz_varweight_bootstrap.txt"
covTsz90 = np.genfromtxt(path)

s150 = np.sqrt(np.diag(covTsz150))
s90 = np.sqrt(np.diag(covTsz90))

s150FromJoint = np.sqrt(np.diag(covTszJoint))[:len(s150)]
s90FromJoint = np.sqrt(np.diag(covTszJoint))[len(s150):]


fig=plt.figure(0)
ax=fig.add_subplot(111)
#
ax.axhline(0.)
ax.plot(s150FromJoint/s150 - 1., label=r'150 from joint / 150 - 1')
ax.plot(s90FromJoint/s90 - 1., label=r'90 from joint / 90 - 1')
#
ax.legend(fontsize='x-small')
ax.set_title(r'tSZ')

plt.show()
















