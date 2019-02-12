import numpy as np
import matplotlib.pyplot as plt

#plt.switch_backend('Agg')
#matplotlib.use('GTK')
plt.switch_backend('GTK')


x = np.linspace(0.,2.*np.pi, 101)

plt.plot(x, np.cos(x))

plt.show()



















#      self.filtMap = np.genfromtxt(self.pathOut+"/filtmap.txt")
#      self.filtMask = np.genfromtxt(self.pathOut+"/filtmask.txt")
#      self.filtVar = np.genfromtxt(self.pathOut+"/filtvar.txt")
#      self.diskArea = np.genfromtxt(self.pathOut+"/diskarea.txt")
#
#
#
#
#
#np.where((ts.Catalog.RA>0.)*(ts.Catalog.RA<10.)*(ts.Catalog.DEC<10.)*(ts.Catalog.DEC>0.))
#iObj = 9770
#
#
#np.where(ts.filtMask[:,0]>=1.e-5)
