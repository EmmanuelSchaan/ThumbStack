import numpy as np
import matplotlib.pyplot as plt

#plt.switch_backend('Agg')
#matplotlib.use('GTK')
#plt.switch_backend('GTK')


x = np.linspace(0.,2.*np.pi, 101)

plt.plot(x, np.cos(x))

plt.show()
