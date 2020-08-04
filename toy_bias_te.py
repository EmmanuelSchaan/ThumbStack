import numpy as np
import matplotlib.pyplot as plt


nSamples = 100000

# tSZ data
yTrue = 16.
sY = 0.1 * yTrue
ySamples = np.random.normal(loc=yTrue, scale=sY, size=nSamples)

# kSZ data
tauTrue = 1.e4
sTau = 0.19 * tauTrue
tauSamples = np.random.normal(loc=tauTrue, scale=sTau, size=nSamples)

# temperature estimator
tTrue = yTrue / tauTrue

tSamples = ySamples / tauSamples

print "mean should be", tTrue
print "the fractional bias is", np.mean(tSamples)/tTrue-1.
