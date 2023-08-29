import numpy as np

a = np.arange(5)
a[np.logical_or.reduce([a[:] < 3])] = 1

print "a", a
