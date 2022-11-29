from BoxFinder import *
import numpy as np
import sys

# Usage: input file, 0 < bandwidth coefficient <= 1, number of iterations

X = np.loadtxt(sys.argv[1])
X = X[:,0:-1]

mins = np.min(X, axis=0)
maxs = np.max(X, axis=0)
size = (maxs - mins)*float(sys.argv[2])

non_fixed = size != 0
tmp_center1 = np.array(X[0,:])
tmp_center2 = np.array(X[0,:])
X = X[:, non_fixed]
size = size[non_fixed]

if size.any():
	b = BoxFinder(X, size=size, n_iter=int(sys.argv[3]))
	max_ctr, max_cov, base_ctr, base_cov = b.fit()
else:
	max_ctr = np.array([])
	base_ctr = np.array([])
	max_cov = 1
	base_cov = 1

tmp_center1[non_fixed] = max_ctr
tmp_center2[non_fixed] = base_ctr
max_ctr = tmp_center1
base_ctr = tmp_center2

# Print center of the maximal box
print("\t".join(map(str, max_ctr)))
#print("\t".join(map(str, base_ctr)))

f = open("log.txt", "w")
f.write("Center of the base box: {}\n".format("\t".join(map(str, base_ctr))))
f.write("Number of points in base box: {}\n".format(base_cov))
f.write("Number of points in optimal box: {}\n".format(max_cov))
f.close()
