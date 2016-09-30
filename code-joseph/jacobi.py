from matplotlib.pyplot import *
from numpy import *
S = loadtxt("eigenvectors.data")
E = loadtxt("eigenvalues.data")
index = argmin(E)
vector = S[:,index]

plot(range(len(vector)),vector*vector)
show()
