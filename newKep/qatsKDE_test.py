from getdbData import qatsResults_3d
from plotQats import meshgrid_3d
from astroML.density_estimation import KDE
from mpl_toolkits.mplot3d import Axes3D
import numpy as num
import pylab

def scatterDens(x, y, z, out):
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(x, y, z, c=out)

def contourDens(x, y, z, out, nBins = 5):
    fig = pylab.figure()
    ax = fig.add_subplot(111, projection='3d')
    bins = num.linspace(out.min(), out.max(), nBins + 1)
    # plot outermost shell where bin criterion is satisfied
    idx = num.where((out > bins[1]) & (out < bins[2]))
    ax.scatter(x[idx], y[idx], z[idx])
    #for bin in range(nBins):
        #idx = num.where(out > bins[bin])
        #color = pylab.matplotlib.colors.ColorConverter().to_rgb(str(bins[bin] / bins.max()))
        #ax.scatter(x[idx], y[idx], z[idx])
        #ax.plot_surface(x[idx], y[idx], z[idx], color=color)

qMatrix = qatsResults_3d()
densEstimator = KDE(h=5)

"""
Need qMatrix shape to be (Npoints, 3) or (Npoints, 4)
Each row will be a coordinate specifying a point in 4-space
We want the density to be largest near signal peaks
(so points with high SNR)
If indices are used to specify the first three coordinates, the main factor in deciding density SHOULD be the snr (since indices are all evenly spaced)
"""
print qMatrix.shape
wIdx = num.arange(qMatrix.shape[0])
dIdx = num.arange(qMatrix.shape[1])
pBinIdx = num.arange(qMatrix.shape[2])
P, D, W = meshgrid_3d(pBinIdx, dIdx, wIdx)
w = W.ravel()
d = D.ravel()
p = P.ravel()
# qMatrix[w, d, p] is equivalent to qMatrix.ravel()
# take the transpose of this array (need each ROW to be a set of coordinates, not each column)
kde_qMatrix = num.vstack( (w, d, p, qMatrix[w, d, p]) ).T

# fit the density estimator to the qats matrix
densEstimator.fit(kde_qMatrix)
#import pdb; pdb.set_trace()

# evaluate the KDE on the qats matrix
# apparently evaluation is extremely expensive...
# so we cut down the array to a manageable size
cheapSlice = 4096
print kde_qMatrix[::cheapSlice].shape

if False:
    print 'evaluating...'
    out = densEstimator.eval(kde_qMatrix[::cheapSlice])
    print 'evaluation complete'
    num.savetxt('densityEstimationOut.txt', out)
else:
    out = num.loadtxt('densityEstimationOut.txt')

x = w[::cheapSlice]
y = d[::cheapSlice]
z = p[::cheapSlice]

scatterDens(x, y, z, out)
contourDens(x, y, z, out)

pylab.show()