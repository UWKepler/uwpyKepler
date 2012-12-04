from getdbData import *
import numpy as num
import pylab
from mpl_toolkits.mplot3d import Axes3D

QATS_F = 0.005

def plotTopQats(kid):
    qMatrix = qatsResults_3d()
    snr = getTopSpectrum(qMatrix)
    periods = reconstructPeriodArray(snr, QATS_F)
    #pylab.plot(num.arange(len(snr)), snr)
    pylab.plot(periods, snr)

def qatsContours(kid):
    qMatrix = qatsResults_3d()
    ax = pylab.gca(projection='3d')
    
    widths = num.arange(qMatrix.shape[0])
    depths = num.arange(qMatrix.shape[1])
    
    #ax.contour(widths, depths, qMatrix)
    
    x = num.linspace(0, 2*num.pi, num=100)
    y = num.linspace(0, 2*num.pi, num=100)
    X, Y = num.meshgrid(x, y)
    Z = num.sin(X) * num.cos(Y)
    ax.plot_surface(X, Y, Z, linewidth=0, alpha=0.5)

# qMatrix contains too much data for practical plotting
def qatsScatter_3d(kid):
    qMatrix = qatsResults_3d()
    ax = pylab.gca(projection='3d')
    
    snr = getTopSpectrum(qMatrix)
    idx = num.where(snr != 0)[0]
    qMatrix = qMatrix[:, :, idx]
    
    widths  = num.arange(qMatrix.shape[0])
    depths  = num.arange(qMatrix.shape[1])
    periods = num.arange(qMatrix.shape[2])
    
    widths, periods, depths = meshgrid_3d(widths, periods, depths)
    ax.scatter(widths, periods, depths, c=qMatrix.ravel(), cmap=pylab.cm.RdBu)

def qatsAnimate(kid):
    qMatrix = qatsResults_3d()
    snr = getTopSpectrum(qMatrix)
    periods = reconstructPeriodArray(snr, QATS_F)
    periods = num.round(periods, decimals=2)
    pylab.ion()
    
    # sub1 for the width-depth map
    sub1 = pylab.subplot(211)
    # sub2 for the top spectrum slice plot
    sub2 = pylab.subplot(212)
    ax1 = sub1.get_axes()
    ax2 = sub2.get_axes()
    # set title and other boring stuff
    ax1.set_title('QATS width-depth map @ period: ' + \
        str(periods[0]))
    ax1.set_xlabel('transit widths')
    ax1.set_ylabel('transit depths')
    
    # set initial image
    im = ax1.imshow(qMatrix[:,:,0], cmap=pylab.cm.RdBu, \
        vmax=qMatrix.max(), vmin=qMatrix.min())
    # set static colorbar
    cbar = pylab.colorbar(im, ax=ax1)
    # plot top spectrum
    ax2.plot(periods, snr)
    ax2.set_xscale('log')
    # initialize slice line
    vertical = ax2.add_line(pylab.matplotlib.lines.Line2D( \
        [periods[0], periods[0]], \
        [-1, qMatrix.max() * (1.4)], linewidth=1, color='r'))
    
    # import time to add delay to frame rate
    import time
    for i in range(1, len(qMatrix[0,0,:])):
        # reset certain values at each iteration (i.e. each slice)
        ax1.set_title('QATS width-depth map @ period: ' + \
            str(periods[i]))
        im.set_array(qMatrix[:,:,i])
        vertical.set_xdata([periods[i]])
        pylab.draw()
        # time delay
        time.sleep(0.25)

# 3d version of numpy.meshgrid
def meshgrid_3d(x, y, z):
    xlen = len(x)
    ylen = len(y)
    zlen = len(z)
    x = x.reshape(1, 1, xlen)
    y = y.reshape(1, ylen, 1)
    z = z.reshape(zlen, 1, 1)
    X = (x.repeat(ylen, axis=1)).repeat(zlen, axis=0)
    Y = (y.repeat(xlen, axis=2)).repeat(zlen, axis=0)
    Z = (z.repeat(xlen, axis=2)).repeat(ylen, axis=1)

    return X, Y, Z