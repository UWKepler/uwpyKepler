import MySQLdb
import sys
import numpy as num
import scipy
import pylab
import uwpyKepler as kep

keplerid = sys.argv[1]

d1 = kep.io.ReadLightCurve(keplerid)
pylab.plot(d1['x'],d1['y'],'ro')
eData = kep.io.getEclipseData(d1)
d2 = kep.io.FlagTransits(d1,eData)
pylab.plot(d2['x'],d2['y'],'b.')
d3 = kep.io.SplitGap(d2,0.05,2,0.05)

d4 = kep.io.FlagOutliers(d3,10,4)

d5 = kep.proc.detrendData(d4,100,7)

pylab.show()