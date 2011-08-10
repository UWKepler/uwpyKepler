import sys
import MySQLdb
import numpy as num
import pylab
import uwpyKepler as kep

keplerID = sys.argv[1]
keplerID = '9775938'
keplerID = '10341831'

d1 = kep.io.ReadLightCurve(keplerID)
pylab.plot(d1['x'],d1['y'],'ro')
dread=kep.io.getEclipseData(d1)
d2 = kep.io.FlagTransits(d1,dread)
pylab.plot(d2['x'],d2['y'],'b.')

pd = kep.io.SplitGap(d2,.1, 10, 0.05)


