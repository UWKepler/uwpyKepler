import sys
import uwpyKepler as kep
import numpy as num
import pylab

kid = '6056992'
#kid = '10341831'
#kid = '6050504'

d1 = kep.io.ReadLightCurve(kid)
eData = kep.io.getEclipseData(d1)
d2 = kep.io.FlagTransits(d1,eData)
pd = kep.io.SplitGap(d2,0.1)
d3 = kep.io.FlagOutliers(pd,10,4)
