import sys
import uwpyKepler as kep
import numpy as num
import pylab

kid = sys.argv[1]

d1 = kep.io.ReadLightCurve(kid)
eData = kep.io.getEclipseData(d1)
d2 = kep.io.FlagTransits(d1,eData)
pd = kep.io.SplitGap(d2,0.05,2,0.05)
d3 = kep.io.FlagOutliers(pd,10,4)
d4 = kep.proc.detrendData(d3,50,3)
d5 = kep.proc.stackPortions(d4)
d6 = kep.proc.cutOutliers(d5)
d7 = kep.proc.bin(d6)