import sys
import MySQLdb
import numpy as num
import pylab
import uwpyKepler as kep

keplerId = sys.argv[1]

d1 = kep.io.ReadLightCurve(keplerId)

dread=kep.io.db(d1)

d2 = kep.io.FlagTransits(d1,dread)

pd = kep.io.SplitGap(d2,.1)

d3 = kep.io.FlagOutliers(pd,10,4)

d4 = kep.proc.detrendData(d3,100,5)