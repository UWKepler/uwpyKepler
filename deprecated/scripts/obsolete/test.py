import sys
import MySQLdb
import numpy as num
import pylab
import uwpyKepler as kep

keplerId = sys.argv[1]

d1 = kep.io.ReadLightCurve(keplerId)
eData = kep.io.getEclipseData(d1)
d2 = kep.io.FlagTransits(d1,eData)


pd = kep.io.SplitGap(d2,0.05,2,0.05)

d3 = kep.io.FlagOutliers(pd,10,4)

d4 = kep.proc.detrendData(d3,100,5)
pylab.show()
d5 = kep.proc.stackPortions(d4)

print d5.keys()

d6 = kep.proc.cutTransits(d5)
#d7 = kep.proc.cutOutliers(d5)

#for portion in d4.keys():
	#pylab.plot(d4[portion]['x'],d4[portion]['y'],'bo')

pylab.plot(d6['x'],d6['y'],'r.')
#pylab.plot(d7['x'],d7['y'],'y.')
pylab.title(keplerId)
pylab.show()