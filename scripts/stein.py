import sys
import MySQLdb
import numpy as num
import pylab
import uwpyKepler as kep

keplerId = sys.argv[1]

#letting user know that it is using the given input
print 'Using', keplerId

d1 = kep.io.ReadLightCurve(keplerId)
print d1.keys(), len(d1['x'])

d2 = kep.io.FlagTransits(d1)
#print d2.keys()

pd = kep.io.SplitGap(d2,.1)

d3 = kep.io.FlagOutliers(pd,8,3)
d4 = kep.proc.detrendData(d3,75,5)

print len(d4.keys()), ' portions'

for portion in d4.keys():
    Mask = 'UnMasked'
    d4 = kep.io.ApplyMask(d4,Mask)
    #print num.ma.count_masked(d4[portion]['x'])
    pylab.plot(d4[portion]['x'],d4[portion]['y'],'b.')
    Mask = 'OutlierMask'
    d4 = kep.io.ApplyMask(d4,Mask)
    #print num.ma.count_masked(d4[portion]['x'])
    pylab.plot(d4[portion]['x'],d4[portion]['y'],'ro')
    Mask = 'TransitMask'
    d4 = kep.io.ApplyMask(d4,Mask)
    #print num.ma.count_masked(d4[portion]['x'])
    pylab.plot(d4[portion]['x'],d4[portion]['y'],'y.')
    Mask = 'OTMask'
    #print num.ma.count_masked(d4[portion]['x'])


pylab.show()