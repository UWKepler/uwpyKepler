


import sys
import MySQLdb
import numpy as num
import pylab
import uwpyKepler as kep

keplerId = sys.argv[1]

#letting user know that it is using the given input
print 'Using', keplerId

d1 = kep.io.ReadLightCurve(keplerId)

pylab.plot(d1['x'],d1['y'],'ro')



d5=kep.io.getEclipseData(d1)



d6=kep.io.FlagTransits(d1,d5)
#pylab.plot(d6['x'],d6['y'],'y.')



d7=kep.io.SplitGap(d6,.1)



d8=kep.io.FlagOutliers(d7,10,4)



d81=kep.io.ApplyMask(d8,'UnMasked')

#for portion in d8.keys():
    #print num.where(d8[portion]['OutlierMask']==True)
    #print len(d8[portion]['OutlierMask'])


d9=kep.proc.stackPortions(d81)

#print num.where(d9['OutlierMask']==True)

#print d9.keys()

#d2=kep.proc.cutTransits(d1)

#pylab.plot(d2['x'],d2['y'],'y.')

#d10=kep.proc.cutOutliers(d9)

d11=kep.proc.onlyOutliers(d9)

#d12=kep.proc.onlyTransits(d9)


pylab.plot(d11['x'],d11['y'],'b.')

#d4 = kep.proc.cutAll(d1)
#pylab.plot(d4['x'],d4['y'],'b.')

pylab.show()
