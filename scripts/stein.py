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

d3 = kep.io.FlagOutliers(pd,10,4)
d4 = kep.proc.detrendData(d3,100,7)

print len(d4.keys()), ' portions'
#sys.exit()
for portion in d4.keys():
    Mask = 'UnMasked'
    d4 = kep.io.ApplyMask(d4,Mask)
    #print num.ma.count_masked(d4[portion]['x'])
    pylab.plot(d4[portion]['x'],d4[portion]['y'],'b.')
    
    Mask = 'OutlierMask'
    d4 = kep.io.ApplyMask(d4,Mask)
    #print num.ma.count_masked(d4[portion]['x'])
    pylab.plot(d4[portion]['x'],d4[portion]['y'],'r.')
    Mask = 'OTMask'
    d4 = kep.io.ApplyMask(d4,Mask)
    #print num.ma.count_masked(d4[portion]['x'])
    pylab.plot(d4[portion]['x'],d4[portion]['y'],'y.')
    pylab.title('Object ' + keplerId)
    pylab.xlabel('Time in BJD')
    pylab.ylabel('Flux')
    pylab.legend(('Outliers', 'Transits','Data'),'lower right', shadow=True)
    #Mask = 'OTMask'
    #print num.ma.count_masked(d4[portion]['x'])
    
pylab.show()

d3 = kep.io.ApplyMask(d3,'UnMasked')

for portion in d4.keys():
	pylab.plot(d3[portion]['x'],d3[portion]['y'],'b.')
	pylab.plot(d4[portion]['x'],d4[portion]['Correction'],'k-',linewidth=3)
	pylab.title('Object ' + keplerId)
	pylab.xlabel('Time in BJD')
	pylab.ylabel('Flux')
	pylab.legend(('Data', 'Correction Function'),
           'upper right', shadow=True)



pylab.show()