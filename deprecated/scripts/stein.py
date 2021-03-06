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

eData=kep.io.getEclipseData(d1)

d2 = kep.io.FlagTransits(d1,eData)

pd = kep.io.SplitGap(d2,.1)

d3 = kep.io.FlagOutliers(pd,10,4)

#for portion in d3.keys():
    #print num.where(d3[portion]['OutlierMask']==True)

d4 = kep.proc.detrendData(d3,100,7)

print len(d4.keys()), ' portions'

#array = []
#arraydt =[]
#a3=[]
#for i in range(len(d1['x'])-1):
	#array.append(d1['y'][i+1]-d1['y'][i])
	#arraydt.append(d1['x'][i+1]-d1['x'][i])
        
#array= num.array(array)
#arraydt=num.array(arraydt)
#for i in range(len(array)-1):
    #a3.append(array[i+1]-array[i])
    
#a3=num.array(a3)

#pylab.plot(d1['x'][:-2],abs(a3),'b.')
#pylab.show()

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
    Mask = 'OTMask'
    print num.ma.count_masked(d4[portion]['x'])
    
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

dTransit = kep.io.ApplyMask(d3, 'UnMasked')

d10 = kep.proc.stackPortions(dTransit)

d11 = kep.proc.cutTransits(d10)

pylab.plot(d10['x'], d10['y'], 'b.')
pylab.plot(d11['x'],d11['y'],'r.')
pylab.show()
