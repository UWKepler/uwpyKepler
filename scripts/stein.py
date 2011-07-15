import sys
import MySQLdb
import numpy as num
import pylab
import uwpyKepler as kep

keplerId = sys.argv[1]
#keplerId = '11853905'
#keplerId = '10666592'
#keplerId = '7907423'


#letting user know that it is using the given input
print 'Using', keplerId

d1 = kep.io.ReadLightCurve(keplerId)
print d1.keys(), len(d1['x'])

d2 = kep.io.FlagTransits(d1)

#print d2.keys()

pd = kep.io.SplitGap(d2,.1)

d3 = kep.io.FlagOutliers(pd,10,4)
#print d3.keys()
#print len(num.where(d3['portion1']['UnMasked']==False)[0])
#print len(num.where(d3['portion1']['TransitMask']==False)[0])
#print len(num.where(d3['portion1']['OutlierMask']==False)[0])
#print len(num.where(d3['portion1']['MaskBoth']==False)[0])

#sys.exit()
#for portion in d3.keys():
    ##print len(d3[portion]['x']), len(num.ma.getdata(d3[portion]['x'])), portion
    ##print len(d3[portion]['x'].mask)
    #d3[portion]['x'].mask = d3[portion]['UnMasked']
    #d3[portion]['y'].mask = d3[portion]['UnMasked']
    #x1 = num.ma.copy(d3[portion]['x'])
    #y1 = num.ma.copy(d3[portion]['y'])
    ##print len(d3[portion]['x'].mask)
    #d3[portion]['x'].mask = d3[portion]['TransitMask']
    #d3[portion]['y'].mask = d3[portion]['TransitMask']
    #x2 = num.ma.copy(d3[portion]['x'])
    #y2 = num.ma.copy(d3[portion]['y'])
    ##print len(d3[portion]['x'].mask)
    #d3[portion]['x'].mask = d3[portion]['UnMasked']
    #d3[portion]['y'].mask = d3[portion]['UnMasked']
    #d3[portion]['x'].mask = d3[portion]['OutlierMask']
    #d3[portion]['y'].mask = d3[portion]['OutlierMask']
    #x3 = num.ma.copy(d3[portion]['x'])
    #y3 = num.ma.copy(d3[portion]['y'])
    ##print len(d3[portion]['x'].mask)
    #pylab.plot(x1,y1,'bo')
    #pylab.plot(x2,y2,'r.')
    #pylab.plot(x3,y3,'y.')
    ##pylab.plot(num.ma.getdata(d3[portion]['x']),num.ma.getdata(d3[portion]['y']),'r.')
    ##pylab.plot(d3[portion]['x'],d3[portion]['y'],'bo')
    #pylab.title(portion)
    #pylab.show()
    
d4 = kep.proc.detrendData(d3,300,7)


for portion in d4.keys():
    Mask = 'UnMasked'
    d4[portion]['x'].mask = d4[portion][Mask]
    d4[portion]['y'].mask = d4[portion][Mask] 
    print num.ma.count_masked(d4[portion]['x'])
    #print num.ma.count(d4[portion]['x']), portion, len(d4[portion]['x']), len(d4[portion]['x'].data)
    pylab.plot(d4[portion]['x'],d4[portion]['y'],'b.')
    Mask = 'OutlierMask'
    d4[portion]['x'].mask = d4[portion][Mask]
    d4[portion]['y'].mask = d4[portion][Mask]
    print num.ma.count_masked(d4[portion]['x'])
    pylab.plot(d4[portion]['x'],d4[portion]['y'],'ro')
    Mask = 'TransitMask'
    d4[portion]['x'].mask = d4[portion][Mask]
    d4[portion]['y'].mask = d4[portion][Mask]
    print num.ma.count_masked(d4[portion]['x'])
    pylab.plot(d4[portion]['x'],d4[portion]['y'],'y.')
    Mask = 'MaskBoth'
    d4[portion]['x'].mask = d4[portion][Mask]
    d4[portion]['y'].mask = d4[portion][Mask]
    #pylab.plot(d4[portion]['x'],d4[portion]['y'],'g.')
    print num.ma.count_masked(d4[portion]['x'])

    #print num.ma.count_masked(d4[portion]['OutlierMask']), num.ma.count_masked(d4[portion]['MaskBoth']), num.ma.count_masked(d4[portion]['TransitMask']), len(d4[portion]['x'].data), portion

pylab.show()

#dout = kep.proc.detrendData(d3, 25, 200, 5)

#for portion in dout.keys():
    #pylab.plot(dout[portion]['x'],dout[portion]['y'],'b.')
    #pylab.plot(num.ma.getdata(dout[portion]['x']),num.ma.getdata(dout[portion]['y']),'r.')

    
#pylab.show()



#x = []
#y = []
#yerr = []
#for portion in d3.keys():
    #for i in range(len(d3[portion]['OutlierFlag'])):
        #if not d3[portion]['OutlierFlag'][i]:
            #x.append(d3[portion]['x'][i])
            #y.append(d3[portion]['y'][i])
        ##yerr.append(d2['TransitFlag']['x'])
    #pylab.plot(x,y,'bo')
    #pylab.plot(d3[portion]['x'],d3[portion]['y'],'g.')
    #pylab.show()
        

#pylab.errorbar(d1['x'], d1['y'], fmt='ro')
#pylab.errorbar(d3['x'], d3['y'], fmt='go')
#pylab.ylim()



#d4 = kep.io.detrend2(d3,'portion7',5)

#pylab.plot(d3['portion7']['x'],d3['portion7']['y'],'ro')
#pylab.plot(d4['x'],d4['y'],'b.')
#pylab.show()


#pylab.show()