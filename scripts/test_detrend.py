
import uwpyKepler as kep
import numpy as num

kid = '6056992'
#kid = '5812701'
d1 = kep.io.ReadLightCurve(kid)
eData = kep.io.getEclipseData(d1)
d2 = kep.io.FlagTransits(d1,eData)
pd = kep.io.SplitGap(d2,0.1)
d3 = kep.io.FlagOutliers(pd,10,4)
d4 = kep.proc.detrendData(d3,100,7)

import pylab
for portion in d4.keys():
    print num.shape(num.where(d3[portion]['OTMask'])), num.shape(num.where(d3[portion]['TransitMask'])),num.shape(num.where(d3[portion]['OutlierMask'])),num.shape(num.where(d3[portion]['UnMasked']))
    pylab.plot(d4[portion]['x'],d4[portion]['y'],'b.')
    pylab.plot(d4[portion]['x'], d4[portion]['Correction'], 'k-')
    
pylab.show()