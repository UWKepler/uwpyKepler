import sys
import uwpyKepler as kep
import numpy as num
import pylab

#kid = '10748390'
kid = '6056992'
d1 = kep.io.ReadLightCurve(kid)
eData = kep.io.getEclipseData(d1)
d2 = kep.io.FlagTransits(d1,eData)
pd = kep.io.SplitGap(d2,0.1)
d3 = kep.io.FlagOutliers(pd,10,4)
d4 = kep.proc.detrendData(d3,120,3)

pylab.plot(d1['x'],d1['y'],'b.')

for portion in d4.keys():
    #pylab.plot(d4[portion]['x'],d4[portion]['y'],'b.')
    pylab.plot(d4[portion]['x'], d4[portion]['Correction'], 'k-',linewidth=3)

pylab.title('Object ' + kid)
pylab.xlabel('Time in BJD')
pylab.ylabel('Flux')
pylab.legend(('Data', 'Correction Function'),'upper right', shadow=True)
pylab.show()


#pylab.plot(d4[portion]['x'],d4[portion]['y'],'y.')

d7 = kep.proc.stackPortions(d4)
#d8 = kep.proc.cutOT(d7)
#pylab.plot(d8['x'],d8['y'],'b.')

pylab.plot(d7['x'],d7['y'],'y.')

d8 = kep.proc.onlyOutliers(d7)

pylab.plot(d8['x'],d8['y'],'b.')

d9 = kep.proc.stackPortions(d4)
d10 = kep.proc.onlyTransits(d9)

pylab.plot(d10['x'],d10['y'],'r.')

pylab.title('Object ' + kid)
pylab.xlabel('Time in BJD')
pylab.ylabel('Flux')
pylab.legend(('Data', 'Outliers','Transits'),'lower right', shadow=True)

pylab.show()
