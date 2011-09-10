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
#d12 = kep.proc.stackPortions(d3)
#d13 = kep.proc.cutTransits(d12)
#pylab.plot(d13['x'],d13['y'],'b.')
#pylab.show()
d4 = kep.proc.detrendData(d3,50,3)
#pylab.show()
d5 = kep.proc.stackPortions(d4)
if eData['bool']==False:
    bind = kep.proc.bin2(d5,[1, 5, 10, 20, 50, 100])
    
    i=0
    colors = ['r','y','b','c','g','m','k']
    
    for binsize in bind.keys():
            pylab.plot(bind[binsize]['x'], bind[binsize]['y'],colors[i]+'o')
            if i < 6:
                    i += 1
            else:
                    i = 0
            
    pylab.show()
else:
    d6 = kep.proc.cutTransits(d5)
    #d7 = kep.proc.cutTransits(d5)
    bind = kep.proc.bin2(d6,[1, 5, 10, 20, 50, 100])
    
    i=0
    colors = ['r','y','b','c','g','m','k']
    
    for binsize in bind.keys():
            pylab.plot(bind[binsize]['x'], bind[binsize]['y'],colors[i]+'o')
            if i < 6:
                    i += 1
            else:
                    i = 0
            
    pylab.show()

#pylab.plot(bind['binsize1']['x'],bind['binsize1']['y'],'r.')
#pylab.plot(bind['binsize5']['x'],bind['binsize5']['y'],'g.')
#pylab.plot(bind['binsize10']['x'],bind['binsize10']['y'],'y.')	
#pylab.plot(bind['binsize15']['x'],bind['binsize15']['y'],'m.')
#pylab.show()