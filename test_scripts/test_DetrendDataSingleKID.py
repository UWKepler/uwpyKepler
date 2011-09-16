import uwpyKepler as kep
import pylab
import numpy as num
import sys

KeplerID = sys.argv[1]
#KeplerID = '1722276'
#KeplerID = '10341831'
#KeplerID = '6850504'

lcData = kep.iodb.ReadLightCurve(KeplerID,selection='LC')
eData = kep.iodb.getEclipseData(KeplerID)
BJDREFI = kep.iodb.getBJDREFI(KeplerID)
lcData = kep.pipeline.FlagKeplerEvents(lcData)
lcData = kep.pipeline.RemoveBadEvents(lcData)
lcData = kep.pipeline.FlagEclipses(lcData,eData,BJDREFI)
lcDataP = kep.pipeline.SplitPortions(lcData,2)
lcDataP = kep.pipeline.FlagOutliers(lcDataP,10,5)
lcDataP = kep.pipeline.DetrendData(lcDataP, 150, 5)

ic = 0
cl = ['b','g','r','y','c']
for portion in lcDataP.keys():
    lcDataP = kep.lcmod.ApplyMaskPortions(lcDataP,'ALLMask',portion)
    pylab.plot(lcDataP[portion]['x'],lcDataP[portion]['y'],cl[ic]+'o')
    pylab.plot(lcDataP[portion]['x'],lcDataP[portion]['correction'],'k-')
    ic += 1
    if ic > 4:
        ic = 0
pylab.show()

lcData = kep.pipeline.StackPortions(lcDataP)

pylab.plot(lcData['x'],lcData['ydt'],'bo')
pylab.show()

