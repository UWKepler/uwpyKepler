
import uwpyKepler as kep
import pylab
import numpy as num
import sys

#KeplerID = '1722276'
#KeplerID = '10341831'
KeplerID = '6850504'

lcData = kep.iodb.ReadLightCurve(KeplerID,selection='LC')
eData = kep.iodb.getEclipseData(KeplerID)
BJDREFI = kep.iodb.getBJDREFI(KeplerID)
lcData = kep.pipeline.FlagKeplerEvents(lcData)
lcData = kep.pipeline.RemoveBadEvents(lcData)
lcData = kep.pipeline.FlagEclipses(lcData,eData,BJDREFI)
lcDataP = kep.pipeline.SplitPortions(lcData,2)
lcDataP = kep.pipeline.FlagOutliers(lcDataP,10,5)
lcDataP = kep.pipeline.DetrendData(lcDataP, 50, 3)
lcData = kep.pipeline.StackPortions(lcDataP)

for typeName in ['all','elc','eonly','o','k','ok']:
    x,y,yerr = kep.lcmod.returnData(lcData,typeName)
    pylab.plot(x,y,'bo')
    pylab.title(typeName)
    pylab.show()