
import uwpyKepler as kep
import pylab
import numpy as num
import sys

#KeplerID = '1722276'
KeplerID = '10341831'
#KeplerID = '6850504'

lcData = kep.iodb.ReadLightCurve(KeplerID,selection='LC')
eData = kep.iodb.getEclipseData(KeplerID)
BJDREFI = kep.iodb.getBJDREFI(KeplerID)
lcData = kep.pipeline.FlagKeplerEvents(lcData)
lcData = kep.pipeline.RemoveBadEvents(lcData)
lcData = kep.pipeline.FlagEclipses(lcData,eData,BJDREFI)
lcDataP = kep.pipeline.SplitPortions(lcData,2)
lcDataP = kep.pipeline.FlagOutliers(lcDataP,10,5)
lcDataP = kep.pipeline.DetrendData(lcDataP,50,4)

for portion in lcDataP.keys():
    lcDataP = kep.lcmod.ApplyMaskPortions(lcDataP,'OKMask',portion)
    x = lcDataP[portion]['x']
    y = lcDataP[portion]['y']
    yerr = lcDataP[portion]['yerr']
    ycor = lcDataP[portion]['Correction']
    pylab.plot(x,y,'yo')
    #pylab.errorbar(x,y,yerr=yerr,ecolor=None)
    pylab.plot(x,ycor,'k.',linewidth=3)
    
pylab.show()