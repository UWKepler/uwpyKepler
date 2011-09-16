
import uwpyKepler as kep
import pylab
import numpy as num
import sys

#KeplerID = '1722276'
#KeplerID = '10341831'
#KeplerID = '6850504'
KeplerID = sys.argv[1]

lcData = kep.iodb.ReadLightCurve(KeplerID,selection='LC')
eData = kep.iodb.getEclipseData(KeplerID)
BJDREFI = kep.iodb.getBJDREFI(KeplerID)
lcData = kep.pipeline.FlagKeplerEvents(lcData)
lcData = kep.pipeline.RemoveBadEvents(lcData)
lcData = kep.pipeline.FlagEclipses(lcData,eData,BJDREFI)
lcDataP = kep.pipeline.SplitPortions(lcData,2)

ic = 0
cl = ['k','b','g','r','y','c']
for portion in lcDataP.keys():
    lcDataP = kep.lcmod.ApplyMaskPortions(lcDataP,'NoMask',portion)
    x = lcDataP[portion]['x']
    y = lcDataP[portion]['y']
    yerr = lcDataP[portion]['yerr']
    pylab.plot(x,y,cl[ic]+'o')
    ic += 1
    if ic > 5:
        ic = 0
    #pylab.errorbar(x,y,yerr=yerr,ecolor=None)
    
pylab.show()