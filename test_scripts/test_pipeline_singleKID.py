
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
lcDataP = kep.pipeline.DetrendData(lcDataP, 70, 5)
lcData = kep.pipeline.StackPortions(lcDataP)

NeMask = kep.func.getNumBool(lcData['eMask'],True)
NKEMask = kep.func.getNumBool(lcData['KEMask'],True)
NOKMask = kep.func.getNumBool(lcData['OKMask'],True)
NOMask = kep.func.getNumBool(lcData['OMask'],True)
NALLMask = kep.func.getNumBool(lcData['ALLMask'],True)

rlist = ['all','elc','eonly','o','k','ok','flat','allflags']
Nplots = len(rlist)
i = 1
for typeName in rlist:
    if i == 1:
        ax1 = pylab.subplot(Nplots,1,i)
    else:
        pylab.subplot(Nplots,1,i, sharex=ax1)
    x,y,yerr = kep.lcmod.returnData(lcData,typeName)
    pylab.plot(x,y,'bo')
    #pylab.errorbar(x,y,yerr=yerr,fmt=None)
    pylab.title(typeName)
    i += 1
    
pylab.show()