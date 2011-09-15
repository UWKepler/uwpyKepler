from matplotlib import pyplot as plt
import uwpyKepler as kep
import pylab
import numpy as num
import sys

KeplerID = sys.argv[1]
#KeplerID = '1722276'
#KeplerID = '10341831'
#KeplerID = '6850504'

lcData1 = kep.iodb.ReadLightCurve(KeplerID,selection='LC')
eData = kep.iodb.getEclipseData(KeplerID)
BJDREFI = kep.iodb.getBJDREFI(KeplerID)
lcData = kep.pipeline.FlagKeplerEvents(lcData1)
lcData = kep.pipeline.RemoveBadEvents(lcData)
lcData = kep.pipeline.FlagEclipses(lcData,eData,BJDREFI)
lcDataP = kep.pipeline.SplitPortions(lcData,2)
lcDataP = kep.pipeline.FlagOutliers(lcDataP,10,5)
lcDataP = kep.pipeline.DetrendData(lcDataP, 50, 3)
lcData = kep.pipeline.StackPortions(lcDataP)

#fluxdt=[]
#yflux=[]
#xflux=[]

#for i in range(len(lcData1['y'])-1):
    #fluxdt.append(abs(lcData1['y'][i+1]-lcData1['y'][i]))
    #if fluxdt[i]>60:
        #idx.append(lcData1['y'][i])

ax1 = plt.subplot(2,1,1)
ax1.plot(lcData1['x'], lcData1['y'], 'bo')
ax2 = plt.subplot(2,1,2, sharex=ax1)
ax2.plot(lcData['x'],lcData['ydt'], 'bo')


pylab.show()