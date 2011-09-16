
import uwpyKepler as kep
import pylab
import numpy as num
import sys, os
import traceback

selectoption = ['LC','SC','']
outlierSigma = [3]
#allKIDs = kep.dbinfo.getKIDsSource()
allKIDs = ['9631995','8845026','3544595']
failcount = 0
passcount = 0
totalcount = len(allKIDs)*len(selectoption)*len(outlierSigma)

os.system('rm -v allKIDs.pipeline.log')
print 'Total Runs = ',totalcount
for KeplerID in allKIDs:
    for selopt in selectoption:
        for sig in outlierSigma:
            try:
                lcData = kep.iodb.ReadLightCurve(KeplerID,selection=selopt)
                eData = kep.iodb.getEclipseData(KeplerID)
                BJDREFI = kep.iodb.getBJDREFI(KeplerID)
                lcData = kep.pipeline.FlagKeplerEvents(lcData)
                lcData = kep.pipeline.RemoveBadEvents(lcData)
                lcData = kep.pipeline.FlagEclipses(lcData,eData,BJDREFI)
                lcDataP = kep.pipeline.SplitPortions(lcData,2)
                lcDataP = kep.pipeline.FlagOutliers(lcDataP,10,sig)
                lcDataP = kep.pipeline.DetrendData(lcDataP,100,7)
                lcData = kep.pipeline.StackPortions(lcDataP)
                passcount += 1
                print 'pass ',KeplerID, selopt, sig
            except:
                failcount += 1
                print 'fail ',KeplerID, selopt, sig
                ofile = open("allKIDs.pipeline.log","a")
                traceback.print_exc(file=open("allKIDs.pipeline.log","a"))
                print >> ofile,'# ',KeplerID, selopt, sig
                continue

print totalcount, passcount, failcount 