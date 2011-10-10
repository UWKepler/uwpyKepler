
import uwpyKepler as kep
import pylab
import numpy as num
import sys, os
import traceback
import optparse

parser = optparse.OptionParser()
parser.add_option('-l','--lcdtwin',type=int,default=100,help='Detrending window size for LC data.')
parser.add_option('-s','--scdtwin',type=int,default=10000,help='Detrending window size for SC data.')
parser.add_option('-o','--order',type=int,default = 7,help='Detrending function polynomial order.')
parser.add_option('-m','--lcoutwin',type=int,default=10,help='Window size for determining median of LC data points.')
parser.add_option('-n','--scoutwin',type=int,default=1000,help='Window size for determining median of SC data points.')
parser.add_option('-t','--threshold',type=int,default=3,help="Sigma threshold for flagging outliers.")
parser.parse_args()

SGNumber = sys.argv[0]
allkids = kep.dbinfo.returnSkyGroupsKIDs(SGNumber)
selectoption = ['LC','SC','']
outlierSigma = 3
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
                if selopt =='LC':
                    outwin = 10
                    dtwin = 100
                else:
                    outwin = 1000
                    dtwin = 10000
                lcDataP = kep.pipeline.FlagOutliers(lcDataP,outwin,sig)
                lcDataP = kep.pipeline.DetrendData(lcDataP,dtwin,7)
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