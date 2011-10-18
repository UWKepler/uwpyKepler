import uwpyKepler as kep
import pylab
import numpy as num
import sys, os
import traceback
import optparse

parser = optparse.OptionParser()
parser.add_option('-l','--lcdtwin',\
                   type=int,\
                   default=100,\
                   help='Detrending window size for LC data. Default=100')
parser.add_option('-s','--scdtwin',\
                   type=int,\
                   default=10000,\
                   help='Detrending window size for SC data. Default=10000')
parser.add_option('-o','--order',\
                   type=int,\
                   default=7,\
                   help='Detrending function polynomial order.\nDefault=7')
parser.add_option('-m','--lcoutwin',\
                   type=int,\
                   default=10,\
                   help='Window size for determining median of LC data points. Default=10')
parser.add_option('-n','--scoutwin',\
                   type=int,\
                   default=1000,\
                   help='Window size for determining median of SC data points. Default=1000')
parser.add_option('-t','--thresholdmax',\
                   type=int,\
                   default=3,\
                   help='Maximum sigma threshold for flagging outliers. Default=3')
parser.add_option('-b','--thresholdmin',\
                   type=int,\
                   default=3,\
                   help='Minimum sigma threshold for flagging outliers. Each KID will be tested from the minimum to the maximum threshold. Default=3')
options, SGNumber = parser.parse_args()

#allkids = kep.dbinfo.returnSkyGroupKIDs(SGNumber[0])
selectoption = ['LC']
#outlierSigma = [3]
outlierSigma = range(options.thresholdmin, options.thresholdmax+1)
print "Outlier Sigmas: ", outlierSigma
#allKIDs = kep.dbinfo.getKIDsSource()
#allKIDs = ['9631995','8845026','3544595']
allKIDs = ['9631995']
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
                    outwin = options.lcoutwin
                    dtwin = options.lcdtwin
                else:
                    outwin = options.scoutwin
                    dtwin = options.scdtwin
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

kep.qatsprep.padLC(lcData)
print totalcount, passcount, failcount 