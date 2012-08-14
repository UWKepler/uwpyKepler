import uwpyKepler as kep
import numpy as num
import pylab
import sys

def whichPeriod(pDict):
    count = 1
    choiceList = [-1.]
    print 'choose a period by choice number'
    for key in pDict.keys():
        if len(pDict[key]) > 0:
            print 'from ' + key + ':'
            for period in pDict[key]:
                print '\t' + str(count) + ' - ' + str(period)
                choiceList.append(period)
                count += 1
    choice = int(raw_input('> '))

    return choiceList[choice]

def getPeriodDict(lc):
    ps = {'database':[], 'eDataDiscovery file':[], \
        'transitFit file': [], 'qats best period': []}
    eData1 = lc.eData
    if eData1['eDataExists']:
        for key in eData1['KOI'].keys():
            ps['database'].append(eData1['KOI'][key]['Period'])
    eData2 = kep.analysis.alleDataFromFile(lc.KID)
    for period in eData2[0]:
        ps['eDataDiscovery file'].append(period)
    eData3 = kep.analysis.geteDataFromModelFile(lc.KID)
    for period in eData3[0]:
        ps['transitFit file'].append(period)
    try:
        eData4 = kep.postqats.getBestPeriodByKID(lc.KID)
        ps['qats best period'].append(eData4)
    except:
        pass

    return ps

def plotT0s(lc, nhat):
    preIdx = num.array(map(int, list(nhat)))
    subIdx = num.where(preIdx < len(lc.lcFinal['x']))[0]
    idx = preIdx[subIdx]
    pylab.plot(lc.lcFinal['x'], lc.lcFinal['ydt'], 'b.')
    pylab.plot(lc.lcFinal['x'][idx], lc.lcFinal['ydt'][idx], 'ro')


kid = sys.argv[1]
lc = kep.quicklc.createLC(kid)
#lc = kep.keplc.keplc(kid)
periodDict = getPeriodDict(lc)
period = whichPeriod(periodDict)
qlc = kep.nhatQatsLC.nhatQatsLC(lc.lcFinal, lc.KID, period)
qlc.padLC()
qlc.runNQATS()
nhat = qlc.nhat
plotT0s(lc, nhat)

pylab.show()



