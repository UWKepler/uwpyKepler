#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
from uwpyKepler.analysis import *
import numpy as num
import pylab
import sys
import optparse

# plots mod.ydata w/ or w/o errorbars
# determined by errBool
def withOrWithoutErr(x, errBool):
    if errBool:
        pylab.errorbar(x, mod.ydata, \
            yerr=mod.yerr, fmt='.', color='b')
    else:
        pylab.plot(x, mod.ydata, 'b.')

# plots x and y data and model with
# "unfold" and "error" user specifications
def foldOrUnfold(unfoldBool, errBool):
    if not unfoldBool:
        pylab.xlabel('phase')
        withOrWithoutErr(mod.phase, errBool)
        pylab.plot(mod.phase, mod.model, 'k.')
    else:
        pylab.xlabel('time in days')
        withOrWithoutErr(mod.xdata, errBool)
        pylab.plot(mod.xdata, mod.model, 'k.')

def normalRun(mod, unfold, error):
    if lc.eData['eDataExists']:
        mod.fitIncArsRprs()
        mod.fitPeriodT0()
        mod.fitu1()
        mod.fitu2()
        mod.finalFit()
    else:
        mod.fitT0() # lines up model box fit and transit.. ideally
        mod.fitPeriodT0() # get decent period and t0 first
        mod.fitIncArsRprs()
        mod.fitPeriodT0()
        mod.fitu1()
        mod.fitu2()
        mod.finalFit()
    
    pylab.title('%s: ' % mod.kid + 'transit fit')
    pylab.ylabel('corrflux')
    foldOrUnfold(unfold, error)

# plots x and y data and model at
# all stages of fit
# using "unfold" and "error" user specifications
def allPlots(mod, unfold, error):
    plots = ['fit to inc, aRs, RpRs', \
        'fit to period, t0', \
        'fit to u1, u2', \
        'final fit']
    # functions called before each plot
    funcs = [(mod.fitIncArsRprs,), \
        (mod.fitPeriodT0,), \
        (mod.fitu1, mod.fitu2), \
        (mod.finalFit,),]
    if lc.eData['eDataExists']:
        for i in range(len(plots)):
            for func in funcs[i]:
                func()
            pylab.figure()
            pylab.title('%s: ' % mod.kid + plots[i])
            pylab.ylabel('corrflux')
            foldOrUnfold(unfold, error)
    else:
        plots.insert(0, 'prelim period and t0 fit')
        funcs.insert(0, (mod.fitPeriodT0,))
        plots.insert(0, 'prelim t0 fit')
        funcs.insert(0, (mod.fitT0,))
        for i in range(len(plots)):
            for func in funcs[i]:
                func()
            pylab.figure()
            pylab.title('%s: ' % mod.kid + plots[i])
            pylab.ylabel('corrflux')
            foldOrUnfold(unfold, error)

# writes best fit parameters 
# to specified file overriding an entry with same KID
def write_NoDuplicates(file, mod):
    idx = None
    newline = str(mod.kid) + ' ' + \
        str(mod.period) + ' ' + \
        str(mod.t0) + ' ' + \
        str(mod.q) + ' ' + \
        str(mod.inc) + ' ' + \
        str(mod.aRs) + ' ' + \
        str(mod.RpRs) + ' ' + \
        str(mod.u1) + ' ' + \
        str(mod.u2)
    ofile = open(file, 'r')
    lines = ofile.readlines()
    for i in range(len(lines)):
        if lines[i].split()[0] == str(mod.kid):
            idx = i
            break
    if idx:
        lines[idx] = newline
    else:
        lines.append(newline)
    ofile.close()
    ofile = open(file, 'w')
    for line in lines:
        print >> ofile, line.strip()

def makeModelLC(lc):
    # estimate RpRs without picking up outliers
    ydt = lc.lcFinal['ydt']
    idx = ydt.argsort()
    iMin = idx[num.floor(0.1585*len(idx))]
    RpRsEst = num.sqrt(1. - ydt[iMin])
    
    # aRs is FINNICKY;
    # fit can change drastically based on guess!!
    # estimate aRs taking a = 1 / sin(q*pi) where q ~ 0.02
    aRsEst = 1. / num.sin(num.pi * 0.02)
    
    # estimate inc to be 90 degrees
    incEst = num.pi / 2
    
    u1Est = 0.1
    u2Est = 0.1
    
    firstGuess = num.array([incEst, aRsEst, RpRsEst, u1Est, u2Est]) 
    lcData = lc.lcFinal
    eData  = lc.eData
    return modelLC(kid, lcData, eData, firstGuess)

# returns ydata indices of all points NOT IN TRANSIT using eData 
# in catalogue and using periods, t0s, and qs from file
# except the last set of these encountered
def excludeKnownTransitIdx(lc, ctype):
    rawlcData = kep.iodb.ReadLightCurve(lc.KID, selection=ctype)
    rawlcData = kep.pipeline.FlagKeplerEvents(rawlcData)
    rawlcData = kep.pipeline.RemoveBadEvents(rawlcData)
    eMask = kep.pipeline.FlagEclipses( \
        rawlcData,lc.eData,lc.BJDREFI)['eMask']
    return num.where(eMask == 0)[0]

def insertEstimatedEData(lc, periods, t0s, qs, i):
    lc.eData['KOI']['fake' + str(i + 1)] = { \
        'Period':periods[i], \
        'T0':t0s[i], \
        'Duration':qs[i]}
    
# only used if no data found in file;
# crude estimates to be determined by fit
def insertCrudeEData(lc, fakeCount):
    print 'no eData found in file; using hardcode guesses'
    period = kep.postqats.getBestPeriodByKID(kid)
    t0 = 60.
    q = 1.5
    lc.eData['KOI']['fake' + str(fakeCount)] = { \
        'Period':period, \
        'T0':t0, \
        'Duration':q}
    
def singleMain(lc, ctype):
    if not lc.eData['eDataExists']:
        lc.eData['eDataExists'] = True
        period, t0, q = geteDataFromFile(kid)
        if period == -1:
            insertCrudeEData(lc, 0)
        else:
            insertEstimatedEData(lc, [period], [t0], [q], 0)
    kw = kep.quicklc.quickKW(ctype=ctype)
    lc.runPipeline(kw)

def multiMain(lc, nPlanets, ctype):
    useCrude = False
    periods, t0s, qs = alleDataFromFile(lc.KID)
    useCrude = len(periods) < 1
    if lc.eData['eDataExists']:
        koiCount = len(lc.eData['KOI'])
        if len(lc.eData['KOI']) > nPlanets or \
        len(lc.eData['KOI']) + len(periods) > nPlanets + 1:
            print 'recorded KOIs exceed expected planet count;'
            print 'exiting'
            exit()
    else:
        koiCount = 0
    # insert eData of known transits to be removed
    fakeCount = 0
    while koiCount < nPlanets:
        if not useCrude:
            insertEstimatedEData(lc, periods, t0s, qs, fakeCount)
            useCrude = fakeCount == len(periods)
        else:
            insertCrudeEData(lc, fakeCount)
        fakeCount += 1
        koiCount  += 1
    idx = excludeKnownTransitIdx(lc, ctype)
    # delete known transit eData
    for key in lc.eData['KOI'].keys():
        del lc.eData['KOI'][key]
    # insert preliminary eData of transit to be modeled
    if not useCrude:
        insertEstimatedEData(lc, periods, t0s, qs, fakeCount)
    else:
        insertCrudeEData(lc, fakeCount)
    # exclude known transits from lightcurve
    kw = kep.quicklc.quickKW(ctype=ctype)
    lc.runPipeline(kw)
    lc.lcFinal['x']      = lc.lcFinal['x'][idx]
    lc.lcFinal['ydt']    = lc.lcFinal['ydt'][idx]
    lc.lcFinal['yerrdt'] = lc.lcFinal['yerrdt'][idx]


if __name__ == '__main__':
    parser = optparse.OptionParser(usage=\
    "%prog\nFits lightcurve with Eric's transit model.\nScript takes a KID as a system argument.")
    parser.add_option('-a','--all',\
                        action='store_true',\
                        dest='all',\
                        default=False,\
                        help='display plots at all'\
                        'stages of fit')
    parser.add_option('-w','--write',\
                        type='string',\
                        dest='write',\
                        default=None,\
                        help='write best fit '\
                        'parameters to given file')
    parser.add_option('-u','--unfold',\
                        action='store_true',\
                        dest='unfold',\
                        default=False,\
                        help='display lightcurve unfolded')
    parser.add_option('-e','--error',\
                        action='store_true',\
                        dest='error',\
                        default=False,\
                        help='plot ydata with errors')
    parser.add_option('-m','--multiplanet',\
                        type=float,\
                        dest='nPlanets',\
                        default=None,\
                        help='input number known planets' \
                        'in system;\nmasks these when calculating' \
                        'fit')
    cChoice = ('LC','SC','')
    parser.add_option('-c','--ctype',\
                        choices=cChoice,\
                        type='choice',\
                        dest='ctype',\
                        default='LC',\
                        help='The expected cadence type('\
                        +', '.join(cChoice)\
                        +') [default: %default]')
    opts, args = parser.parse_args()

kid = sys.argv[1]
lc = kep.keplc.keplc(kid)

# exclude known planets from lightcurve
if opts.nPlanets:
    multiMain(lc, opts.nPlanets, opts.ctype)
else:
    singleMain(lc, opts.ctype)

mod = makeModelLC(lc)
if opts.all:
    allPlots(mod, opts.unfold, opts.error)
else:
    normalRun(mod, opts.unfold, opts.error)

if opts.write:
    write_NoDuplicates(opts.write, mod)

print 'BEST FIT PARAMETERS:'
print 'period:', mod.period
print 't0:', '\t', mod.t0
print 'inc:', '\t', mod.inc
print 'aRs:', '\t', mod.aRs
print 'RpRs:', '\t', mod.RpRs
print 'u1:', '\t', mod.u1
print 'u2:', '\t', mod.u2

pylab.show()

