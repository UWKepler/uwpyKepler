import uwpyKepler as kep
from uwpyKepler.analysis import *
import numpy as num
import pylab
import scipy
import sys
import os
import optparse

def normalRun(mod, unfold):
    if lc.eData['eDataExists']:
        mod.fitIncArsRprs()
        mod.fitPeriodT0()
        mod.fitu1()
        mod.fitu2()
        mod.finalFit()
    else:
        mod.fitPeriodT0() # get decent period and t0 first
        mod.fitIncArsRprs()
        mod.fitPeriodT0()
        mod.fitu1()
        mod.fitu2()
        mod.finalFit()
    
    pylab.title('%s: ' % mod.kid + 'transit fit')
    pylab.ylabel('corrflux')
    if not unfold:
        pylab.xlabel('phase')
        pylab.plot(mod.phase, mod.ydata, 'b.')
        pylab.plot(mod.phase, mod.model, 'k.')
    else:
        pylab.xlabel('time in days')
        pylab.plot(mod.xdata, mod.ydata, 'b.')
        pylab.plot(mod.xdata, mod.model, 'k.')

def allPlots(mod, unfold):
    plots = ['fit to inc, aRs, RpRs', \
        'fit to period, t0', \
        'fit to u1, u2', \
        'final fit']
    funcs = [mod.fitIncArsRprs, \
        mod.fitPeriodT0, \
        mod.fitu1, \
        mod.fitu2, \
        mod.finalFit]
    if lc.eData['eDataExists']:
        for i in range(len(plots)):
            if i != 3:
                funcs[i]()
            if i == 2:
                funcs[i + 1]()
            pylab.figure()
            pylab.title('%s: ' % mod.kid + plots[i])
            pylab.ylabel('corrflux')
            if not unfold:
                pylab.xlabel('phase')
                pylab.plot(mod.phase, mod.ydata, 'b.')
                pylab.plot(mod.phase, mod.model, 'k.')
            else:
                pylab.xlabel('time in days')
                pylab.plot(mod.xdata, mod.ydata, 'b.')
                pylab.plot(mod.xdata, mod.model, 'k.')
    else:
        plots.insert(0, 'prelim period and t0 fit')
        funcs.insert(0, mod.fitPeriodT0)
        for i in range(len(plots)):
            if i != 4:
                funcs[i]()
            if i == 3:
                funcs[i + 1]()
            pylab.figure()
            pylab.title('%s: ' % mod.kid + plots[i])
            pylab.ylabel('corrflux')
            if not unfold:
                pylab.xlabel('phase')
                pylab.plot(mod.phase, mod.ydata, 'b.')
                pylab.plot(mod.phase, mod.model, 'k.')
            else:
                pylab.xlabel('time in days')
                pylab.plot(mod.xdata, mod.ydata, 'b.')
                pylab.plot(mod.xdata, mod.model, 'k.')

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

# insert fake eData
if not lc.eData['eDataExists']:
    if eDataInFile(kid):
        period, t0, q = geteDataFromFile(kid)
        lc.eData['eDataExists'] = True
        lc.eData['KOI'] = {'fake':\
                          {'Period':period, 'T0':t0, 'Duration':q}}
    else:
        # crude estimates to be determined by fit
        period = kep.postqats.getBestPeriodByKID(kid)
        t0 = 60.
        q = 1.5
        lc.eData['KOI'] = {'fake':\
                          {'Period':period, 'T0':t0, 'Duration':q}}

kw = kep.quicklc.quickKW(ctype=opts.ctype)
lc.runPipeline(kw)
RpRsEst = num.sqrt(1. - min(lc.lcFinal['ydt']))
firstGuess = num.array([num.pi / 2, 5., RpRsEst, 0.1, 0.1])

mod = modelLC(kid, lc, firstGuess)

if opts.all:
    allPlots(mod, opts.unfold)
else:
    normalRun(mod, opts.unfold)

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

