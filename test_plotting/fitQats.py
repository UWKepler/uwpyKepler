import uwpyKepler as kep
import sys
import numpy as num
import pylab
import scipy
import os
import pdb
import time

def firstFareyValMod(periods, p0):
    Qs = []
    for i in range(len(periods) - 1):
        m, n = kep.postqats.firstFareyVal(p0, periods[i], periods[i + 1])
        Qs.append(1. / num.sqrt(m * n))
    Qs.append(Qs[-1])
    return num.array(Qs)

kid = sys.argv[1]
sgStr = str(kep.dbinfo.getSkyGroup(kid)).zfill(3)
dFileName = '/astro/store/student-scratch1/'\
            'johnm26/SPRING_BREAK_RUNS/SG' + sgStr + \
            '/signal.SG' + sgStr + \
            '.unflipped.' + kid + '.data'
dFile = open(dFileName, 'r')

order = 3
periods, snr, snrLC, snrFlat = kep.postqats.getQatsData(dFile)
periods = num.array(periods)
snr = num.array(snr)

allCoeffs = []
chiSqrs = []
fits = []
for p0 in periods:
    t0 = time.time()
    findQs = lambda Ps: firstFareyValMod(Ps, p0)
    funcs = map(lambda n: lambda x: x**n, range(order + 1)[::-1])
    funcs.append(findQs)
    #funcs = [3, 2, 1, 0, findQs]
    t0 = time.time()
    coeffs, fareyVals = \
    kep.linLeastSquares.linLeastSquares\
    (periods, snr, funcs, order + 2)
    allCoeffs.append(coeffs)
    #print 'solver time: ', time.time() - t0
    polynom = scipy.polyval(coeffs[:-1], periods)
    t0 = time.time()
    fit = polynom + coeffs[-1] * fareyVals
    #print 'findQs time: ', time.time() - t0
    fits.append(fit)
    sqrs = (snr - fit)**2 / fit
    chiSqr = sqrs.sum()
    chiSqrs.append(chiSqr)
    #print 'finished period: ', p0, time.time() - t0
minSqr = min(chiSqrs)
minIdx = chiSqrs.index(minSqr)

bestPeriod = periods[minIdx]
print 'best period from fit: ', bestPeriod
qatsBestPeriod = kep.postqats.getBestPeriodByKID(kid)
print 'best period from qats: ', qatsBestPeriod

bestIdx = periods.tolist().index(qatsBestPeriod)
bestFit = fits[bestIdx]

print 'chiSqr[minIdx]: ' + str(chiSqrs[minIdx])
print 'chiSqr[bestIdx]: ' + str(chiSqrs[bestIdx])
print 'coefficients: ' + str(allCoeffs[minIdx])

pylab.semilogx()
#a = pylab.subplot(211)
pylab.plot(periods[minIdx], snr[minIdx], 'ro')
pylab.plot(periods, fits[minIdx], 'c-')
#pylab.plot(periods, bestFit, 'r-')
#pylab.plot(periods, allQs[bestIdx], 'b-')

#pylab.plot(periods, allMs[bestIdx], 'b.')
#pylab.plot(periods, allNs[bestIdx], 'r.')
pylab.plot(periods[minIdx], fits[minIdx][minIdx], 'ro')
pylab.plot(periods, snr, 'k-')
#b = pylab.subplot(212, sharex=a)
#pylab.plot(periods, fits[chiSqrs.index(max(chiSqrs))])
#pylab.plot(minPeriods, minSnr, 'bo')
#pylab.plot(periods, diff, 'ro')

pylab.show()