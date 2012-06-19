#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import sys
import numpy as num
import pylab
import scipy
import os

def firstFareyValMod(periods, p0):
    Qs = []
    for i in range(len(periods) - 1):
        m, n = kep.postqats.firstFareyVal(p0, periods[i], periods[i + 1])
        Qs.append(1. / num.sqrt(m * n))
    Qs.append(Qs[-1])
    return num.array(Qs)

kid = sys.argv[1]
dFile = open(kep.postqats.getdFileName(kid), 'r')

order = 3
periods, snr, snrLC, snrFlat = kep.postqats.getQatsData(dFile)
periods = num.array(periods)
snr = num.array(snr)

allCoeffs = []
chiSqrs = []
fits = []
for p0 in periods:
    findQs = lambda Ps: firstFareyValMod(Ps, p0)
    funcs = map(lambda n: lambda x: x**n, range(order + 1)[::-1])
    funcs.append(findQs)
    coeffs, fareyVals = \
    kep.linLeastSquares.linLeastSquares\
    (periods, snr, funcs, order + 2)
    allCoeffs.append(coeffs)
    polynom = scipy.polyval(coeffs[:-1], periods)
    fit = polynom + coeffs[-1] * fareyVals
    fits.append(fit)
    sqrs = (snr - fit)**2 / fit
    chiSqr = sqrs.sum()
    chiSqrs.append(chiSqr)
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

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.semilogx()
#a = pylab.subplot(211)
ax.plot(periods, snr, 'k-',\
        periods, fits[minIdx], 'c-',\
        periods[minIdx], snr[minIdx], 'ro',\
        periods[minIdx], fits[minIdx][minIdx], 'ro')
ax.legend(('qats SNR', 'fit to qats', 'fit best period'))
#pylab.plot(periods, bestFit, 'r-')
#pylab.plot(periods, allQs[bestIdx], 'b-')

#pylab.plot(periods, allMs[bestIdx], 'b.')
#pylab.plot(periods, allNs[bestIdx], 'r.')
#b = pylab.subplot(212, sharex=a)
#pylab.plot(periods, fits[chiSqrs.index(max(chiSqrs))])
#pylab.plot(minPeriods, minSnr, 'bo')
#pylab.plot(periods, diff, 'ro')

pylab.ylabel('SNR')
pylab.xlabel('periods')

pylab.show()