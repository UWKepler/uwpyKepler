#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import sys
import numpy as num
import pylab
import scipy
import os

kid = sys.argv[1]
dFile = open(kep.postqats.getdFileName(kid), 'r')
order = 3
periods, snr, snrLC, snrFlat = kep.postqats.getQatsData(dFile)
bestPeriod, bestFit, minSqr,\
qatsBestPeriod, qatsPeriodFit, qatsPeriodSqr\
= kep.postqats.fitQats(kid, periods, snr, order,\
                       qats_best_period_fit=True)
bestPeriodIdx = periods.tolist().index(bestPeriod)
qatsPeriodIdx = periods.tolist().index(qatsBestPeriod)

print 'best period from fit: ', bestPeriod
print 'best period from qats: ', qatsBestPeriod
print 'best fit chi sqr: ' + str(minSqr)
print 'best qats period chi sqr: ' + str(qatsPeriodSqr)

fig = pylab.figure()
ax = fig.add_subplot(111)
ax.semilogx()
ax.plot(periods, snr, 'k-',\
        periods, bestFit, 'c-',\
        periods[bestPeriodIdx], snr[bestPeriodIdx], 'ro',\
        periods[bestPeriodIdx], bestFit[bestPeriodIdx], 'ro',\
        periods[qatsPeriodIdx], snr[qatsPeriodIdx], 'bo',\
        periods[qatsPeriodIdx], bestFit[qatsPeriodIdx], 'bo')
ax.legend(('qats SNR', 'fit to qats', \
           'fit best period', 'qats best period'))

pylab.title('KID: %s' % kid)
pylab.ylabel('SNR')
pylab.xlabel('periods')

pylab.show()