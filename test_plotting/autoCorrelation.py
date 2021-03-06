#!/astro/apps/pkg/python64/bin//python

import pylab as py
import numpy as num
import os
import uwpyKepler as kep
import sys

kid = sys.argv[1]

#getting data
#keplc = kep.keplc.keplc(kid)
#eData = keplc.eData
#BJDREFI = keplc.BJDREFI
#kw = kep.quicklc.quickKW(ctype='LC')
lc = kep.keplc.keplc(kid)
kw = kep.quicklc.quickKW()
lc.runPipeline(kw)
qlc = kep.qats.qatslc(lc.lcFinal, kid)
qlc.padLC()

x = qlc.lcData['x']
y = qlc.lcData['y']

#get lightcurve data after the pipeline has run on it
#lcData = kep.keplc.lcData(kid,eData,BJDREFI,kw).lcData
#x = lcData['x']+BJDREFI-2454900e0
#y = lcData['ydt']

#running autocorrelation function
crf = kep.func.autoCor(y)

#creating indexes of days and cadences to plot based off of, as well as chopping off silly end bits.
index = num.arange(len(crf))
daysIndex = (29.425*num.array(index))/(60.0*24.0)
chopped = crf[5:-5]
choppedIdx = index[5:-5]
daysChoppedIdx = (29.425*num.array(choppedIdx))/(60.0*24.0)
realx = num.array(x) - x[0]
xChopped = realx[5:-5]

#plotting results
py.title('Correllation function of KID: %s' % kid)
py.xlabel('Frequency')
py.ylabel('Correlation function')
py.plot(realx, crf, '-r')
py.plot(xChopped, chopped, '-b')
py.grid(True)
py.show()

