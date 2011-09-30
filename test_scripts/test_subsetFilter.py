import uwpyKepler as kep
import pylab
import numpy as num
import sys

stars = [6056992,7849854,8013419,6448431,6448574,6448651,6448299,6448304]

filteredList = kep.dbinfo.returnKIDsNonKOI(stars)

print 'unfiltered list', stars
print 'filtered list', filteredList


