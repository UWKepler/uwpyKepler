import sys
import uwpyKepler as kep
import numpy as num
import pylab
import MySQLdb

#KeplerID = sys.argv[1]
KeplerID = 8478994
#KeplerID = 11295426
KeplerID = 10341831

db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
cursor = db.cursor()
foo    = 'select * from source where (KEPLERID = %s)' % (KeplerID)
cursor.execute(foo)
results = cursor.fetchall()

# reading time, corrected flux and flux errors
time    = num.ma.array([x[2] for x in results])
f1 = num.ma.array([x[5] for x in results])
e1  = num.ma.array([x[6] for x in results])
f2 = num.ma.array([x[7] for x in results])
e2  = num.ma.array([x[8] for x in results])

btime = num.where(time <= 0)
bcf = num.where(f2 <= 0)
bcerr = num.where(e2 <= 0)

good = num.where( (f2 > 0) & (e2 > 0))

print time[btime]
print f2[bcf]
print e2[bcerr]

print btime
print bcf
print bcerr

pylab.plot(time[good],f2[good]-num.median(f2),'b.')
pylab.plot(time[bcf],f2[bcf],'r.')
pylab.show()
