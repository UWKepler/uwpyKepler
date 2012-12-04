import MySQLdb
import numpy as num
import pylab

dBhost = 'tddb.astro.washington.edu'
dBuser = 'tddb'
dBpass = 'tddb'
dBname = 'Kepler'
db     = MySQLdb.connect(host=dBhost, user=dBuser, passwd=dBpass, db=dBname)
cursor = db.cursor()
 
kid  = 10004519

sql  = "SELECT runs.width, runs.depth, results.qatsVals FROM"
sql += "   qatsRuns AS runs"
sql += " JOIN"
sql += "   qatsResults AS results"
sql += " ON"
sql += "   runs.qatsValues = results.qatsId"
sql += " WHERE"
sql += "   runs.KID = %s" % (kid)
print sql

cursor.execute(sql)
results = cursor.fetchall()

widths = list(set([x[0] for x in results]))
widths.sort()
depths = list(set([x[1] for x in results]))
depths.sort()
nqats  = num.array(results[0][2].split(',')).shape[0]

qatsMatrix = num.ndarray((len(widths), len(depths), nqats))

for result in results:
    width = result[0]
    depth = result[1]
    qats  = num.array(result[2].split(',')).astype(num.float)
    if qats.shape[0] != nqats:
        print "UH, MIXING DIFFERENT SHAPE QATS PLOTS"
        continue
    widx  = widths.index(width)
    didx  = depths.index(depth)
    qatsMatrix[widx][didx] = qats

pylab.imshow(qatsMatrix[:,:,99], aspect="auto", origin="lower", extent=(widths[0], widths[-1], depths[0], depths[-1]))
pylab.show()
import pdb; pdb.set_trace()

