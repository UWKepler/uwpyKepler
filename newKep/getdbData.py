import MySQLdb
import numpy as num

dBhost = 'tddb.astro.washington.edu'
dBuser = 'tddb'
dBpass = 'tddb'
dBname = 'Kepler'

# query db; return widths, depths, corresponding qats values
def qatsResults_all():
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
    
    cursor.execute(sql)
    results = cursor.fetchall()
    
    widths   = num.array([x[0] for x in results])
    depths   = num.array([x[1] for x in results])
    qatsVals = [x[2].split(',') for x in results]
    for i in range(len(qatsVals)):
        qatsVals[i] = num.array(qatsVals[i]).astype(num.float)
    qatsVals = num.array(qatsVals)
    
    return widths, depths, qatsVals

# create qatsMatrix:
# 1st dimension is transit width
# 2nd dimension is transit depth
# 3rd dimension is qats values
# full_output returns transit widths and depths corresponding
# to indices of qatsMatrix
def qatsResults_3d(full_output = False):
    widths, depths, qatsVals = qatsResults_all()
    #widthDepthMap = dict(zip(qatsVals, zip(widths, depths)))
    widthsSort = list(set(widths))
    depthsSort = list(set(depths))
    widthsSort.sort()
    depthsSort.sort()
    qatsMatrix = num.zeros((len(widthsSort), len(depthsSort), len(qatsVals[0])))
    for i in range(len(qatsVals)):
        widx = widthsSort.index(widths[i])
        didx = depthsSort.index(depths[i])
        qatsMatrix[widx][didx] = qatsVals[i]
    
    if full_output:
        widthsSort = num.array(widthsSort)
        depthsSort = num.array(depthsSort)
        return qatsMatrix, widthsSort, depthsSort
    else:
        return qatsMatrix

# determine and return which qats spectrum has the highest snr
def getTopSpectrum(qMatrix):
    maxes = qMatrix.max(axis=2)
    imax  = maxes.argmax()
    snr = qMatrix[imax / maxes.shape[1], imax % maxes.shape[1]]
    # trim off trailing zeros
    return snr[snr != 0]

# accurate to at least first decimal
def reconstructPeriodArray(snr, f, cadence = 29.42):
    # uses p0 * (1 + f) for the first bin
    pmins, pmaxes = periodSynthesis(snr, f, cadence=cadence)

    # we need the midpoints of the period bins
    periods = (pmins + pmaxes) / 2

    return periods

def periodSynthesis(snr, f, cadence = 29.42):
    # trim snr array to avoid unnecessary computations
    snr = snr[snr != 0]
    # use the standard Kepler value for 29.42-minute cadences
    cadence_in_days = cadence / 60 / 24
    # define a function to round up to nearest cadence
    cadRound = lambda p : num.ceil(p / cadence_in_days) * cadence_in_days
    
    pmins = [1.]
    # the first period bin uses 'f' instead of 'f/2'
    pmins.append( cadRound(pmins[-1] * (1 + f)) )
    for idx in range(len(snr) - 1):
        pmins.append( cadRound(pmins[-1] * (1 + f/2)) )
    pmins = num.array(pmins)
    pmaxes = pmins[1:]
    pmins = pmins[:-1]
    
    return pmins, pmaxes
