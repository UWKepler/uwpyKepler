import MySQLdb
import sys
import numpy as num

    
def ReadLightCurve(KeplerID):
    """ This function reads the Kepler database and returns the 
    corrected lightcurve in a dictionary
    """

    db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
    cursor = db.cursor()
    foo    = 'select * from source where (KEPLERID = %s)' % (KeplerID)
    cursor.execute(foo)
    results = cursor.fetchall()
    
    time    = num.array([x[2] for x in results])
    rawflux = num.array([x[5] for x in results])
    rawerr  = num.array([x[6] for x in results])
    corflux = num.array([x[7] for x in results])
    corerr  = num.array([x[8] for x in results])
    
    if len(time) == 0:
        print 'No data on object %s' % (KeplerID)
        sys.exit(1)
    
    idx = []
    for i in range(len(time)):
        if corflux[i] != 0 and corerr[i] != 0: 
            idx.append(i)
    time    = time[idx]
    corflux = corflux[idx]
    corerr  = corerr[idx]
    
    return {'x':time,'y':corflux,'yerr':corerr}

#def ReadKOITable(input parameters):
    #""" This function reads the Kepler database and returns ... """


#def checkifKnownPlanet(input parameters):
    #""" returns True is Kepler ID matches a planet in the datatable ... """
    
    