import MySQLdb
import numpy as num
#import pylab

dBhost = 'tddb.astro.washington.edu'
dBuser = 'tddb'
dBpass = 'tddb'
dBname = 'Kepler'
 
def dbConnect(dBhost,dBuser,dBpass,dBname):
    """ Connect to server and return cursor """
    
    db     = MySQLdb.connect(host=dBhost, user=dBuser, passwd=dBpass, db=dBname)
    cursor = db.cursor()
    return cursor
    
def inSource(KeplerID):
    """ Checks if a certain KID exists in the source database. """
    
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    foo    = 'select KEPLERID from source where (KEPLERID = %s)' % (KeplerID)
    cursor.execute(foo)
    results = cursor.fetchall()
    Exist = False
    if len(results) > 0:
        Exist = True
    
    return Exist
        
def inKEPPC(KeplerID):
    """ Checks if a certain KID exists in the KEPPC database. """
    
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    foo1 = 'select KID from KEPPC where (KID = %s)' % (KeplerID)
    cursor.execute(foo1)
    results = cursor.fetchall()
    Exist = False
    if len(results) > 0:
        Exist = True
    
    return Exist

def inUWKOI(KeplerID):
    """ Checks if a certain KID exists in the UWKOI database. """
    
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    foo1 = 'select KID from UWKOI where (KID = %s)' % (KeplerID)
    cursor.execute(foo1)
    results = cursor.fetchall()
    Exist = False
    if len(results) > 0:
        Exist = True
    
    return Exist
    
def inKEPFP(KeplerID):
    """ Checks if a certain KID exists in the KEPFP database. """
    
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    foo1 = 'select KID from KEPFP where (KID = %s)' % (KeplerID)
    cursor.execute(foo1)
    results = cursor.fetchall()
    Exist = False
    if len(results) > 0:
        Exist = True
    
    return Exist

def getKOI(KeplerID):
    """ returns all KOI IDs (Kepler Object of Interest IDs) """
    
    # check for ID in Planet candidate and False positive table
    if inKEPPC(KeplerID):
        cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
        foo1 = 'select KOI from KEPPC where (KID = %s)' % (KeplerID)
        cursor.execute(foo1)
        results = cursor.fetchall()
    elif inKEPFP(KeplerID):
        cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
        foo1 = 'select KOI from KEPFP where (KID = %s)' % (KeplerID)
        cursor.execute(foo1)
        results = cursor.fetchall()
    elif inUWKOI(KeplerID):
        cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
        foo1 = 'select KOI from UWKOI where (KID = %s)' % (KeplerID)
        cursor.execute(foo1)
        results = cursor.fetchall()
    else:
        results = None
    
    return results

def ReadLightCurve(KeplerID, **kwargs):
    """ This function reads the Kepler database and returns the 
    corrected lightcurve in a dictionary
    
    Input = KeplerID
    Output = Dictionary with time, flux and fluxerror data
    """
    
    # the selection statement
    addition = ''
    for key in kwargs:
        if key == 'selection':
            if kwargs[key] == 'LC':
                addition = ' and LCFLAG = 1'
            elif kwargs[key] == 'SC':
                addition = ' and LCFLAG = 0'
            else:
                #setting default to LC
                addition = ' and LCFLAG = 1'
        else:
            #print key+' not recognized. using default'
            #setting default to LC
            addition = ' and LCFLAG = 1'
            continue
            
    foo = 'select * from source where (KEPLERID = %s' % (KeplerID)
    foo += addition+');'
    if inSource(KeplerID):
        cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
        cursor.execute(foo)
        results = cursor.fetchall()
        
        # reading time, corrected flux and flux errors
        cadence = num.array([x[1] for x in results])
        if len(cadence) == 0:
            #print 'Data not found for %s in Kepler.source for' % (KeplerID)
            #print ' %s option' % kwargs['selection']
            return
        else:
            time    = num.ma.array([x[2] for x in results])
            corflux = num.ma.array([x[5] for x in results])                                                                                                        
    	    corerr  = num.ma.array([x[6] for x in results])
            qflag  = num.array([x[7] for x in results]) 
            sortindex = cadence.argsort()
            time = time[sortindex]
            corflux = corflux[sortindex]
            corerr = corerr[sortindex]
            qflag = qflag[sortindex]
            cadence = cadence[sortindex]
            return {'kid':KeplerID,'x':time,'y':corflux,\
            'yerr':corerr,'qflag':qflag,'cadence':cadence}
    else:
        #print 'Kepler ID %s not found in Kepler.source' % (KeplerID)
        return

def getEclipseData(KeplerID):
    """ Queries the database for Transit or Eclipse data """
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    eData = {'KOI':{}}
    if inKEPPC(KeplerID):
        foo1 = 'select KOI, Period, Dur, Epoch from KEPPC where (KID = %s)' % KeplerID
        cursor.execute(foo1)
        r1 = cursor.fetchall()
	for i in range(len(r1)):
            eData['KOI'][format(r1[i][0],'.2f')] =\
            {'Period':r1[i][1],'Duration':r1[i][2],'T0':r1[i][3]}
        eData['eDataExists'] = True
    elif inKEPFP(KeplerID):
        foo1 = 'select KOI, Period, Duration, Epoch from KEPFP where (KID = %s)' % KeplerID
        cursor.execute(foo1)
        r1 = cursor.fetchall()
	for i in range(len(r1)):
            eData['KOI'][format(r1[i][0],'.2f')] =\
            {'Period':r1[i][1],'Duration':r1[i][2],'T0':r1[i][3]}
        eData['eDataExists'] = True
    elif inUWKOI(KeplerID):
        foo1 = 'select KOI, Period, Dur, Epoch from UWKOI where (KID = %s)' % KeplerID
        cursor.execute(foo1)
        r1 = cursor.fetchall()
	for i in range(len(r1)):
            eData['KOI'][format(r1[i][0],'.2f')] =\
            {'Period':r1[i][1],'Duration':r1[i][2],'T0':r1[i][3]}
        eData['eDataExists'] = True
        print "UWKOI", eData
    else:
        #print 'Kepler ID not found in Kepler.KEPPC or Kepler.KEPFP'
        eData['eDataExists'] = False
        
    return eData

def getBJDREFI(KeplerID):
    """ Queries the database for Transit or Eclipse data """
    
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    foo1 = 'select DISTINCT BJDREFI from object where (KEPLERID = %s)' % KeplerID
    cursor.execute(foo1)
    r1 = cursor.fetchall()
    if len(r1) == 1:
        return r1[0]
    elif len(r1) > 1:
        #print 'Multiple BJDREFIs found'
        return None
    else:
        #print 'No BJDREFIs found'
        return None
    
