
dBhost = 'tddb.astro.washington.edu'
dBuser = 'tddb'
dBpass = 'tddb'
dBname = 'KeplerNew'
dBname0 = 'Kepler'
import uwpyKepler as kep
from uwpyKepler import iodb.dbConnect as dbConnect

def returnKIDsNonKOI(starList):
    filteredList = []
    for i in range(len(starList)):
        FPflag = kep.iodb.inKEPFP(starList[i])
        PCflag =kep.iodb.inKEPPC(starList[i])
        if (FPflag == False) & (PCflag == False):
            filteredList.append(starList[i])

    return filteredList

def getKIDsFP():
    """
    Gets all distinct KIDs from False positive table
    """
    
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname0)
    foo1 = 'select DISTINCT KID from KEPFP'
    cursor.execute(foo1)
    results = cursor.fetchall()
    r1 = ["%s" % el[0] for el in results]
    
    return r1

def getKIDsPC():
    """
    Gets all distinct KIDs from the Planet Candidates table
    """
    
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname0)
    foo1 = 'select DISTINCT KID from KEPPC'
    cursor.execute(foo1)
    results = cursor.fetchall()
    r1 = ["%s" % el[0] for el in results]

    return r1

def getKIDsSource():
    """
    Gets all distinct KIDs from the source table
    """
    
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    foo1 = 'select distinct KEPLERID from source'
    cursor.execute(foo1)
    results = cursor.fetchall()
    r1 = ["%s" % el[0] for el in results]
    
    return r1

def returnSkyGroupKIDs(SGNumber):
    """
    Gets all distinct KIDs from the source table from \
    a given Sky Group (SGNumber)
    Kepler SGNumbers go from 1 to 84
    """

    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    command = 'select distinct KEPLERID from \
    object where SKYGROUP = %s' % long(SGNumber)
    cursor.execute(command)
    results = cursor.fetchall()
    r1 = ["%s" % el[0] for el in results]
    
    return r1

def retrunCoordKID(KIDlist):
    """
    Returns RA and DEC for a given list of KIDs
    """
    
    RA = []
    DEC = []
    for KID in KIDlist:
        cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
        command = 'select distinct RA_OBJ, DEC_OBJ from object where \
        KEPLERID = %s' % KID
        cursor.execute(command)
        r1 = cursor.fetchall()
        RA.append(r1[0][0])
        DEC.append(r1[0][1])

    return RA, DEC
    
def returnKIDsInKEPFP(KIDlist):

    """1. iterate through KIDs in list
    2. keep a list of only those KIDs which are are in KEPFP
    return list """

    KIDFP = []
    for kid in KIDlist:
        result = kep.iodb.inKEPFP(kid)
        if result == True:
            KIDFP.append(kid)
    return KIDFP


def returnKIDsInKEPPC(KIDlist):

    """1. iterate through KIDs in list
    2. keep a list of only those KIDs which are are in KEPPC
    return list """

    KIDPC = []
    for kid in KIDlist:
        result = kep.iodb.inKEPPC(kid)
        if result == True:
            KIDPC.append(kid)
    return KIDPC

