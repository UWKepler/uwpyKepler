dBhost = 'tddb.astro.washington.edu'
dBuser = 'tddb'
dBpass = 'tddb'
dBname = 'Kepler'

from iodb import dbConnect, inKEPFP, inKEPPC

def getKIDsFP():
    """
    Gets all distinct KIDs from False positive table
    """
    
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    foo1 = 'select DISTINCT KID from KEPFP'
    cursor.execute(foo1)
    results = cursor.fetchall()
    r1 = ["%s" % el[0] for el in results]
    
    return r1

def getKIDsPC():
    """
    Gets all distinct KIDs from the Planet Candidates table
    """
    
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
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

def returnSkyGroupKIDs(SGNumber, **kwargs):
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
    
    for kw in kwargs:
	if kw == 'dict_output':
	    if kwargs[kw]:
		return {'SG' + str(SGNumber).zfill(3): r1}
    return r1

def returnCoordKID(KIDlist):
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
        #print r1, KID
        # temporary fix to ISSUE26
        if len(r1)> 0:
            RA.append(r1[0][0])
            DEC.append(r1[0][1])

    return RA, DEC
    
def returnKIDsInKEPFP(KIDlist):
    """ Iterates through KIDs in list
        and returns a list of KIDs in KEPFP.
    """

    KIDFP = []
    for kid in KIDlist:
        result = inKEPFP(kid)
        if result == True:
            KIDFP.append(kid)
    return KIDFP

def returnKIDsInKEPPC(KIDlist):
    """ Iterates through KIDs in list
        and returns a list of KIDs in KEPPC.
    """

    KIDPC = []
    for kid in KIDlist:
        result = inKEPPC(kid)
        if result == True:
            KIDPC.append(kid)
    return KIDPC

def returnKIDsNonKOI(KIDlist):
    """ Iterates through KIDs in list
        and returns list of KIDs not in KEPPF and KEPPC.
    """
    
    filteredList = []
    for kid in KIDlist:
        FPflag = inKEPFP(kid)
        PCflag = inKEPPC(kid)
        if (FPflag == False) & (PCflag == False):
            filteredList.append(kid)

    return filteredList
        
def returnRstar(KID):
    """
    Returns Rstar for a given KID
    """
    
    #Rstar = []
    #for KID in KIDlist:
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    command = 'select distinct RADIUS from object where \
    KEPLERID = %s' % KID
    cursor.execute(command)
    r1 = cursor.fetchall()
    # temporary fix to ISSUE26
    if len(r1)> 0:
        #Rstar.append(r1[0])
        Rstar = r1[0][0]
        
    return float(Rstar)
        
def returnLOGG(KID):
    """
    Returns LOGG for a given KID
    """
    
    #LOGG = []
    #for KID in KIDlist:
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    command = 'select distinct LOGG from object where \
    KEPLERID = %s' % KID
    cursor.execute(command)
    r1 = cursor.fetchall()
    # temporary fix to ISSUE26
    if len(r1)> 0:
        #LOGG.append(r1[0])
        LOGG = r1[0][0]

    return float(LOGG)
    
def getSkyGroup(KID):
    cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    command = 'select skygroup from object where' \
    ' KEPLERID = %s' % KID
    cursor.execute(command)
    r1 = cursor.fetchall()
    
    return r1[0][0]
