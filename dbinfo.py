
dBhost = 'tddb.astro.washington.edu'
dBuser = 'tddb'
dBpass = 'tddb'
dBname = 'KeplerNew'
dBname0 = 'Kepler'

from uwpyKepler import iodb.dbConnect as dbConnect

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