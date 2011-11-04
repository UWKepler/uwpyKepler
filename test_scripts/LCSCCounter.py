#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import sys
import numpy as num

dBhost = 'tddb.astro.washington.edu'
dBuser = 'tddb'
dBpass = 'tddb'
dBname = 'KeplerNew'
def ReturnCadence(KeplerID, **kwargs):
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
     
     
     
     
    #cursor = dbConnect(dBhost,dBuser,dBpass,dBname)
    #command = 'select distinct KEPLERID from \
    #object where SKYGROUP = %s' % long(SGNumber)
    #cursor.execute(command)
    #results = cursor.fetchall()
    #r1 = ["%s" % el[0] for el in results]
    
    #return r1
    
    foo = 'select * from source where (KEPLERID = %s' % (KeplerID)
    foo += addition+');'
    if kep.iodb.inSource(KeplerID):
        cursor = kep.iodb.dbConnect(dBhost,dBuser,dBpass,dBname)
        cursor.execute(foo)
        results = cursor.fetchall()
        
        # reading time, corrected flux and flux errors
        cadence = num.array([x[1] for x in results])
        
        return cadence
    

SGNumber = sys.argv[1]
LCkids = []
SCkids = []
SCLCkids = []
#skygroup = kep.dbinfo.returnSkyGroupKIDs(SGNumber)

if int(SGNumber) < 10:
    strSGNumber = '0'+str(SGNumber)
else:
    strSGNumber = SGNumber
skygroup = open('SG0'+strSGNumber+'_NoKOIs.list', 'r')

KIDs = skygroup.readlines()
for KID in KIDs:
    #if count/200. == count//200. and count != 0:
        #print 200+count*200, " KIDs checked"
    KID = int(KID)
    LCcadence = ReturnCadence(KID, selection='LC')
    if len(LCcadence) > 0:
        LCkids.append(KID)
    SCcadence = ReturnCadence(KID, selection='SC')
    if len(SCcadence) > 0:
        SCkids.append(KID)
for KID in SCkids:
    try:
        LCkids.index(KID)
    except:
        continue
    else:
        SCLCkids.append(KID)

file = open('countresults.txt', 'a')
print >> file, SGNumber,' | ',len(KIDs),' | ',len(LCkids),' | ', len(SCkids),' | ',len(SCLCkids)
    
