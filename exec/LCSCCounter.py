#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import sys
import numpy as num

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

    dBhost = 'tddb.astro.washington.edu'
    dBuser = 'tddb'
    dBpass = 'tddb'
    dBname = 'KeplerNew'
    foo = 'select DISTINCT KEPLERID from source where (KEPLERID = %s' % (KeplerID)
    foo += addition+');'
    if kep.iodb.inSource(KeplerID):
        cursor = kep.iodb.dbConnect(dBhost,dBuser,dBpass,dBname)
        cursor.execute(foo)
        results = cursor.fetchall()
        return results
    
if __name__ == '__main__':

    KIDListFile = sys.argv[1]
    skygroup = open(KIDListFile, 'r')
    
    #count = 0
    LCkids = []
    SCkids = []
    SCLCkids = []
    KIDs = skygroup.readlines()
    for KID in KIDs:
        #count += 1
        #if float(count)/100. - float(count)//100. == 0:
            #print count
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
    
    print KIDListFile
    print 'N KIDs = ', len(KIDs)
    print 'N LC KIDs = ', len(LCkids)
    print 'N SC KIDs = ', len(SCkids)
    print 'N LC & SC KIDs = ', len(SCLCkids) 
    #file = open('countresults.txt', 'a')
    #print >> file, SGNumber,' | ',len(KIDs),' | ',len(LCkids),' | ', len(SCkids),' | ',len(SCLCkids)
    
