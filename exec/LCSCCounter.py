#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep
import sys
import numpy as num

def returnLCSClists(KIDlistfile, KIDlist):
    dBhost = 'tddb.astro.washington.edu'
    dBuser = 'tddb'
    dBpass = 'tddb'
    dBname = 'KeplerNew'
    cursor = kep.iodb.dbConnect(dBhost,dBuser,dBpass,dBname)
    #sg = eval(KIDlistfile[2:5])
    #skycall = 'select distinct KEPLERID from\
               #object where SKYGROUP = %s and LCFLAG = 0' % (sg)
    #cursor.execute(skycall)
    #results = cursor.fetchall()
    #print len(results)
    LCkids = []
    SCkids = []
    LCSCkids = []
    count = 0

    for KID in KIDlist:
        count += 1
        if float(count)%100 == 0:
            print count
        foo0 = 'select KEPLERID from source where\
           (KEPLERID = %s and LCFLAG = 0)' % (KID)
        foo1 = 'select KEPLERID from source where\
           (KEPLERID = %s and LCFLAG = 1)' % (KID)
        cursor.execute(foo0)
        scresults = cursor.fetchall()
        cursor.execute(foo1)
        lcresults = cursor.fetchall()
        if len(scresults) > 0:
            SCkids.append(KID)
        if len(lcresults) > 0:
            LCkids.append(KID)
        if len(scresults) > 0 and len(lcresults) > 0:
            LCSCkids.append(KID)
        if count == 50:
            break

    return SCkids, LCkids, LCSCkids



if __name__ == '__main__':
    print 'format:'
    print 'KIDListFile | KIDs | LC KIDs | SC KIDs | LCSC KIDs'

    KIDListFile = sys.argv[1]
    skygroup = open(KIDListFile, 'r')
    KIDs = skygroup.readlines()
    SCkids, LCkids, LCSCkids = returnLCSClists(KIDListFile, KIDs)

    print KIDListFile,'|',len(KIDs),'|',len(LCkids),'|',len(SCkids),'|',len(LCSCkids)






