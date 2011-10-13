#!/astro/apps/pkg/python64/bin//python

import uwpyKepler as kep

if __name__ == '__main__':
    
    SGList = range(1,85,1)
    ALL_KIDsFP = []
    ALL_KIDsPC = []
    ALL_KIDs = []
    
    FileALLKIDs = open('ALL_KIDs.list','w')
    FileALLFPs = open('ALL_FPs.list','w')
    FileALLPCs = open('ALL_PCs.list','w')
    FileALLSinglePCs = open('ALL_SinglePCs.list','w')
    FileALLMultiPCs = open('ALL_MultiPCs.list','w')
    for SG in SGList:
        print 'printing SG # = ',SG 
        KIDlist = kep.dbinfo.returnSkyGroupKIDs(SG)
        FileSGNoKOI = open('SG'+str(SG).zfill(3)+'_NoKOIs.list','w')
        for KID in KIDlist:
            print KID, SG
            print >> FileALLKIDs, KID
            if kep.iodb.inKEPFP(KID):
                print >> FileALLFPs, KID
            elif kep.iodb.inKEPPC(KID):
                KOIs = kep.iodb.getKOI(KID)
                if len(KOIs) > 1:
                    print >> FileALLMultiPCs, KID
                elif len(KOIs) == 1:
                    print >> FileALLSinglePCs, KID
                else:
                    pass
            else:
                print >> FileSGNoKOI, KID
        FileSGNoKOI.close()
