

import optparse

class keplc:
    
    def __init__(self,KeplerID):
        self.KID = KeplerID
        self.inSource = iodb.inSource(KeplerID)
        self.inKEPPC = iodb.inKEPPC(KeplerID)
        self.inKEPFP = iodb.inKEPFC(KeplerID)
        self.eData = iodb.getEclipseData(KeplerID)
        self.BJDREFI = iodb.getBJDREFI(KeplerID)
    
    def runPipeline(self, kw):
        
        x = lcData(kw)
        
        self.lcData = x
        self.lcFinal = x.lcData
        
class lcData:
    
    def __init__(self, kw):
        self.lcData0 = iodb.ReadLightCurve(self.KID,selection=kw.ctype)
        if not lcData0 == None:
            self.lcData1 = pipeline.FlagKeplerEvents(self.lcData0)
            self.lcData2 = pipeline.RemoveBadEvents(self.lcData1)
            self.lcData3 = pipeline.FlagEclipses(self.lcData2,self.eData,self.BJDREFI)
            self.lcData4 = kep.pipeline.SplitPortions(self.lcData3,kw.gapSize)
            self.lcData5 = kep.pipeline.FlagOutliers(self.lcData4,kw.oWin,kw.oThreshold)
            self.lcData6 = kep.pipeline.DetrendData(self.lcData5,kw.dWin,kw.dPolyorder)
            self.lcData = kep.pipeline.StackPortions(self.lcData6)
        else:
            raise NameError('No lightcurve data found for '+str(self.KID))
        
class kw:
    
    def __init__(self, **kwargs):
        
        self.agap = None
        self.durf = None
        for key in kwargs:
            if key == 'ctype':
                self.ctype = kwargs[key]
            elif key == 'gsize':
                self.gapSize = kwargs[key]
            elif key == 'owin':
                self.oWin = kwargs[key]
            elif key == 'othresh':
                self.oThreshold = kwargs[key]
            elif key == 'dwin':
                self.dWin = kwargs[key]
            elif key == 'polyorder':
                self.polyorder = kwargs[key]
            elif key == 'agap':
                self.aGap = kwargs[key]
            elif key == 'durf':
                self.durf = kwargs[key]
                