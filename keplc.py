
import optparse
import iodb, pipeline

class keplc:
    
    def __init__(self,KeplerID):
        self.KID = KeplerID
        self.inSource = iodb.inSource(KeplerID)
        self.inKEPPC = iodb.inKEPPC(KeplerID)
        self.inKEPFP = iodb.inKEPFP(KeplerID)
        self.eData = iodb.getEclipseData(KeplerID)
        self.BJDREFI = iodb.getBJDREFI(KeplerID)
    
    def runPipeline(self, kw):
        
        x = lcData(self.KID,self.eData,self.BJDREFI,kw)
        self.lcData = x
        self.lcFinal = x.lcData
        
class lcData:
    
    def __init__(self, KID,eData,BJDREFI, kw):
        self.lcData0 = iodb.ReadLightCurve(KID,selection=kw.ctype)
        if not self.lcData0 == None:
            self.lcData1 = pipeline.FlagKeplerEvents(self.lcData0)
            self.lcData2 = pipeline.RemoveBadEvents(self.lcData1)
            self.lcData3 = pipeline.FlagEclipses(self.lcData2,eData,BJDREFI)
            self.lcData4 = pipeline.SplitPortions(self.lcData3,kw.gapSize)
            self.lcData5 = pipeline.FlagOutliers(self.lcData4,kw.oWin,kw.oThreshold)
            self.lcData6 = pipeline.DetrendData(self.lcData5,kw.dWin,kw.dPolyorder)
            self.lcData = pipeline.StackPortions(self.lcData6)
        else:
            raise NameError('No lightcurve data found for '+str(KID))
        
class kw:
    
    def __init__(self, **kwargs):
        
        self.agap = None
        self.durfac = None
        printString = ''
        for key in kwargs:
            if key == 'ctype':
                self.ctype = kwargs[key]
                printString += ' ctype='+str(kwargs[key])+','
            elif key == 'gsize':
                self.gapSize = kwargs[key]
                printString += ' gsize='+str(kwargs[key])+','
            elif key == 'owin':
                self.oWin = kwargs[key]
                printString += ' owin='+str(kwargs[key])+','
            elif key == 'othresh':
                self.oThreshold = kwargs[key]
                printString += ' othreshold='+str(kwargs[key])+','
            elif key == 'dwin':
                self.dWin = kwargs[key]
                printString += ' dwin='+str(kwargs[key])+','
            elif key == 'polyorder':
                self.dPolyorder = kwargs[key]
                printString += ' polyorder='+str(kwargs[key])+','
            elif key == 'agap':
                self.agap = kwargs[key]
                printString += ' agap='+str(kwargs[key])+','
            elif key == 'durationfactor':
                self.durfac = kwargs[key]
                printString += ' durationfactor='+str(kwargs[key])+','

        self.printString = printString[:-1] 