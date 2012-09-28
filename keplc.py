
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
        self.lcData0 = iodb.ReadLightCurve(KID,selection=kw.ctype,data=kw.data)
        if not self.lcData0 == None:
            self.lcData1 = pipeline.FlagKeplerEvents(self.lcData0)
            self.lcData2 = pipeline.RemoveBadEvents(self.lcData1)
            self.lcData3 = pipeline.FlagEclipses(self.lcData2,eData,BJDREFI)
            self.lcData4 = pipeline.SplitPortions(self.lcData3,kw.gapSize)
            self.lcData5 = pipeline.FlagOutliers(self.lcData4,kw.oWin,kw.oThreshold)
            self.lcData6 = pipeline.DetrendChoice(self.lcData5,kw.dWin,kw.dPolyorder, kw.detChoice)
            self.lcData7 = pipeline.StackPortions(self.lcData6)
            self.lcData  = pipeline.CorrectNegVals(self.lcData7)
        else:
            raise NameError('No lightcurve data found')
        
class kw:
    
    def __init__(self, **kwargs):
        
        self.agap = None
        self.durfac = None
        self.detChoice = 'polynomial'
        self.data = 'pdc'
        self.maske = 'false'
        printString = ''
        for key in kwargs:
            if key == 'ctype': #cadence type, 'LC' = long cadence, 'SC' = short cadence data 
                self.ctype = kwargs[key]
                printString += ' ctype='+str(kwargs[key])+','
            elif key == 'gsize': #minimum gapsize (in days) to create a new data portion 
                self.gapSize = kwargs[key]
                printString += ' gsize='+str(kwargs[key])+','
            elif key == 'owin': #window, in cadences, used for computing the median for sigma clipping outliers
                self.oWin = kwargs[key]
                printString += ' owin='+str(kwargs[key])+','
            elif key.lower().startswith('othresh'): #sigma clipping threshold for outlier rejection
                self.oThreshold = kwargs[key]
                printString += ' othreshold='+str(kwargs[key])+','
            elif key == 'dwin': #detrending window size in cadences
                self.dWin = kwargs[key]
                printString += ' dwin='+str(kwargs[key])+','
            elif key == 'polyorder': #polynomial order to be used in detrend fitting
                self.dPolyorder = kwargs[key]
                printString += ' polyorder='+str(kwargs[key])+','
            elif key == 'agap': #unused?
                self.agap = kwargs[key]
                printString += ' agap='+str(kwargs[key])+','
            elif key == 'durationfactor': #unused?
                self.durfac = kwargs[key]
                printString += ' durationfactor='+str(kwargs[key])+','
            elif key == 'detChoice': #detrending method, a few to choose from, not much difference in results? pipeline.py has 'em
                self.detChoice = kwargs[key]
                printString += ' detChoice='+str(kwargs[key])+','
            elif key == 'data': #either use the 'pdc' (corrected) fluxes or the raw 'sap', pdc by default
                self.data = kwargs[key]
                printString += ' data='+str(kwargs[key])+','
            elif key == 'maske': #'true' = mask all eclipse times associated with KOIs in the database, 'false' = don't mask 'em
                self.maske = kwargs[key]
                printString += ' maske='+str(kwargs[key])+','
        self.printString = printString[:-2] 