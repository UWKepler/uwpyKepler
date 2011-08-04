
import uwpyKepler.io as io

class keplerlc:
    def __init__(self, KIDstring):
        self.kid = KIDstring
        self.dataFlag = io.inSource(KIDstring)
        self.PCFlag = io.inKEPPC(KIDstring)
        self.FPFlag = io.inKEPFP(KIDstring)
        self.KOI = io.getKOI(KIDstring) # getKOI needs to work with io.getEclipseData 
    def data(self):
        # fix this -- 
        data1 = io.ReadLightCurve(self.kid)
        self.xraw = data['x']
        self.yraw = data['y']
        self.yerrraw = data['yerr']
        return data1
    def eData(self):
        eData = io
        
d1 = kep.io.ReadLightCurve(kid)
eData = kep.io.getEclipseData(d1)
d2 = kep.io.FlagTransits(d1,eData)
pd = kep.io.SplitGap(d2,0.1)
d3 = kep.io.FlagOutliers(pd,10,4)
