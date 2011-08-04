
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

