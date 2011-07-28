
import uwpyKepler.io as io

class keplerlc:
    def __init__(self, KIDstring):
        self.kid = KIDstring
        self.dataFlag = io.inSource(KIDstring)
        self.PCFlag = io.inKEPPC(KIDstring)
        self.FPFlag = io.inKEPFP(KIDstring)
        self.KOI = io.getKOI(KIDstring)
        
