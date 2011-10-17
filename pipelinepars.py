import itertools
import numpy as num

def readPipelinePars(File):
    """
    reads Pipeline Parameter Files
    """

    DataFile = open(File,'r')
    ParDict = {}
    for line in DataFile.readlines():
        if not line.startswith('#'):
            DataArray = map(str, line.split('|'))
            key = DataArray[0].strip()
            par = DataArray[1].strip('\n').strip()
            ParDict[key] = par

    return ParDict

def SortParOptions(OptionString):
    """
    Convert Option string to a list of parameters
    """
    
    List = []
    if OptionString.strip().startswith('ar'):
        Dummy_array = map(str,OptionString.split(':'))
        ANumbers = map(long,Dummy_array[1].split(','))
        try:
            List = num.arange(ANumbers[0],ANumbers[1],ANumbers[2])
        except:
            raise NameError("Bad arange number options '"+OptionString+"'")
    else:
        List = map(long,OptionString.split(','))
        if not isinstance(List,list):
            List = [List]
    
    return list(set(List))

def computeNTrials(ParObject):
    NTrials = 1
    for attr in ParObject.ParList:
        NTrials *= len(getattr(ParObject,attr))

class PipelinePars:

    def __init__(self, FileName):
        ParDict = readPipelinePars(FileName)
        ParList = []
        for key in ParDict.keys():
            if key != 'parsetName' and key != 'ctype':
                ListData = SortParOptions(ParDict[key])
                setattr(self.__class__,key,ListData)
                ParList.append(key)
            elif key == 'ctype':
                ListData = map(str,ParDict[key].split(','))
                if not isinstance(ListData,list):
                    ListData = [ListData]
                setattr(self.__class__,key,ListData)
                ParList.append(key)
            else:
                setattr(self.__class__,key,ParDict[key])
        self.ParList = ParList
        self.NKeyWordOptions = computeNTrials(self)
    
    def Trials(self, KIDList):
        self.NKIDs = len(KIDList)
        UberList = {'KIDList':KIDList}
        for parName in self.ParList:
            UberList[parName] = getattr(self,parName)
        Trials = TrialList(UberList)
        self.NTrials = len(Trials)
        
        return Trials

def TrialList(UberList):
    
    keyList = []
    InList = []
    TList = []
    for key in UberList.keys():
        keyList.append(key)
        InList.append(UberList[key])

    ComboList = itertools.product(*InList)
    for combo in ComboList:
        TList.append(trial(dict(zip(keyList,combo))))
        
    return TList

class trial:
    
    def __init__(self,dict):
        self.kid = dict['KIDList']
        del dict['KIDList']
        self.kw = dict