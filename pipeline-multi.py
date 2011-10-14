import numpy as num

def readPipelinePars(File):
    """
    reads Pipeline Parameter Files
    """

    DataFile = open(File,'r')
    ParDict = {}
    for line in DataFile.readline():
        if not line.startswith('#'):
            DataArray = map(str, line.split('|'))

    return ParDict

def SortParOption(OptionString):
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
    
    return List

class PipelinePars:

    def __init__(self, FileName):
        
        ParDict = readPipelinePars(FileName)
        for key in ParDict.keys():
            if key != 'parsetName' or key != 'cChoice':
                ListData = SortParOptions(ParDict[key])
                setattr(self,key,ListData)
            elif key == 'cChoice':
                ListData = map(str,ParDict[key].split(','))
                if not isinstance(ListData,list):
                    ListData = [ListData]
                setattr(self,key,ListData)
            else:
                setattr(self,key,ParDict[key])
                 