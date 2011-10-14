
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

def genPipelineParLists(ParDict):
    """
    generates a key word argument string for runPipeline
    """
    
    for key in ParDict.keys():
        if key != 'parsetName' or key != 'cChoice':


class PipelinePars:

    def __init__(self, FileName):
        
        
