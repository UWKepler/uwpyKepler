

def readKIDsFromFile(File):
    """
    reads a list of KIDS from File
    """
    
    fileObj = open(File,'r')
    Out = {'File':File,'KIDList':[]}
    for line in fileObj.readlines():
        if len(line) > 0:
            Out['KIDList'].append(long(line.strip('\n')))
    
    return Out
        
        
def passFileToTrialList(PassLogFile):
    """
    Convert the pass entries from the pass-log file to a list of trails
    """
    
    FileObj = open(PassLogFile,'r')
    FileLines = FileObj.readlines()
    TrialList = []
    keyList = ['KIDList']
    for line in FileLines:
        Dict = {}
        Split1 = map(str, line.split('|'))
        valueList = [long(Split1[0].strip())]
        Split2 = map(str,Split1[1].split(','))
        for el in Split2:
            keypar = map(str,el.split('='))
            keyList.append(keypar[0].strip())
            try:
                valueList.append(long(keypar[1].strip()))
            except:
                valueList.append(keypar[1].strip())
                    
        TrialList.append(dict(zip(keyList,valueList)))
    
    return TrialList