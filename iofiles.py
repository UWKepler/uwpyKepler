

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