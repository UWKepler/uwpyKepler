

def readKIDsFromFile(File):
    """
    reads a list of KIDS from File
    """
    
    fileObj = open(File,'r')
    Out = {'File':File,'KIDlist':[]}
    for line in fileObj.readline():
        Out['KIDlist'].append(map(float, line.strip('\n')))
    
    return Out