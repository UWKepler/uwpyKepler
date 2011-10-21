#import difflib
import sys, os

def getUniqueTraceback(logFile):
    """
    sorts though exception log files and
    prints unique exceptions to separate files
    with a list of KIDs and input parameters.
    """
    
    file = open(logFile,'r')
    fileLines = file.readlines()
    FileTag = logFile[:-9]
    
    erri = 0
    StartStack = False
    errDict = {}
    objDict = {}
    stkDict = {}
    Stack = ''
    StackList = []

    # Read through each line in the log
    for i in range(len(fileLines)):
        # Look for the start of a new error message
        # and start writing the lines of this message
        if fileLines[i].startswith('Traceback'):
            errDict[erri] = Stack
            StackList.append(Stack)
            erri += 1
            StartStack = True
            Stack = ''
        # Look for the start of the KID and 
        # input parameter. note and write these lines
        if fileLines[i].startswith('#'):
            objDict[erri] = fileLines[i]
            StartStack = False
        #Stack error lines
        if StartStack:
            Stack += fileLines[i]

    # Get the set of unique bugs and write to file
    UniqueBugs = list(set(StackList[1:]))
    FOList = []
    for i in range(len(UniqueBugs)):
        # Failed object (KID) list file
        bFileName = FileTag+'.'+str(i+1).zfill(5)+'.Bug.log'
        os.system('rm -v %s' % bFileName)
        FOList.append(open(bFileName,'w'))
        print >> FOList[i],'%'+'-'*20+'%'
        print >> FOList[i],'%'+UniqueBugs[i]
        print >> FOList[i],'%'+'-'*20+'%'
        
    # remove useless key-value pair
    del errDict[0]
    
    for i in errDict.keys():
        for j in range(len(UniqueBugs)):
            if errDict[i] == UniqueBugs[j]:
                print >> FOList[j], objDict[i].strip('\n')
