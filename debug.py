#import difflib
import sys, os

def getUniqueTraceback(logFile):
    
    file = open(logFile,'r')
    fileLines = file.readlines()
    
    erri = 0
    StartStack = False
    errDict = {}
    objDict = {}
    stkDict = {}
    Stack = ''
    StackList = []
    for i in range(len(fileLines)):
        if fileLines[i].startswith('Traceback'):
            errDict[erri] = Stack
            stkDict[erri] = Stack
            StackList.append(Stack)
            erri += 1
            StartStack = True
            Stack = ''
        if fileLines[i].startswith('#'):
            objDict[erri] = fileLines[i]
            StartStack = False
        if StartStack:
            Stack += fileLines[i]

    print len(set(StackList[1:])), len(StackList[1:])
    UniqueBugs = list(set(StackList[1:]))
    UBugFile = open('UniqueBug.log','w')
    FOList = []
    for i in range(len(UniqueBugs)):
        print >> UBugFile, '# '+str(i+1)
        print >> UBugFile, UniqueBugs[i]
        bFileName = 'OBJ_BUG.'+str(i+1)+'.log'
        os.system('rm -v %s' % bFileName)
        FOList.append(open(bFileName,'a'))
        
    UBugFile.close()
    del errDict[0]
    
    for i in errDict.keys():
        for j in range(len(UniqueBugs)):
            if errDict[i] == UniqueBugs[j]:
                print >> FOList[j], objDict[i].strip('\n')
