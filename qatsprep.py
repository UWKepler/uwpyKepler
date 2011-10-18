import numpy as num
from lcmod import returnData

def getListIndicies(Array,ListValues):
    """
    Return a list of indices where a list of values exists
    in a given array
    """ 
    
    lArray = Array.tolist()
    lLV = ListValues.tolist()
    
    indX = [lArray.index(x) if x in lLV else None for x in lLV]
    
    return indX

def padLC(lcData,FlagIds):
    """
    Pads missing cadences with ones.
    Returns the padded lightcurve.
    """
    
    # the eclipse lightcurve data are the best place to start
    x,y,yerr,cad = returnData(lcData,'elc')
    
    if len(FlagIds) > 0:
        if len(x)-1 < max(FlagIds):
            raise NameError("len(x) < max(FlagIds) "+\
            str(len(x))+" < "+str(len(FlagIds)) )

    x0 = min(x)
    x1 = max(x)
    c0 = min(cad)
    c1 = max(cad)
    diffx = num.diff(x)
    Tcad = num.median(diffx)
    xcomplete = num.arange(x0,x1,Tcad)
    cadcomplete = num.arange(c0,c1+1,1)
    
    if len(xcomplete) != len(cadcomplete):
        print len(xcomplete), len(cadcomplete)
        raise NameError("mismatch in cadence and time intervals") 
    
    #missing = num.array(list(set.difference(set(cadcomplete),set(cad))))
    #missingIDX = getListIndicies(cadcomplete,missing)
    existingIDX = getListIndicies(cadcomplete,cad)
    print len(missingIDX),len(existingIDX),len(cadcomplete),len(missingIDX)+len(existingIDX)

    zeros = num.zeros(len(xcomplete))
    yerrcomplete = zeros
    ycomplete = zeros+1e0                    #to be used for y
    
    xcomplete[existingIDX] = x
    ycomplete[existingIDX] = y
    yerrcomplete[existingIDX] = yerr
    
    print x0, min(xcomplete)
    print x1, max(xcomplete)
    print c0, min(cadcomplete)
    print c1, max(cadcomplete)

    
    return
