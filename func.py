import numpy as num
import pylab
"contains functions"

def getNumBool(List,bool):
    """ 
    Returns the total number of booleans in a list
    bool should be True or False
    
    e.g.
    print len(getNumBool(List,True)), will list the number
    of elements which are True in List.
    """
    
    List = num.array(List)
    id0 = num.where(List == bool)[0]
    #print len(id0)
    return len(id0)

def makeDTwindows(nsize,window):
    """ Given a portion length (nsize) and window size (window),
        this function creates appropriate indicies
        used to sclice windows for detrending.
    """
    
    # two sets of slices for the cos**2 and sin**2 weights
    set1 = num.arange(0,nsize,window)
    set2 = num.arange(0+window/2,nsize,window)
    set2 = num.hstack( (0,set2) )
    
    outwin = {}
    #diff, computes excess points near the edge of a portion
    diff1 = nsize - set1[-1]
    diff2 = nsize - set2[-1]
    if diff1 < window/2:
        set1[-1] = set1[-1]+diff1	#adding excess points to index
    else:
        set1 = num.hstack((set1,nsize)) #stacking slice index pairs
    if diff2 < window/2:
        set2[-1] = set2[-1]+diff2	#adding excess points to index
    else:
        set2 = num.hstack((set2,nsize)) #stacking slice index pairs

    # for the special case when only 1 portion needs to be fit
    if len(set1) == 1:
        set1 = num.hstack((0,set1))
    if len(set2) == 1:
        set2 = num.hstack((0,set2))

    # make pairs of the slicing indices for easier use
    outSet1 = makeRangePairs(set1)
    outSet2 = makeRangePairs(set2)
    idict = {1:outSet1,2:outSet2}
 
    return idict

def makeRangePairs(indexList):
    """
    Given a list of indicies make a list of pairs
    using sequential indicies.
    """
    
    set = []
    for i in range(len(indexList)-1):
        set.append([indexList[i],indexList[i+1]])

    return set

def movingMedian(data,window):
    """
    Returns the median smoothed version of given data.
    data is the data (usually time series)
    window is the size of the window used for smoothing
    """
    
    mvavg = []
    npoints = len(data)
    for i in range(npoints-1):
        j1  = max(0,i-window)
        j2 = min(i+window,npoints-1)
        mvavg.append(num.median(data[j1:j2]))
        
    return num.array(mvavg)
    
def compute1Sigma(data):
    """
    return 1-sigma
    """
    
    dsort = num.sort(data)
    npts = len(dsort)
    sigma = (dsort[.8415*npts]-dsort[.1585*npts])/2
    return sigma
                   

def foldPhase(xdata,t0,period):
    """ enter duration dur in hours and t0 is center time of transit"""
    
    t0 = t0 + 2454900e0
    phase = (xdata-t0)/period - (xdata-t0)//period

    #pylab.plot(phase,data['ydt'],'b.')
    #pylab.show()

    return phase
