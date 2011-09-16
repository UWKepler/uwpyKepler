import numpy as num


def getNumBool(List,bool):
    """ 
    Returns the total number of certain booleans in a list
    """
    List = num.array(List)
    id0 = num.where(List == bool)[0]
    #print len(id0)
    return len(id0)

def makeDTwindows(nsize,window):
    """ Given a portion length and window size,
        this function creates appropriate for detrending.
    """
    
    set1 = num.arange(0,nsize,window)
    set2 = num.arange(0+window/2,nsize,window)
    set2 = num.hstack( (0,set2) )
    x = range(nsize)
    outwin = {}
    diff1 = nsize - set1[-1]
    diff2 = nsize - set2[-1]
    if diff1 < window/2:
        set1[-1] = set1[-1]+diff1
    else:
        set1 = num.hstack((set1,nsize))
    if diff2 < window/2:
        set2[-1] = set2[-1]+diff2
    else:
        set2 = num.hstack((set2,nsize))

    if len(set1) == 1:
        set1 = num.hstack((0,set1))
    if len(set2) == 1:
        set2 = num.hstack((0,set2))

    outSet1 = makeRangePairs(set1)
    outSet2 = makeRangePairs(set2)
    idict = {1:outSet1,2:outSet2}
 
    return idict

def makeRangePairs(indexList):
    """
    Given a list of indicies make a list of pairs.
    """
    
    set = []
    for i in range(len(indexList)-1):
        set.append([indexList[i],indexList[i+1]])

    return set

def movingMedian(data,window):
    """
    Returns the median smoothed version of a given function.
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
                   
def bin2(data,binlist):
	""" bins data and returns dictionary with an x and y list
	    for every binsize """
	
	npts = len(data['x'])
	#lst = [(bininc*i) for i in range(binnum)]
	#lst[0] = lst[0] + 1
	lst = binlist

	bind = {}

	for nbr in lst: 
		
		xlist = []
		ylist = []
		
		for i in range(npts-nbr):
			
			xlist.append(data['x'][i:i+nbr].mean() )
			ylist.append(data['y'][i:i+nbr].mean() )
		
		xlist = num.array(xlist)
        	ylist = num.array(ylist)
        	ylist /= ylist.mean()
	
		bind['binsize' + str(nbr)] = {'x':xlist,'y':ylist}
	
	return bind
