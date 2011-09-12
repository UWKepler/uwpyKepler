import numpy as num

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
