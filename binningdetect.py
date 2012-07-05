import numpy as num
import pdb

def binFromList(data,binlist):
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
			
			xlist.append(data['x'][i:i+nbr].mean())
			ylist.append(data['ydt'][i:i+nbr].mean())
		
		xlist = num.array(xlist)
        	ylist = num.array(ylist)
        	ylist /= ylist.mean()
	
		bind['binsize' + str(nbr)] = {'x':xlist,'y':ylist}
	
	return bind

# bins data given binsize in cadence
def bin(x, y, binsize):
    npts = len(x)
    nbins = npts / int(binsize)
    bx = num.zeros(nbins)
    by = num.zeros(nbins)
    for n in range(nbins):
        i1 = n * binsize
        i2 = (n + 1) * binsize
        bx[n] = num.mean(x[i1:i2])
        by[n] = num.mean(y[i1:i2])
    if npts % binsize > 0:
        bx = num.hstack( (bx, num.zeros(1)) )
        by = num.hstack( (by, num.zeros(1)) )
        bx[-1] = num.mean(x[nbins:-1])
        by[-1] = num.mean(y[nbins:-1])
    
    return bx, by

