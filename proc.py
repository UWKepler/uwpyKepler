import sys
import MySQLdb
import numpy as num
import scipy
import pylab
import uwpyKepler as kep
#num.warning.warn(action = 'ignore')
import warnings
warnings.simplefilter('ignore', num.RankWarning)
from uwpyKepler.io import ApplyMask
 
def detrendData(data, eData, window, polyorder):
    """Detrends the data"""
    
    dout = {}
    # loop through portions
    
    for portion in data.keys():
        
	data = ApplyMask(data,'OTMask')
        
        nsize = len(data[portion]['x'])
        points = num.arange(0,nsize,0.5*window)
        points= num.hstack( (points,nsize))
        dtfunc1 = num.array([])
        dtfunc2 = num.array([])
        set1 = []
        set2 = []
        weight1 = (num.cos(num.pi*num.arange(nsize)/window))**2
        weight2 = (num.sin(num.pi*num.arange(nsize)/window))**2
        
        for i in range(len(points)):
            if i < len(points)-1:
                i1  = long(max(0,points[i]-window/2))
            else:
                i1 = long(points[i-2]+window/2)
            i2 = long(min(points[i]+window/2,points[-1]))
            
            unmasked = num.where(data[portion]['x'][i1:i2].mask == False)[0].tolist()
            xdata = num.array(data[portion]['x'][i1:i2][unmasked])
            ydata = num.array(data[portion]['y'][i1:i2][unmasked])
            
            # find the fit
            coeff = scipy.polyfit(xdata, ydata, polyorder)

            # unmask data and apply the polynomial
            all_data_x = num.ma.getdata(data[portion]['x'][i1:i2])
            
            outy = scipy.polyval(coeff,all_data_x)

            if i%2 == 0:
                set1.append( (i1,i2) )
                dtfunc1 = num.hstack((dtfunc1,outy))
                #pylab.plot(all_data_x,outy,'r-',linewidth=3)
            else:
                set2.append( (i1,i2) )
                dtfunc2 = num.hstack((dtfunc2,outy))
                #pylab.plot(all_data_x,outy,'c-',linewidth=3)
            
        mergedy = weight1*dtfunc1 + weight2*dtfunc2
        #pylab.plot(data[portion]['x'],mergedy,'k-')

        # apply correction
        newarr = num.ma.getdata(data[portion]['y'])/mergedy
        newerr = num.ma.getdata(data[portion]['yerr'])/mergedy
        
        data = ApplyMask(data,'UnMasked')
        if eData['bool']==False:
            dout[portion] = {'kid':data[portion]['kid'],'x':data[portion]['x'],'y':newarr,'yerr':newerr,'OTMask':data[portion]['OTMask'],'OutlierMask':data[portion]['OutlierMask'],'UnMasked':data[portion]['UnMasked'],'Correction':mergedy}
        else:
            dout[portion] = {'kid':data[portion]['kid'],'x':data[portion]['x'],'y':newarr,'yerr':newerr,'TransitMask':data[portion]['TransitMask'],'OTMask':data[portion]['OTMask'],'OutlierMask':data[portion]['OutlierMask'],'UnMasked':data[portion]['UnMasked'],'Correction':mergedy}
        
        #pylab.plot(data[portion]['x'],(weight1*1.6e4) + num.median(mergedy) + 20000,'r-')
        #pylab.plot(data[portion]['x'],(weight2*1.6e4) + num.median(mergedy) + 20000,'c-')
    return dout

def cutOutliers(data):
    """ This function cuts out outliers. 
        Inputs - data = data dictionary
               - medwin = the window size used to compute the median
               - threshold = the sigma-clipping factor (suggested, 3 or greater)
        Outputs - the x data now only contains times that don't correspond to outliers. 
    """

    # tagging outliers
    idx=num.where(data['OutlierMask']==False)
    
    xnew = []
    ynew = []
    yerrnew = []

    for el in idx:
        xnew.append(data['x'][el])
        ynew.append(data['y'][el])
        yerrnew.append(data['yerr'][el])
    
    dout= {'kid':data['kid'],'x':num.array(xnew).ravel(),'y':num.array(ynew).ravel(),'yerr':num.array(yerrnew).ravel()}
        
    return dout

def cutTransits(dTransit):
    	""" This function cuts out points within a tranit.
         
         Input = data dictionary
         Output = data dictionary without points in transits.
        """
    	
	xnew = []
	ynew = []
	yerrnew = []
	
	idx=num.where(dTransit['TransitMask'] == False)
        for element in idx:
		xnew.append(dTransit['x'][element])
		ynew.append(dTransit['y'][element])
		yerrnew.append(dTransit['yerr'][element])
		
	dout= {'kid':dTransit['kid'],'x':num.array(xnew).ravel(),'y':num.array(ynew).ravel(),'yerr':num.array(yerrnew).ravel()}

        return dout

def cutOT(data):
	""" This function cuts outliers and transits from the data and returns x, y and yerr"""

	idx = num.where((data['OutlierMask']==False) & (data['TransitMask']==False))
	
	idx2 = num.where(data['TransitMask'] != False)
	idx3 = num.where(data['OutlierMask'] != False)
	
	#print len(idx2[0])
	#print len(idx3[0])
	#print len(idx[0])
		
	xnew = []
	ynew = []
	yerrnew = []
	
	for el in idx[0]:
		xnew.append(data['x'][el])
		ynew.append(data['y'][el])
		yerrnew.append(data['yerr'][el])
		#print el
	
	dout= {'kid':data['kid'],'x':num.array(xnew).ravel(),'y':num.array(ynew).ravel(),'yerr':num.array(yerrnew).ravel()}
	
	return dout

def onlyOutliers(data):
    """ This function singles out outliers.
        Inputs - data = data dictionary
               - medwin = the window size used to compute the median
               - threshold = the sigma-clipping factor (suggested, 3 or greater)
        Outputs - the x data now only contains times that correspond to outliers.
    """

    # tagging outliers

    idx=num.where(data['OutlierMask']==True)
      
    xnew = []
    ynew = []
    yerrnew = []

    for el in idx:
        xnew.append(data['x'][el])
        ynew.append(data['y'][el])
        yerrnew.append(data['yerr'][el])
   
    dout= {'kid':data['kid'],'x':num.array(xnew).ravel(),'y':num.array(ynew).ravel(),'yerr':num.array(yerrnew).ravel()}
       
    return dout

def onlyTransits(dTransit):
        """ This function singles out points within a tranit.
        
         Input = data dictionary
         Output = data dictionary without points in transits.
        """

	xnew = []
	ynew = []
	yerrnew = []
	
	idx=num.where(dTransit['TransitMask'] == True)
        for element in idx:
		xnew.append(dTransit['x'][element])
		ynew.append(dTransit['y'][element])
		yerrnew.append(dTransit['yerr'][element])
		
		
	dout= {'kid':dTransit['kid'],'x':num.array(xnew),'y':num.array(ynew),'yerr':num.array(yerrnew)}

        return dout
	
def stackPortions(data, eData):
    
    if eData['bool']==False:
        """rejoins/stacks all portions in the dictionary into one."""
        xarr=num.array([])
        yarr=num.array([])
        yerrarr=num.array([])
        TransitMask=num.array([])
        OutlierMask=num.array([])
        OTMask=num.array([])
        UnMasked=num.array([])

        
        for portion in data.keys():
            xarr=num.hstack((xarr,data[portion]['x']))
            yarr=num.hstack((yarr,data[portion]['y']))
            yerrarr=num.hstack((yerrarr,data[portion]['yerr']))
            OutlierMask=num.hstack((OutlierMask,data[portion]['OutlierMask']))
            OTMask=num.hstack((OTMask,data[portion]['OTMask']))
            UnMasked=num.hstack((UnMasked,data[portion]['UnMasked']))
            
            #print len(data[portion]['x']), len(xarr), portion
            kid=data[portion]['kid']
        #kid=kid1
        pd={'OTMask':OTMask,'OutlierMask':OutlierMask,'UnMasked':UnMasked,'yerr':yerrarr,'y':yarr,'x':xarr,'kid':kid}
    else:
        """rejoins/stacks all portions in the dictionary into one."""
        xarr=num.array([])
        yarr=num.array([])
        yerrarr=num.array([])
        TransitMask=num.array([])
        OutlierMask=num.array([])
        OTMask=num.array([])
        UnMasked=num.array([])
        print data.keys()
        
        for portion in data.keys():
            xarr=num.hstack((xarr,data[portion]['x']))
            yarr=num.hstack((yarr,data[portion]['y']))
            yerrarr=num.hstack((yerrarr,data[portion]['yerr']))
            TransitMask=num.hstack((TransitMask,data[portion]['TransitMask']))
            OutlierMask=num.hstack((OutlierMask,data[portion]['OutlierMask']))
            OTMask=num.hstack((OTMask,data[portion]['OTMask']))
            UnMasked=num.hstack((UnMasked,data[portion]['UnMasked']))
            
            #print len(data[portion]['x']), len(xarr), portion
            kid=data[portion]['kid']
        #kid=kid1
    
        pd={'OTMask':OTMask,'TransitMask':TransitMask,'OutlierMask':OutlierMask,'UnMasked':UnMasked,'yerr':yerrarr,'y':yarr,'x':xarr,'kid':kid}
    return pd

def bin(data):
    """bins data, doesn't return anything, just graphs"""
    npts = len(data['x'])
    for nperbin in [1, 5, 10, 20, 50, 100]:
    
        binX    = []
        binY    = []
        for i in range(npts-nperbin):
            binX.append( data['x'][i:i+nperbin].mean() )
            binY.append( data['y'][i:i+nperbin].mean() )
        binX    = num.array(binX)
        binY    = num.array(binY)
        binY   /= binY.mean()
    
        if nperbin == 1: symb='ro'
        if nperbin == 5: symb='mo'
        if nperbin == 10: symb='yo'
        if nperbin == 20: symb='go'
        if nperbin == 50: symb='co'
        if nperbin == 100: symb='bo'
        pylab.errorbar(binX, binY, fmt=symb)
        pylab.ylim()
        pylab.title('binsize = %d' % (nperbin))
        
    pylab.show()

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
