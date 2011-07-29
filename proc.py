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
 
def detrendData(data, window, polyorder):
    """Detrends the data"""
    dout = {}
    # loop through portions
    
    for portion in data.keys():
	data = ApplyMask(data,'OTMask')
        nsize = len(data[portion]['x'])
        nfullwindows = int(num.floor(nsize/window))
        leftover = nsize - nfullwindows*window
        i1 = 0
        i2 = window+1
        newarr = num.ma.masked_array([])
        newx = num.ma.masked_array([])
        newerr = num.ma.masked_array([])
	correction = num.array([])
        if nfullwindows == 0:
            i1 = 0
            i2 = nsize
            nfullwindows = 1
        for i in range(nfullwindows):
            if i == nfullwindows-1:
                i2 = i2 + leftover
            xdata = data[portion]['x'][i1:i2][num.where(data[portion]['x'][i1:i2].mask == False)]
            ydata = data[portion]['y'][i1:i2][num.where(data[portion]['y'][i1:i2].mask == False)]
            
            # find the fit
            coeff = scipy.polyfit(xdata,ydata, polyorder)

            # unmask data and apply the polynomial
            data[portion]['x'][i1:i2].mask = data[portion]['UnMasked'][i1:i2]
            data[portion]['y'][i1:i2].mask = data[portion]['UnMasked'][i1:i2]
            data[portion]['yerr'][i1:i2].mask = data[portion]['UnMasked'][i1:i2]
            
            outy = scipy.polyval(coeff,data[portion]['x'][i1:i2])
            
            # apply correction
            d1 = data[portion]['y'][i1:i2]/outy
            d2 = data[portion]['yerr'][i1:i2]/outy
            
            newarr = num.ma.hstack( (newarr,d1))
            newerr = num.ma.hstack( (newerr,d2))
	    correction = num.hstack((correction,outy))
            i1 = i2
            i2 = i2+window
        
        data = ApplyMask(data,'UnMasked')
        dout[portion] = {'x':data[portion]['x'],'y':newarr,'yerr':newerr,'TransitMask':data[portion]['TransitMask'],'OTMask':data[portion]['OTMask'],'OutlierMask':data[portion]['OutlierMask'],'UnMasked':data[portion]['UnMasked'],'Correction':correction}
        
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
    
    print len(data['x']), len(data['OutlierMask']), len(data['TransitMask']),len(data['UnMasked']), len(data['OTMask'])
    
    xnew = []
    ynew = []
    yerrnew = []

    for el in idx:
        xnew.append(data['x'][el])
        ynew.append(data['y'][el])
        yerrnew.append(data['yerr'][el])

    print num.shape(num.array(xnew).ravel())
    
    dout= {'kid':data['kid'],'x':num.array(xnew).ravel(),'y':num.array(ynew).ravel(),'yerr':num.array(yerrnew).ravel()}
        
    return dout

def cutTransits(pd):
    	""" This function cuts out points within a tranit.
         
         Input = data dictionary
         Output = data dictionary without points in transits.
        """

        # reading planetary data from database
        db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
        cursor = db.cursor()
        foo1    = 'select Period, Epoch, Dur from KEPPC where (KID = %s)' % (pd['kid'])
        cursor.execute(foo1)
        results = cursor.fetchall()
        period, t0, dur = results[0][0], results[0][1], results[0][2]
        dur = (1.2*dur/24e0)
        t0 = t0 + 54900e0
        # defining start and end time lists
        width = dur/period
        maxphase=1-width/2
        minphase=width/2
        phase= (pd['x']-t0)/period-(pd['x']-t0)//period
        idx=num.where((phase<maxphase)&(phase>minphase))
        #import pdb; pdb.set_trace()
        #mask0=num.ma.getmaskarray(pd['x'])
        
        xnew=pd['x'][idx]
        ynew=pd['y'][idx]
        
        pd['x']=xnew
        pd['y']=ynew
        #mask1=num.ma.copy(pd['x'].mask)
        #pd['TransitMask']=mask1
        #pd['UnMasked']=mask0

        return pd

def cutOT(data):
    d2=kep.proc.cutTransits(data)
    d3=kep.proc.cutOutliers(d2,10,4)
    return d3

def stackPortions(data):
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
        print type(data[portion]['TransitMask'])
        TransitMask=num.hstack((TransitMask,data[portion]['TransitMask']))
        OutlierMask=num.hstack((OutlierMask,data[portion]['OutlierMask']))
        OTMask=num.hstack((OTMask,data[portion]['OTMask']))
        UnMasked=num.hstack((UnMasked,data[portion]['UnMasked']))
        
        #print len(data[portion]['x']), len(xarr), portion
    kid=data[portion]['kid']
    
    
    pd={'OTMask':OTMask,'TransitMask':TransitMask,'OutlierMask':OutlierMask,'UnMasked':UnMasked,'yerr':yerrarr,'y':yarr,'x':xarr,'kid':kid}
    return pd