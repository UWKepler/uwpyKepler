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
    data = ApplyMask(data,'OTMask')
    for portion in data.keys():
        nsize = len(data[portion]['x'])
        nfullwindows = int(num.floor(nsize/window))
        leftover = nsize - nfullwindows*window
        i1 = 0
        i2 = window+1
        #print portion, nsize, nfullwindows
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
            
            outx = scipy.polyval(coeff,data[portion]['x'][i1:i2])
         
            d1 = data[portion]['y'][i1:i2]/outx
            d2 = data[portion]['yerr'][i1:i2]/outx
            print num.ma.count(xdata), num.ma.count(ydata), num.ma.count(d1), num.ma.count(d2)
            
            newarr = num.ma.hstack( (newarr,d1))
            newerr = num.ma.hstack( (newerr,d2))
	    correction = num.hstack((correction,outx))
            i1 = i2
            i2 = i2+window
        
        data = ApplyMask(data,'UnMasked')
        dout[portion] = {'x':data[portion]['x'],'y':newarr,'yerr':newerr,'TransitMask':data[portion]['TransitMask'],'OTMask':data[portion]['OTMask'],'OutlierMask':data[portion]['OutlierMask'],'UnMasked':data[portion]['UnMasked'],'Correction':correction}
        
    return dout

def cutOutliers(data,medwin,threshold):
    """ This function cuts out outliers. 
        Inputs - data = data dictionary
               - medwin = the window size used to compute the median
               - threshold = the sigma-clipping factor (suggested, 3 or greater)
        Outputs - the x data now only contains times that don't correspond to outliers. 
    """
    
    dout = {}
    # cycling through portions
    #for portion in data.keys():
        #data[portion]['x'].mask = data[portion]['UnMasked']
        #data[portion]['y'].mask = data[portion]['UnMasked']
        #data[portion]['yerr'].mask = data[portion]['UnMasked']
    npts = len(data['x'])
    
    # defining the window
    medflux = []
    medhalf = (medwin-1)/2

    # placing the window and computing the median
    for i in range(npts):
        i1 = max(0,i-medhalf)
        i2 = min(npts, i + medhalf)
        medflux.append(num.median(data['y'][i1:i2]))
    
    # finding outliers
    medflux = num.array(medflux)
    outliers = data['y'] - medflux
    
    outliers.sort()
    sigma = (outliers[.8415*npts]-outliers[.1585*npts])/2
    outliers = data['y'] - medflux
    
    # tagging outliers
    idx=num.where( (abs(num.array(outliers))<threshold*sigma) )
    
    xnew=data['x'][idx]
    ynew=data['y'][idx]


    data['x']=xnew
    data['y']=ynew
    
    #creating the outlier mask
    #data[portion]['x'].mask = data[portion]['UnMasked']
    #data[portion]['x'][idx[0]] = num.ma.masked
    
    #mask2 = num.ma.copy(data[portion]['x'].mask)
    
    #data[portion]['OutlierMask']=mask2
    
    # creating the outlier + transit mask
    #mask3 = num.ma.mask_or(data[portion]['TransitMask'],mask2)
    
    dout= {'kid':data['kid'],'x':data['x'],'y':data['y'],'yerr':data['yerr']}
        
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

def cutAll(data):
    d2=kep.proc.cutTransits(data)
    d3=kep.proc.cutOutliers(d2,10,4)
    return d3

def stackPortions(data):
    """rejoins/stacks all portions in the dictionary into one."""
    xarr=num.ma.array([])
    yarr=num.ma.array([])
    yerrarr=num.ma.array([])
    for portion in data.keys():
        xarr=num.ma.hstack((xarr,data[portion]['x']))
        yarr=num.ma.hstack((yarr,data[portion]['y']))
        yerrarr=num.ma.hstack((yerrarr,data[portion]['yerr']))
        #print len(data[portion]['x']), len(xarr), portion
    
    #xarr=num.ma.hstack([(data[portion]['x'], data[portion+1]['x']) for i in range(len(data)-1)])
    #print len(data[portion]['x']), portion
    #print 'here', len(xarr)
    pd={'x':xarr,'y':yarr,'yerr':yerrarr}
    return pd