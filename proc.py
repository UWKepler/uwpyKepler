import sys
import numpy as num
import scipy
import pylab
#num.warning.warn(action = 'ignore')
import warnings
warnings.simplefilter('ignore', num.RankWarning)
from uwpyKepler.io import ApplyMask
 
def detrendData(data, window, polyorder):

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
        