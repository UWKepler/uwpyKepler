import sys
import numpy as num
import scipy
import pylab
#num.warning.warn(action = 'ignore')
import warnings
warnings.simplefilter('ignore', num.RankWarning)

def removeOutliers(d):
    """ returns array without outliers """
    
    dout = {}
    for key in d.keys():
        x = []
        y = []
        yerr = []
        outliers = []
        transitflag = []
        for i in range(len(d[key]['OutlierFlag'])):
            if not d[key]['OutlierFlag'][i]:
                x.append(d[key]['x'][i])
                y.append(d[key]['y'][i])
                yerr.append(d[key]['yerr'][i])
                outliers.append(d[key]['OutlierFlag'][i])
                transitflag.append(d[key]['TransitFlag'][i])
        dout[key] = {'x':x,'y':y,'yerr':yerr,'kid':d[key]['kid'],'OutlierFlag':outliers,'TransitFlag':transitflag}
        
    return dout

def removeTransits(d):
    """ excludes Transit data """

    dout = {}
    for key in d.keys():
        x = []
        y = []
        yerr = []
        outliers = []
        transitflag = []
        for i in range(len(d[key]['TransitFlag'])):
            if not d[key]['TransitFlag'][i]:
                x.append(d[key]['x'][i])
                y.append(d[key]['y'][i])
                yerr.append(d[key]['yerr'][i])
                outliers.append(d[key]['OutlierFlag'][i])
                transitflag.append(d[key]['TransitFlag'][i])
        dout[key] = {'x':x,'y':y,'yerr':yerr,'kid':d[key]['kid'],'OutlierFlag':outliers,'TransitFlag':transitflag}

    return dout

def detrendData(data, window, polyorder):

    dout = {}
    #print data.keys()
    for portion in data.keys():
        data[portion]['x'].mask = data[portion]['MaskBoth']
        data[portion]['y'].mask = data[portion]['MaskBoth']
        data[portion]['yerr'].mask = data[portion]['MaskBoth']
        nsize = len(data[portion]['x'])
        nfullwindows = int(num.floor(nsize/window))
        leftover = nsize - nfullwindows*window
        i1 = 0
        i2 = window+1
        #print portion, nsize, nfullwindows
        newarr = num.ma.masked_array([])
        newx = num.ma.masked_array([])
        newerr = num.ma.masked_array([])
        if nfullwindows == 0:
            i1 = 0
            i2 = nsize
            nfullwindows = 1
        for i in range(nfullwindows):
            if i == nfullwindows-1:
                i2 = i2 + leftover
            xdata = num.ma.copy(data[portion]['x'][i1:i2])
            ydata = num.ma.copy(data[portion]['y'][i1:i2])
            #print len(xdata), len(ydata), i1, i2, portion
            # find the fit
            coeff = scipy.polyfit(xdata,ydata, polyorder)
            
            #unmask data and apply the polynomial
            data[portion]['x'][i1:i2].mask = data[portion]['UnMasked'][i1:i2]
            data[portion]['y'][i1:i2].mask = data[portion]['UnMasked'][i1:i2]
            data[portion]['yerr'][i1:i2].mask = data[portion]['UnMasked'][i1:i2]
            
            outx = scipy.polyval(coeff,data[portion]['x'][i1:i2])
        
            #pylab.plot(xdata,ydata,'bo')
            #pylab.plot(data[portion]['x'][i1:i2],outx,'r-')
            #pylab.plot(data[portion]['x'][i1:i2],data[portion]['y'][i1:i2],'y.')
            
            d0 = data[portion]['x'][i1:i2]
            d1 = data[portion]['y'][i1:i2]/outx
            d2 = data[portion]['yerr'][i1:i2]/outx
            #print len(outx), len(data[portion]['y'][i1:i2]), len(newarr), len(d1), len(d2)
            #print num.shape(outx), num.shape(data[portion]['y'][i1:i2]), num.shape(newarr), num.shape(d1), num.shape(d2)
            newx = num.ma.hstack((newarr,d0))
            newarr = num.ma.hstack( (newarr,d1))
            newerr = num.ma.hstack( (newerr,d2))
            i1 = i2
            i2 = i2+window
        dout[portion] = {'x':data[portion]['x'],'y':newarr,'yerr':newerr,'TransitMask':data[portion]['TransitMask'],'MaskBoth':data[portion]['MaskBoth'],'OutlierMask':data[portion]['OutlierMask'],'UnMasked':data[portion]['UnMasked']}
        
    #pylab.show()
    return dout
        