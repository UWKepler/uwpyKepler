import sys
import numpy as num
import scipy
import pylab

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

def detrendData(data, minwindow, maxwindow, polyorder):

    #dout = data.copy()
    print data.keys()
    for portion in data.keys():
        nsize = len(data[portion]['x'])
        nfullwindows = int(num.floor(nsize/maxwindow))
        leftover = nsize - nfullwindows*maxwindow
        i1 = 0
        i2 = maxwindow+1
        #print portion, nsize, nfullwindows
        if nfullwindows == 0:
            i1 = 0
            i2 = nsize
            nfullwindows = 1
        for i in range(nfullwindows):
            xdata = data[portion]['x'][i1:i2]
            ydata = data[portion]['y'][i1:i2]
            print len(xdata), len(ydata)
            coeff = scipy.polyfit(xdata,ydata, polyorder)
            outx = scipy.polyval(coeff,num.ma.getdata(data[portion]['x'][i1:i2]))
            print len(num.ma.getdata(data[portion]['y'][i1:i2])), type(num.ma.getdata(data[portion]['y'][i1:i2]))
            print len(outx), type(outx), len(data[portion]['y'][i1:i2])
            print data[portion]['y'][i1]/outx[0]
            data[portion]['y'][i1:i2] = data[portion]['y'][i1:i2]/outx 
            #pylab.plot(xdata,ydata,'bo')
            #pylab.plot(xdata,outx,'r.')
            #pylab.show()
            #dout[portion]['y'] = outx 
            i1 = i2
            i2 = i2+maxwindow
        #pylab.plot(num.ma.getdata(data[portion]['x']),num.ma.getdata(data[portion]['y']), 'g.')
    #pylab.show()
        #pylab.plot(d[portion]['x'],d[portion]['y'],'b.')
        #pylab.plot(d3[portion]['x'],out,'b.')
        #pylab.show()
        #dout[portion] = {'x':d3[portion]['x'],'y':out}
    return data
        