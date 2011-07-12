
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

def detrendData(data, windowsize, polyorder):

    d = removeOutliers(removeTransits(data3))
    dout = {}
    for portion in d.keys():
        if windowsize < len(d[portion]['x'])/2e0:
            i1 = max(0,i-medhalf)
            i2 = min(npts, i + medhalf)
            xdata = d[portion]['x'][i1:i2]
            ydata = d[portion]['y'][i1:i2]
            coeff = scipy.polyfit(xdata,ydata, polyorder)
            outx = scipy.polyval(coeff,data[portion]['x'])
        import pylab
        pylab.plot(d[portion]['x'],d[portion]['y'],'b.')
        pylab.plot(d3[portion]['x'],out,'b.')
        pylab.show()
        dout[portion] = {'x':d3[portion]['x'],'y':out}
        return 