import numpy as num
from func import *
from lcmod import ApplyMaskPortions
import scipy
import warnings
warnings.simplefilter('ignore', num.RankWarning)

def FlagKeplerEvents(lcData):
    """ Flag the Kepler event flags as outliers. """
    
    # Creating the UnMaskedArray
    mask0 = num.ma.getmaskarray(lcData['x'])
    lcData['NoMask'] = mask0
    
    # Regions around discontinuities noted by Kepler
    # are marked with False (-99) values
    # will be detected later as marker of different portions
    id1 = num.where( (lcData['qflag'] == 1) | \
    (lcData['qflag'] == 1024) )[0]
    
    for el in id1:
        for j in [-1,0,1]:
            lcData['y'][el+j] = -99
            lcData['qflag'][el+j] = 0
     
    # Flagging events detected by Kepler as outliers
    id1 = num.where( (lcData['qflag'] == 64) | \
    (lcData['qflag'] == 128))[0]
    
    lcData['x'][id1] = num.ma.masked
    mask1=num.ma.copy(lcData['x'].mask)
    lcData['KEMask'] = mask1
    #print len(num.where(mask1))
    lcData['x'].mask = lcData['NoMask']
    
    return lcData

def RemoveBadEvents(lcData):
    """ Removes bad events from Kepler lightcurve """
    
    id1 = num.where( (lcData['y'] > 0) & \
    ( (lcData['qflag'] == 64) | (lcData['qflag'] == 128) | \
    (lcData['qflag'] == 0) ) )[0]
    #num.where( (lcData['qflag'] != 64) & \
    #(lcData['qflag'] != 128) & (lcData['qflag'] != 0))[0]
    lcData['x'] = lcData['x'][id1]
    lcData['y'] = lcData['y'][id1]
    lcData['yerr'] = lcData['yerr'][id1]
    lcData['KEMask'] = lcData['KEMask'][id1]
    lcData['NoMask'] = lcData['NoMask'][id1]
    lcData['cadence'] = lcData['cadence'][id1]
    del lcData['qflag']
    
    return lcData
    
def FlagEclipses(lcData,eclipseData,BJDREFI):
    """ This function flags points within a tranit and
    applies a mask.
        Input = data dictionary
              = eclipseData, a list of period, dur, and t0 for eclipses and transits
        Output = data dictionary with addition of the keys
        'TransitMask' and 'NoMask'
    """
    
    i = 0
    # use eclipse data to create the transit mask
    if eclipseData['eDataExists']:
        for koi in eclipseData['KOI'].keys():
            period = eclipseData['KOI'][koi]['Period']
            t0 = eclipseData['KOI'][koi]['T0']
            dur = eclipseData['KOI'][koi]['Duration']
            dur = (1.2*dur/24e0)
            t0 = t0 + 2454900e0
            width = dur/period
            maxphase=1-width/2
            minphase=width/2
            phase= (lcData['x']+BJDREFI-t0)/period\
            -(lcData['x']+BJDREFI-t0)//period
            idx=num.where((phase>maxphase)|(phase<minphase))[0]
            lcData['x'][idx]= num.ma.masked
            mask1=num.ma.copy(lcData['x'].mask)
            if i == 0:
                lcData['eMask']=mask1
            else:
                lcData['eMask']=\
                num.ma.mask_or(mask1,lcData['eMask'])
            i +=1
    else:
        lcData['eMask'] = lcData['NoMask']

    return lcData

def SplitPortions(lcData,gapsize):
    """
    This function finds gaps and splits data into portions.
    
    Input =  data - data dictionary
                gapsize - size of gap in days
    
    Output = new data dictionary with data split into portions
    """
    
    # defining new empty lists and stuff
    pcount=0
    istamps=[]
    outData = {}
    
    lcData['x'].mask = lcData['NoMask']
    lcData['y'].mask = lcData['NoMask']
    lcData['yerr'].mask = lcData['NoMask']
    
    # The grand master loop >=}
    # to make portion slices
    for i in range(len(lcData['x'])-1):
        diff =  lcData['cadence'][i+1]-lcData['cadence'][i]
        if pcount == 0:
            i0 = 0
        if pcount > 0:
            i0 = i1+1
        if diff > gapsize:
            i1 = i
            istamps.append([i0,i1])
            pcount += 1
    i1 = i+1
    istamps.append([i0,i1])
    # Applying slices
    for j in range(len(istamps)):
        #print istamps[j][0], istamps[j][1]
        outData['portion' + str(j+1)] =\
        {'kid':lcData['kid'],\
        'x':lcData['x'][istamps[j][0]:istamps[j][1]+1],\
        'y':lcData['y'][istamps[j][0]:istamps[j][1]+1],\
        'yerr':lcData['yerr'][istamps[j][0]:istamps[j][1]+1],\
        'NoMask':lcData['NoMask'][istamps[j][0]:istamps[j][1]+1],\
        'eMask':lcData['eMask'][istamps[j][0]:istamps[j][1]+1],\
        'KEMask':lcData['KEMask'][istamps[j][0]:istamps[j][1]+1],\
        'cadence':lcData['cadence'][istamps[j][0]:istamps[j][1]+1]}

    return outData
    
def FlagOutliers(lcData,medwin,threshold):
    """ This function flags outliers. 
        Inputs - data = data dictionary
               - medwin = the window size used to compute the median
               - threshold = the sigma-clipping factor (suggested, 3 or greater)
        Outputs - the data dictionary now contains mask arrays named
                'OutlierMask' and 'OTMask'
    """

    for portion in lcData.keys():
        npts = len(lcData[portion]['x'])
        # defining the window
        medflux = []
        medhalf = (medwin-1)/2
        # placing the window and computing the median
        for i in range(npts):
            i1 = max(0,i-medhalf)
            i2 = min(npts, i + medhalf)
            try:
                if num.ma.median(lcData[portion]['y'][i1:i2]).mask:
                    medflux.append(medflux[-1])
            except:
                medflux.append(num.ma.median(lcData[portion]['y'][i1:i2]))
        # finding outliers
        medflux = num.array(medflux)
        outliers = num.ma.getdata(lcData[portion]['y']) - medflux
        sigma = compute1Sigma(outliers)
        outliers = lcData[portion]['y'] - medflux
        idx=num.where( (abs(num.array(outliers))>threshold*sigma))[0]
        # creating the outlier mask
        lcData[portion]['x'].mask = lcData[portion]['NoMask']
        lcData[portion]['x'][idx] = num.ma.masked
        mask2 = num.ma.copy(lcData[portion]['x'].mask)
        lcData[portion]['OMask']=mask2
        mask3 = num.ma.mask_or(mask2,lcData[portion]['eMask'])
        mask4 = num.ma.mask_or(mask3,lcData[portion]['KEMask'])
        mask5 = num.ma.mask_or(mask2,lcData[portion]['KEMask'])
        lcData[portion]['ALLMask'] = mask4
        lcData[portion]['OKMask'] = mask5
        #print len(mask2), len(mask3), len(mask4), len(mask5)
        #if len(mask2) == 103:
            #print mask5
            #print len(mask2), len(lcData[portion]['KEMask'])
    
    return lcData
                    
def DetrendData(lcData, window, polyorder):
    """Detrends the data"""
    
    for portion in lcData.keys():
	lcData = ApplyMaskPortions(lcData,'ALLMask',portion)
        nsize = len(lcData[portion]['x'])
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
            
            unmasked = \
            num.where(lcData[portion]['x'][i1:i2].mask == False)[0]
            xdata =\
            num.array(lcData[portion]['x'][i1:i2][unmasked])
            ydata =\
            num.array(lcData[portion]['y'][i1:i2][unmasked])
            
            # find the fit
            coeff = scipy.polyfit(xdata, ydata, polyorder)

            # unmask data and apply the polynomial
            all_data_x =\
            num.ma.getdata(lcData[portion]['x'][i1:i2])
            
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
        newarr = num.ma.getdata(lcData[portion]['y'])/mergedy
        newerr = num.ma.getdata(lcData[portion]['yerr'])/mergedy
        
        lcData = ApplyMaskPortions(lcData,'NoMask',portion)
        lcData[portion]['correction'] = mergedy
        lcData[portion]['ydt'] = newarr
        lcData[portion]['yerrdt'] = newerr
        
    return lcData

def StackPortions(lcData):
    """rejoins/stacks all portions in the dictionary into one."""
    
    xarr=num.array([])
    yarr=num.array([])
    yerrarr=num.array([])
    ydtarr=num.array([])
    yerrdtarr=num.array([])
    corrarr = num.array([])
    cadarr = num.array([])
    eMask=num.array([])
    OMask=num.array([])
    NoMask=num.array([])
    KEMask=num.array([])
    OKMask=num.array([])
    ALLMask=num.array([])
    oData = {}
        
    for portion in lcData.keys():
        xarr=num.hstack((xarr,lcData[portion]['x']))
        yarr=num.hstack((yarr,lcData[portion]['y']))
        corrarr=num.hstack((corrarr,lcData[portion]['correction']))
        cadarr=num.hstack((cadrarr,lcData[portion]['cadence']))
        yerrarr=num.hstack((yerrarr,lcData[portion]['yerr']))
        ydtarr=num.hstack((ydtarr,lcData[portion]['ydt']))
        yerrdtarr=num.hstack((yerrdtarr,lcData[portion]['yerrdt']))
        OMask=num.hstack((OMask,lcData[portion]['OMask']))
        eMask=num.hstack((eMask,lcData[portion]['eMask']))
        NoMask=num.hstack((NoMask,lcData[portion]['NoMask']))
        KEMask=num.hstack((KEMask,lcData[portion]['KEMask']))
        OKMask=num.hstack((OKMask,lcData[portion]['OKMask']))
        ALLMask=num.hstack((ALLMask,lcData[portion]['ALLMask']))
    
    kid=lcData[portion]['kid']
    
    oData={'kid':kid,'x':xarr,'y':yarr,'yerr':yerrarr, \
    'ydt':ydtarr,'yerrdt':yerrdtarr,'correction':corrarr,\
    'cadence':cadarr,'eMask':eMask,'OMask':OMask,\
    'KEMask':KEMask,'OKMask':OKMask,'ALLMask':ALLMask,\
    'NoMask':NoMask}
    
    return oData