import MySQLdb
import sys
import numpy as num
import scipy
    
def ReadLightCurve(KeplerID):
    """ This function reads the Kepler database and returns the 
    corrected lightcurve in a dictionary
    """

    db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
    cursor = db.cursor()
    foo    = 'select * from source where (KEPLERID = %s)' % (KeplerID)
    cursor.execute(foo)
    results = cursor.fetchall()
    
    # reading time, corrected flux and flux errors
    time    = num.ma.array([x[2] for x in results])
    corflux = num.ma.array([x[7] for x in results])
    corerr  = num.ma.array([x[8] for x in results])
    
    if len(time) == 0:
        print 'No data found for KID  %s' % (KeplerID)
        sys.exit(1)
    
    idx = num.where((corflux>0)&(corerr>0))
    
    time    = time[idx]
    corflux = corflux[idx]
    corerr  = corerr[idx]
    
    return {'kid':KeplerID,'x':time,'y':corflux,'yerr':corerr}

def FlagTransits(pd):
    	""" This function flags data arrays for known transits """

        # reading planetary data from database
        db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
        cursor = db.cursor()
        foo1    = 'select Period, Epoch, Dur from KEPPC where (KID = %s)' % (pd['kid'])
        cursor.execute(foo1)
        results = cursor.fetchall()
        period, t0, dur = results[0][0], results[0][1], results[0][2]
        dur = (dur/24e0)
        t0 = t0 + 54900e0
        # defining start and end time lists
        width = dur/period
        maxphase=1-width/2
        minphase=width/2
        phase= (pd['x']-t0)/period-(pd['x']-t0)//period
        idx=num.where((phase>maxphase)|(phase<minphase))
        #import pdb; pdb.set_trace()
        pd['x'][idx]= num.ma.masked
        pd['y'][idx]= num.ma.masked
        pd['yerr'][idx]= num.ma.masked

        return pd

def SplitGap(data,gapsize):
	"""
        This function finds gaps and splits data into portions.
        The gapsize must be provided in days.
	"""
	
	#defining new empty lists and stuff
	pcount=0
        istamps=[]
        #transitflag = []
	pd={}
	
	# minimum sized gap that we are flagging, 2.4 hours
	# The grand master loop >=}
	for i in range(len(data['x'])-1):
		dt =  data['x'][i+1]- data['x'][i]
                if pcount == 0:
                    i0 = 0
                if pcount > 0:
                    i0 = i1+1
		if dt > gapsize:
                    i1 = i
                    istamps.append([i0,i1])
                    pcount += 1
        i1 = i+1
        istamps.append([i0,i1])
        for j in range(len(istamps)):
            pd['portion' + str(j+1)] = {'kid':data['kid'],'x':data['x'][istamps[j][0]:istamps[j][1]+1], 'y':data['y'][istamps[j][0]:istamps[j][1]+1], 'yerr':data['yerr'][istamps[j][0]:istamps[j][1]+1]}
             
        return pd
                
def FlagOutliers(data,medwin,threshold):
    """ This function flags outliers. 
        medwin is the window size used to compute the median
        threshold is the sigma-clipping factor (suggested, 3 or greater)
    """
    
    dout = {}
    for portion in data.keys():
        
        npts = len(data[portion]['x'])
        medflux = []
        medhalf = (medwin-1)/2
        OutlierFlag = []
        
        for i in range(npts):
            i1 = max(0,i-medhalf)
            i2 = min(npts, i + medhalf)
            medflux.append(num.median(data[portion]['y'][i1:i2]))
        
        medflux = num.array(medflux)
        outliers = data[portion]['y'] - medflux
        
        outliers.sort()
        sigma = (outliers[.8415*npts]-outliers[.1585*npts])/2
        outliers = data[portion]['y'] - medflux
        
        #check = abs(outliers)<threshold*sigma

        idx=num.where(abs(num.array(outliers))>threshold*sigma)
        #import pdb; pdb.set_trace()
        data[portion]['x'][idx] = num.ma.masked
        data[portion]['y'][idx] = num.ma.masked
        data[portion]['yerr'][idx] = num.ma.masked
                
        dout[portion] = {'kid':data[portion]['kid'],'x':data[portion]['x'],'y':data[portion]['y'],'yerr':num.array(data[portion]['yerr'])}
        
    return dout
