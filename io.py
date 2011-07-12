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
    time    = num.array([x[2] for x in results])
    corflux = num.array([x[7] for x in results])
    corerr  = num.array([x[8] for x in results])
    
    if len(time) == 0:
        print 'No data found for KID  %s' % (KeplerID)
        sys.exit(1)
    
    idx = []
    for i in range(len(time)):
        if corflux[i] != 0 and corerr[i] != 0: 
            idx.append(i)
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
        dur = (dur/24e0)/2e0
        t0 = t0 + 54900e0
        # defining start and end time lists
        transitslist = []
	transitelist = []
	
        tstime0 = t0 - dur
	tetime0 = t0 + dur
	
	n = 0
	tstime = tstime0 + n*period
        while pd['x'][0] <= tstime < pd['x'][-1]:
                tstime = tstime0 + n*period
                transitslist.append(tstime)
                tsmtime = tstime0 - n*period
                transitslist.append(tsmtime)
                n += 1
	
	m = 0
	tetime = tetime0 + m*period
        while pd['x'][0] <= tetime < pd['x'][-1]:
                tetime = tetime0 + m*period
                transitelist.append(tetime)
                temtime = tetime0 - m*period
                transitelist.append(temtime)
                m += 1
	
	tsnumber = len(transitslist)
	tenumber = len(transitelist)
	
	if tsnumber != tenumber:
                print 'Number of Start times does not equal number of End times'
                print tsnumber, tenumber
		sys.exit(1)

        TransitFlag = []
        for i in range(len(pd['x'])):
                withinTransit = False
                for k in range(len(transitslist)):
                        if pd['x'][i] > transitslist[k] and pd['x'][i] < transitelist[k]:
                                withinTransit = True
                TransitFlag.append(withinTransit)
        pd['TransitFlag'] = TransitFlag
            
        return pd

def SplitGap(data,gapsize):
	"""
        This function finds gaps and splits data into portions.
        The gapsize must be provided in days.
	"""
	
	#defining new empty lists and stuff
        DT = []
	pcount=0
	xarr=[]
	yarr=[]
	yerrarr=[]
	tstart=[]
	tend=[]
        transitflag = []
	pd={}
	
	# minimum sized gap that we are flagging, 2.4 hours
	# The grand master loop >=}
	for i in range(len(data['x'])- 1):
		dt =  data['x'][i + 1]- data['x'][i]
		if dt > gapsize:
		        # Collecting the time and magnitude 
                        # of gap every time we detect one.
			DT.append(dt)
			tstart.append(data['x'][i])
			tend.append(data['x'][i + 1])
			pcount += 1
			pd['portion' + str(pcount)] = {'kid':data['kid'],'x':xarr, 'y':yarr, 'yerr':yerrarr,'TransitFlag':transitflag}
			xarr=[]
			yarr=[]
			yerrarr=[]
                        transitflag=[]
		else:
			xarr.append(data['x'][i])
			yarr.append(data['y'][i])
			yerrarr.append(data['yerr'][i])
                        transitflag.append(data['TransitFlag'][i])
	
	pcount += 1
	pd['portion' + str(pcount)] = {'kid':data['kid'],'x':xarr, 'y':yarr, 'yerr':yerrarr,'TransitFlag':transitflag}
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
        
        check = abs(outliers)<threshold*sigma
        
        #outliers = outliers*(abs(outliers)<threshold*sigma)+medflux
        for el in check:
            if el == 1:
                OutlierFlag.append(False)
            if el == 0:
                OutlierFlag.append(True)
                
        dout[portion] = {'kid':data[portion]['kid'],'x':data[portion]['x'],'y':data[portion]['y'],'yerr':num.array(data[portion]['yerr']),'TransitFlag':data[portion]['TransitFlag'],'OutlierFlag':OutlierFlag}
        
    return dout
