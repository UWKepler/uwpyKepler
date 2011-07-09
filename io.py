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
    
    time    = num.array([x[2] for x in results])
    #rawflux = num.array([x[5] for x in results])
    #rawerr  = num.array([x[6] for x in results])
    corflux = num.array([x[7] for x in results])
    corerr  = num.array([x[8] for x in results])
    
    if len(time) == 0:
        print 'No data on object %s' % (KeplerID)
        sys.exit(1)
    
    idx = []
    for i in range(len(time)):
        if corflux[i] != 0 and corerr[i] != 0: 
            idx.append(i)
    time    = time[idx]
    corflux = corflux[idx]
    corerr  = corerr[idx]
    
    return {'kid':KeplerID,'x':time,'y':corflux,'yerr':corerr}

def splitgap(d1,gapsize):
	"""
        This function finds gaps and splits data into portions.
	"""
	
	#defining new empty lists and stuff
        DT = []
	pcount=0
	xarr=[]
	yarr=[]
	yerrarr=[]
        transitflag = []
	tstart=[]
	tend=[]
	pd={}
	
	# minimum sized gap that we are flagging, 2.4 hours
	# The grand master loop >=}
	for i in range(len(d1['x'])- 1):
		dt =  d1['x'][i + 1]- d1['x'][i]
		if dt > gapsize:
		# Collecting the time and magnitude of gap every time we detect one.
			DT.append(dt)
			tstart.append(d1['x'][i])
			tend.append(d1['x'][i + 1])
			pcount += 1
			pd['portion' + str(pcount)] = {'kid':d1['kid'],'x':xarr, 'y':yarr, 'yerr':yerrarr,'TransitFlag':transitflag}
			xarr=[]
			yarr=[]
			yerrarr=[]
                        transitflag=[]
		else:
			xarr.append(d1['x'][i])
			yarr.append(d1['y'][i])
			yerrarr.append(d1['yerr'][i])
                        transitflag.append(d1['TransitFlag'][i])
	
	pcount += 1
	pd['portion' + str(pcount)] = {'kid':d1['kid'],'x':xarr, 'y':yarr, 'yerr':yerrarr,'TransitFlag':transitflag}
        #pd['kid'] = 
	return pd
	#return {'kid':d1['kid'],'x':xarr, 'y':yarr, 'yerr':yerrarr}


def FlagLightCurve(pd):
    	""" This function flags element ids for known transits """

	# have the KOI database stuff in here
	# 
	# outputs are Period, t0, and dur

	for line in open("/astro/users/martincj/kepler_kois2.csv"):
		if '%s' % (pd['kid']) in line:
			element = line.split('|')
			period = float(element[7])
			t0 = float(element[5]) + 54900
			dur = float(element[2])
			
	transitslist = []
	transitelist = []
	
	d = (dur/24e0)/2e0
	
	tstime0 = t0 - d
	tetime0 = t0 + d
	
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
	
	#return {'kid':pd['kid'],'x':num.array(outtime),'y':num.array(outflux),'yerr':num.array(outerr)}

def FlagOutliers(d2,medwin,threshold):
    """ This function flags outliers ... """
    
    dout = {}
    for portion in d2.keys():
        
        npts = len(d2[portion]['x'])
        medflux = []
        medhalf = (medwin-1)/2
        OutlierFlag = []
        
        for i in range(npts):
                i1 = max(0,i-medhalf)
                i2 = min(npts, i + medhalf)
                medflux.append(num.median(d2[portion]['y'][i1:i2]))
        
        medflux = num.array(medflux)
        
        outliers = d2[portion]['y'] - medflux
        
        outliers.sort()
        sigma = (outliers[.8415*npts]-outliers[.1585*npts])/2
        outliers = d2[portion]['y'] - medflux
        
        check = abs(outliers)<threshold*sigma
        
        #outliers = outliers*(abs(outliers)<threshold*sigma)+medflux
        for el in check:
            if el == 1:
                OutlierFlag.append(False) 
            if el == 0:
                OutlierFlag.append(True)
                
        dout[portion] = {'kid':d2[portion]['kid'],'x':num.array(d2[portion]['x']),'y':num.array(outliers),'yerr':num.array(d2[portion]['yerr']),'TransitFlag':d2[portion]['TransitFlag'],'OutlierFlag':OutlierFlag}
        
    return dout
	
def detrend(d3, window, order):
    
    dout = {}
    for portion in d3.keys():
        
        
        npts = len(d3[portion]['x'])
        
        t0= d3[portion]['x']
        y0= d3[portion]['y']
        dtreflux = []
        for i in range(npts):
            tlow_start = t0[j] - win1
            tlow_end = t0[j] - win2
            thigh_start = t0[j] + win2
            thigh_end = t0[j] + win1
            t_test =
            if t_test > tlow_start and  t_test < tlow_end:
                indwin.append(t_test.ij)
                #need to add outlier flags so we can check
                if ((t0[j] > t0[i] - win2) and (t0[j] < t0[i] - win1)) or ((t0[j] > t0[i]+ win1) and (t0[j] < t0[i] + win2)):
                
                    indwin = []
        
            #print len(indwin)
            #print indwin
                #indwin = num.where(((t0 > t0[i] - win2) and (t0 < t0[i] - win1)) or ((t0 > t0[i]+ win1) and (t0 < t0[i] + win2)))
                #indwin = (((t0 > t0[i] - win2) and (t0 < t0[i] - win1)) or ((t0 > t0[i]+ win1) and (t0 < t0[i] + win2))).nonzero()
                #print indwin
            twin = t0[indwin] - t0[i]
            fwin = y0[indwin]
                
            coeff = scipy.polyfit(twin, fwin, order)
            dtreflux.append(coeff[0])
        dout[portion] = {'kid':d3[portion]['kid'],'x':num.array(d3[portion]['x']),'y':num.array(dtreflux),'yerr':num.array(d3[portion]['yerr']),'TransitFlag':num.array(d3[portion]['TransitFlag']), 'OutlierFlag':num.array(d3[portion]['OutlierFlag'])}
        
    return dout
                
def detrend2(d3, portion, order):
    
    
    
    d_noout = rejectOutliers(d3)
    d = rejectTransit(d_noout)
    
    coeff = scipy.polyfit(d[portion]['x'], d[portion]['y'], order)
    #dtreflux.append(coeff[0])
    out = scipy.polyval(coeff,d3[portion]['x'])
    import pylab
    pylab.plot(d[portion]['x'],d[portion]['y'],'b.')
    pylab.plot(d3[portion]['x'],out,'b.')
    pylab.show()
    return {'x':d3[portion]['x'],'y':out}

                
def rejectOutliers(d):
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
            
def rejectTransit(d):
    """ returns array without outliers """
    
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