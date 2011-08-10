import MySQLdb
import sys
import numpy as num
import scipy
import pylab
import uwpyKepler as kep

def inSource(KeplerID):
    """ Checks if a certain KID exists in the source database. """
    
    db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
    cursor = db.cursor()
    foo    = 'select %s from source where (KEPLERID = %s)' % (KeplerID,KeplerID)
    cursor.execute(foo)
    results = cursor.fetchall()
    Exist = False
    if len(results) > 0:
        Exist = True
    
    return Exist
        
def inKEPPC(KeplerID):
    """ Checks if a certain KID exists in the KEPPC database. """
    
    db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
    cursor = db.cursor()
    foo1 = 'select %s from KEPPC where (KID = %s)' % (KeplerID,KeplerID)
    cursor.execute(foo1)
    results = cursor.fetchall()
    Exist = False
    if len(results) > 0:
        Exist = True
    
    return Exist
    
def inKEPFP(KeplerID):
    """ Checks if a certain KID exists in the KEPFP database. """
    
    db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
    cursor = db.cursor()
    foo1 = 'select %s from KEPFP where (KID = %s)' % (KeplerID,KeplerID)
    cursor.execute(foo1)
    results = cursor.fetchall()
    Exist = False
    if len(results) > 0:
        Exist = True
    
    return Exist

def getKIDsFP():
	"""
	 Gets a list of unique KID from False positive table
	"""
	db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
	cursor = db.cursor()
	foo1 = 'select KID from KEPFP'
	cursor.execute(foo1)
	results = cursor.fetchall()
	r1 = ([el[0] for el in results])
	results=set(r1)
    
	return results

def getKIDsPC():
	"""
	Gets a list of unique KID from Planet Candidates table
	"""
	db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
	cursor = db.cursor()
	foo1 = 'select KID from KEPPC'
	cursor.execute(foo1)
	results = cursor.fetchall()
	r1 = ([el[0] for el in results])
	results=set(r1)

	return results

def getKIDsSource():
	"""
	Gets a list of unique KID from source table
	"""
	
	db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
    	cursor = db.cursor()
    	foo1 = 'select distinct KEPLERID from source'
    	cursor.execute(foo1)
    	results = cursor.fetchall()
    	r1 = ["%s" % el[0] for el in results]
	
	return r1

def getKOI(KeplerID):
    """ returns all KOI IDs (Kepler Object of Interest IDs) """
    
    if inKEPPC(KeplerID):
        db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
        cursor = db.cursor()
        foo1 = 'select KOI from KEPPC where (KID = %s)' % (KeplerID)
        cursor.execute(foo1)
        results = cursor.fetchall()
    elif inKEPFP(KeplerID):
        db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
        cursor = db.cursor()
        foo1 = 'select KOI from KEPFP where (KID = %s)' % (KeplerID)
        cursor.execute(foo1)
        results = cursor.fetchall()
    else:
        results = None
    
    return results

def ReadLightCurve(KeplerID):
    """ This function reads the Kepler database and returns the 
    corrected lightcurve in a dictionary
    
    Input = KeplerID
    Output = Dictionary with time, flux and fluxerror data
    """
    
    if inSource(KeplerID):
        db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
        cursor = db.cursor()
        foo    = 'select * from source where (KEPLERID = %s)' % (KeplerID)
        cursor.execute(foo)
        results = cursor.fetchall()
        
        # reading time, corrected flux and flux errors
        time    = num.ma.array([x[2] for x in results])
        corflux = num.ma.array([x[7] for x in results])
        corerr  = num.ma.array([x[8] for x in results])
        
        idx = num.where((corflux>0)&(corerr>0))
        
        time    = time[idx]
        corflux = corflux[idx]
        corerr  = corerr[idx]
        
        return {'kid':KeplerID,'x':time,'y':corflux,'yerr':corerr}
    else:
        print 'Kepler ID %s not found in Kepler.source' % (KeplerID)
        return

def getEclipseData(data):
    """queries the database for Transit or Eclipse data """
    
    db     = MySQLdb.connect(host='tddb.astro.washington.edu', user='tddb', passwd='tddb', db='Kepler')
    cursor = db.cursor()
    r1 = ()
    d1 = {}
    if inKEPPC(data['kid']):
        foo1 = 'select KOI, Period, Dur, Epoch from KEPPC where (KID = %s)' % (data['kid'])
        cursor.execute(foo1)
        r1 = cursor.fetchall()
	for i in range(len(r1)):
		d1[format(r1[i][0],'.2f')] = {'period':r1[i][1],'duration':r1[i][2],'t0':r1[i][3]}
	
    elif inKEPFP(data['kid']):
        foo1 = 'select KOI, Period, Duration, Epoch from KEPFP where (KID = %s)' % (data['kid'])
        cursor.execute(foo1)
        r1 = cursor.fetchall()
	for i in range(len(r1)):
		d1[format(r1[i][0],'.2f')] = {'period':r1[i][1],'duration':r1[i][2],'t0':r1[i][3]}
	
    else:
        print 'Kepler ID %s not found in Kepler.KEPPC or Kepler.KEPFP' % (KeplerID)
        
    return d1

def FlagTransits(data,eclipseData):
    """ This function flags points within a tranit and
    applies a mask.
        
        Input = data dictionary
              = eclipseData, a list of period, dur, and t0 for eclipses and transits
        Output = data dictionary with addition of the keys
        'TransitMask' and 'UnMasked'
    """
    
    mask0=num.ma.getmaskarray(data['x'])
    data['UnMasked']=mask0
    i = 0
    for koi in eclipseData.keys():
	period = eclipseData[koi]['period']
	t0 = eclipseData[koi]['t0']
        dur = eclipseData[koi]['duration']
        dur = (1.2*dur/24e0)
        t0 = t0 + 54900e0
        width = dur/period
        maxphase=1-width/2
        minphase=width/2
        phase= (data['x']-t0)/period-(data['x']-t0)//period
        idx=num.where((phase>maxphase)|(phase<minphase))
        data['x'][idx]= num.ma.masked
        mask1=num.ma.copy(data['x'].mask)
        if i == 0:
            data['TransitMask']=mask1
        else:
            data['TransitMask']=num.ma.mask_or(mask1,data['TransitMask'])
	i +=1

    return data
        

def SplitGap(data,gapsize,medwin,fluxdiff):
    """
    This function finds gaps and splits data into portions.
    
    Input =  data - data dictionary
                gapsize - size of gap in days
            
    Output = new data dictionary with data split into portions
    """
    
    # defining new empty lists and stuff
    pcount=0
    istamps=[]
    outData={}
    
    data['x'].mask = data['UnMasked']
    data['y'].mask = data['UnMasked']
    data['yerr'].mask = data['UnMasked']
    
    # median smoothing the lightcurve
    mvavg1 = movingMedian(data['y'],medwin)
    mvavg1 = num.append(mvavg1,mvavg1[-1])
    mvavg1 = data['y']
    # first derivative of smoothed lightcurve
    diff1 = num.diff(mvavg1)
    diff1 = num.hstack((diff1,diff1[-1]))
    
    # second derivative of smoothed lightcurve
    diff2 = num.diff(diff1)
    diff2 = num.hstack((diff2[-1],diff2))

    # compute ourlier resistant sigma
    sig = compute1Sigma(diff1)
    
    pylab.plot(diff1,'g.')
    pylab.plot([0,6000],[5*sig,5*sig],'k-')
    pylab.plot([0,6000],[3*sig,3*sig],'k-')
    pylab.plot([0,6000],[1*sig,1*sig],'k-')
    pylab.show()

    # The grand master loop >=}
    # to make portion slices
    for i in range(len(data['x'])-1):
        dt =  data['x'][i+1]- data['x'][i]
        j1 = max(0,i-medwin)
        j2 = i + medwin
        if pcount == 0:
            i0 = 0
        if pcount > 0:
            i0 = i1+1
        if dt > gapsize:
            i1 = i
            istamps.append([i0,i1])
            pcount += 1
        if num.abs(diff1[i]) > 5*sig:
            i1 = i
            istamps.append([i0,i1])
            pcount += 1
            #print num.abs(diff1[i]/data['y'][i]), diff1[i], data['y'][i], diff1[i+1], data['y'][i+1]
            #print i, ' test flux gap'
    i1 = i+1
    istamps.append([i0,i1])
        
    # Applying slices
    for j in range(len(istamps)):
        print istamps[j][0], istamps[j][1]
        outData['portion' + str(j+1)] = {'kid':data['kid'],'x':data['x'][istamps[j][0]:istamps[j][1]+1], 'y':data['y'][istamps[j][0]:istamps[j][1]+1], 'yerr':data['yerr'][istamps[j][0]:istamps[j][1]+1], 'TransitMask':data['TransitMask'][istamps[j][0]:istamps[j][1]+1],'UnMasked':data['UnMasked'][istamps[j][0]:istamps[j][1]+1]}
        
    return outData
    
def FlagOutliers(data,medwin,threshold):
    """ This function flags outliers. 
        Inputs - data = data dictionary
               - medwin = the window size used to compute the median
               - threshold = the sigma-clipping factor (suggested, 3 or greater)
        Outputs - the data dictionary now contains mask arrays named
                'OutlierMask' and 'OTMask'
    """
    
    dout = {}
    # cycling through portions
    for portion in data.keys():
        data[portion]['x'].mask = data[portion]['TransitMask']
        data[portion]['y'].mask = data[portion]['TransitMask']
        data[portion]['yerr'].mask = data[portion]['TransitMask']
        npts = len(data[portion]['x'])
        
        # defining the window
        medflux = []
        medhalf = (medwin-1)/2

        # placing the window and computing the median
        for i in range(npts):
            i1 = max(0,i-medhalf)
            i2 = min(npts, i + medhalf)
            try:
                if num.ma.median(data[portion]['y'][i1:i2]).mask:
                    medflux.append(medflux[-1])
            except:
                medflux.append(num.ma.median(data[portion]['y'][i1:i2]))
            
            #if len(data[portion]['y'][i1:i2][num.where(data[portion]['y'][i1:i2].mask == True)]) > 0:
                #print i1, i2, npts, data[portion]['y'][i1:i2], medflux[-1]
                #print type(medflux[-1])
            
        # finding outliers
        medflux = num.array(medflux)
        outliers = num.ma.getdata(data[portion]['y']) - medflux
        
        sigma = compute1Sigma(outliers)
        
        outliers = data[portion]['y'] - medflux
        
        # tagging outliers (which are not part of the transit)
        idx=num.where( (abs(num.array(outliers))>threshold*sigma) & (data[portion]['TransitMask'] == False) )

        # creating the outlier mask
        data[portion]['x'].mask = data[portion]['UnMasked']
        data[portion]['x'][idx[0]] = num.ma.masked

        mask2 = num.ma.copy(data[portion]['x'].mask)
        
        data[portion]['OutlierMask']=mask2
        
        # creating the outlier + transit mask
        mask3 = num.ma.mask_or(data[portion]['TransitMask'],mask2)
        
        dout[portion] = {'kid':data[portion]['kid'],'x':data[portion]['x'],'y':data[portion]['y'],'yerr':data[portion]['yerr'],'TransitMask':data[portion]['TransitMask'],'UnMasked':data[portion]['UnMasked'],'OutlierMask':data[portion]['OutlierMask'],'OTMask':mask3}
        
    return dout
                    
def movingMedian(data,window):
    """
    Returns the median smoothed version of a given function.
    """
    
    mvavg = []
    npoints = len(data)
    for i in range(npoints-1):
        j1  = max(0,i-window)
        j2 = min(i+window,npoints-1)
        mvavg.append(num.median(data[j1:j2]))
        
    return num.array(mvavg)
    

def compute1Sigma(data):
    """
    
    """
    
    dsort = num.sort(data)
    npts = len(dsort)
    sigma = (dsort[.8415*npts]-dsort[.1585*npts])/2
    return sigma
                   
def ApplyMask(data,mask):
    """ This function applies a given mask """
    
    # loop through portions
    for portion in data.keys():
        # match data keys and apply mask    
        for key in data[portion].keys():
            if key in 'xyerr':
                if mask != 'UnMasked':
                    data[portion][key].mask = data[portion]['UnMasked']
                data[portion][key].mask = data[portion][mask]
	
    return data