import uwpyKepler as kep
import numpy as num
import nhat_qats as qats
from uwpyKepler.lcmod import returnData
from uwpyKepler import func
import sys
import pdb
import pylab as py
import scipy

def getPer():
    eDatafileDir = '/astro/store/student-scratch1/martincj/condor/discoveries/'
    name = 'discoveries.txt'
    dFile = open(eDatafileDir + name, 'r')
    lines = dFile.readlines()
    for line in lines:
        if line.split()[0] == str(self.kid):
            line   = line.split()
            period = float(line[1])
            t0     = float(line[2])
            q      = float(line[3])
    return period
###older qatspars
class qatslc:
    def __init__(self, lc):
        self.lcData = lc.lcFinal
        self.kid = lc.KID
    
    def addNoise(self,**kwargs):
        """
        add noise to padded data
        """
        
        # Default Sigma
        Sigma = self.sigma
        for key in kwargs:
            if key.lower() == 'noise':
                Sigma = kwargs[key]
            else:
                continue
            
        NoiseIDs = num.where(self.lcData['yerr'] == 0e0)[0]
        self.lcData['y'][NoiseIDs] += Sigma*num.random.randn(len(NoiseIDs))
        self.status = 'Noise Added'
    
    def padLC(self,**kwargs):
        """
        Pads missing cadences with ones.
        Returns the padded lightcurve.
        
        FlagIDs must corresepond to input lcData indices
        i.e. indices of the array from the pipline
        """
        
        #Default = no flaged points
        FlagIDs = []
        for key in kwargs:
            if key.lower() == 'flagids':
                FlagIDs = kwargs[key]

        # the eclipse lightcurve data are the best place to start
        
        if len(FlagIDs) > 0:
            x,y,yerr,cad = returnData(self.lcData,'all')
            if len(x)-1 < max(FlagIDs):
                raise NameError("len(x) < max(FlagIds) "+\
                str(len(x))+" < "+str(max(FlagIDs)) )
            x = x[FlagIDs]
            y = y[FlagIDs]
            yerr = yerr[FlagIDs]
            cad = cad[FlagIDs]
        else:
            x,y,yerr,cad = returnData(self.lcData,'elc')
            pass

        x0 = min(x)
        x1 = max(x)
        c0 = min(cad)
        c1 = max(cad)
        dx = num.median(x[1L:len(x)]\
                            -x[0L:len(x)-1])
        Tcad = num.median(dx)
        NCad = c1-c0+1
        xcomplete = x0 + Tcad*num.arange(NCad)
        cadcomplete = c0 + num.arange(NCad)
        

	test = num.arange(len(xcomplete))
        #print test[0], test[-1]

        #print len(xcomplete)
        #print len(cadcomplete)
        
        if len(xcomplete) != len(cadcomplete):
            print len(xcomplete), len(cadcomplete)
            raise NameError("Mismatch in Cadence and Time Intervals")
        
        #missing = num.array(list(set.difference(set(cadcomplete),set(cad))))
        #missingIDX = getListIndicies(cadcomplete,missing)
        existingIDX = getListIndicies(cadcomplete,cad)
    
        zeros = num.zeros(len(xcomplete))
        yerrcomplete = zeros                #padding errors with 0
        padflag = zeros
        ycomplete = zeros+1e0               #padding lc with 1
        
        # padded datasets
        # re-using original time-stamps for existing cadences
        xcomplete[existingIDX] = x
        ycomplete[existingIDX] = y
        yerrcomplete[existingIDX] = yerr
        sigma = func.compute1Sigma(ycomplete[existingIDX])
        flatgauss = 1e0 + sigma*num.random.randn(len(ycomplete))
        self.sigma = sigma
        qflag = zeros[existingIDX] = 1
        self.lcData = {'x':xcomplete,\
                       'y':ycomplete,\
                       'yerr':yerrcomplete,\
                       'flat':flatgauss,\
                       'padflag':padflag}
        self.status = 'Padded Lightcurve'
    
    def runQATS(self, **kwargs):
        NPoints = len(self.lcData['x'])
        #gets differences between x values
        self.dt = num.median(self.lcData['x'][1L:NPoints]\
                            -self.lcData['x'][0L:NPoints-1])
    
        # default pmin, pmax
        pmin = long(num.floor(1.3e0/self.dt))
        pmax = long(num.ceil(100e0/self.dt))
        #print pmin, pmax
        f = 0.005e0
        b = 0e0
        rho_s = qats.stellar_dens(self.kid)
        if rho_s == -99: rho_s = 1.4  #average solar density in g/cc
    
        self.pmin = pmin
        self.pmax = pmax
        self.f = f
        self.nperiod = long(num.log(pmax/pmin)/num.log((1+f/2)/(1-f/2)))
        period0 = pmin/(1+f/2)*(1-f/2)
        
        ds = {'rho_s':rho_s,'b':b,'dt':self.dt}
        self.ds = ds
        period, snr0, nhat0 = period_snr(self.lcData['y'],period0,f,self.nperiod,ds)
        period, snr1, nhat1 = period_snr(self.lcData['flat'],period0,f,self.nperiod,ds)
        
        self.periods = period
        self.snr0 = snr0
        self.snr1 = snr1
        self.nhat0 = nhat0
        self.nhat1 = nhat1

    def runQATS2(self, **kwargs):
        
        NPoints = len(self.lcData['x'])
        self.dt = num.median(self.lcData['x'][1L:NPoints]\
                            -self.lcData['x'][0L:NPoints-1])

        # default pmin, pmax
        pmin = long(num.floor(1.3e0/self.dt))
        pmax = long(num.ceil(100e0/self.dt))
        #print pmin, pmax
        f = 0.005e0
        b = 0e0
        rho_s = stellar_dens(self.kid)
        if rho_s == -99: rho_s = 1.4  #average solar density in g/cc

        self.pmin = pmin
        self.pmax = pmax
        self.f = f
        self.nperiod = long(num.log(pmax/pmin)/num.log((1+f/2)/(1-f/2)))
        period0 = pmin/(1+f/2)*(1-f/2)
        self.ds = {'rho_s':rho_s,'b':b,'dt':self.dt}
        pDict, periods, qlist = mkpDict(period0,self.nperiod,f,NPoints,self.ds)
        data0 = self.lcData['y']-num.mean(self.lcData['y'])
        data1 = self.lcData['flat']-num.mean(self.lcData['flat'])

        dcData0 = qats_struct.qpt_convolve_multi(data0,qlist)
        dcData1 = qats_struct.qpt_convolve_multi(data1,qlist)
        
        snr0 = qats_struct.qpt_detect_multi(dcData0,pDict,NPoints)
        snr1 = qats_struct.qpt_detect_multi(dcData1,pDict,NPoints)
        
        self.periods = periods
        self.snr0 = snr0
        self.snr1 = snr1

    def runNQATS(self, **kwargs):
        NPoints = len(self.lcData['x'])
        #gets differences between x values
        self.dt = num.median(self.lcData['x'][1L:NPoints]\
                            -self.lcData['x'][0L:NPoints-1])
    
        # default pmin, pmax
        pmin = long(num.floor(1.3e0/self.dt))
        pmax = long(num.ceil(100e0/self.dt))
        #print pmin, pmax
        f = 0.02e0
        b = 0e0
        rho_s = qats.stellar_dens(self.kid)
        if rho_s == -99: rho_s = 1.4  #average solar density in g/cc
    
        self.pmin = pmin
        self.pmax = pmax
        self.f = f
        self.nperiod = long(num.log(pmax/pmin)/num.log((1+f/2)/(1-f/2)))

        #getting resolved period from discoveries.txt or just going with default

        period0 = kep.analysis.geteDataFromFile(self.kid)[0]
        if period0 == -1:
            print 'object has no data entered in discoveries.txt'
            print 'using default QATS best period'
            period0 = kep.postqats.getBestPeriodByKID(self.kid)
        print 'using period of:', period0, ' days'
        # convert period0 to cadence
        period0 = period0*24.0*60.0/29.425
        
        ds = {'rho_s':rho_s,'b':b,'dt':self.dt}
        self.ds = ds
        period, snr0, nhat0, tmin, tmax = period_snr_nhat(self.lcData['y'],period0,f,self.nperiod,ds)
        period, snr1, nhat1, tmin, tmax = period_snr_nhat(self.lcData['flat'],period0,f,self.nperiod,ds)
        
        self.periods = period
        self.snr0 = snr0
        self.snr1 = snr1
        self.nhat0 = nhat0
        self.nhat1 = nhat1
        self.snrLC = snr0/self.sigma
        self.snrFLAT = snr1/self.sigma
        self.SignalPower = snr0/snr1
     
    def runNQATSHD(self, tmin, tmax, **kwargs):
        NPoints = len(self.lcData['x'])
        #gets differences between x values
        self.dt = num.median(self.lcData['x'][1L:NPoints]\
                            -self.lcData['x'][0L:NPoints-1])
    
        # default pmin, pmax
        pmin = long(num.floor(1.3e0/self.dt))
        pmax = long(num.ceil(100e0/self.dt))
        #print pmin, pmax
        f = 0.02e0
        b = 0e0
        rho_s = qats.stellar_dens(self.kid)
        if rho_s == -99: rho_s = 1.4  #average solar density in g/cc
    
        self.pmin = pmin
        self.pmax = pmax
        self.f = f
        self.nperiod = long(num.log(pmax/pmin)/num.log((1+f/2)/(1-f/2)))

        #getting resolved period from discoveries.txt or just going with default

        period0 = kep.analysis.geteDataFromFile(self.kid)[0]
        if period0 == -1:
            print 'object has no data entered in discoveries.txt'
            print 'using default QATS best period'
            period0 = kep.postqats.getBestPeriodByKID(self.kid)
        print 'using period of:', period0, ' days'
        # convert period0 to cadence
        period0 = period0*24.0*60.0/29.425
        
        ds = {'rho_s':rho_s,'b':b,'dt':self.dt}
        self.ds = ds
        period, snr0, nhat0, tmin1, tmax1 = period_snr_nhatHD(self.lcData['y'],period0,f,self.nperiod,ds, tmin, tmax)
        period, snr1, nhat1, tmin2, tmax2 = period_snr_nhatHD(self.lcData['flat'],period0,f,self.nperiod,ds, tmin, tmax)
        
        self.periods = period
        self.snr0 = snr0
        self.snr1 = snr1
        self.nhat0 = nhat0
        self.nhat1 = nhat1
        self.snrLC = snr0/self.sigma
        self.snrFLAT = snr1/self.sigma
        self.SignalPower = snr0/snr1
        
### period_snr
def period_snr(data,period0,f,nperiod,ds):
    
    period = []
    speriod = []
    mbest = []
    qbest = []
    
    for i in range(nperiod):
    #for i in range(1):
        print period0
        spmax = 0e0
        period0 = (period0*(1+f/2)/(1-f/2))
        period.append(period0)
        tmin = num.floor(period0*(1e0-f/2e0))
        tmax = num.ceil(period0*(1e0+f/2e0))
        q = num.floor(qats.tdur(ds['rho_s'],ds['b'],period0)/ds['dt'])
        MM, nhat, smax, dc = \
        qats.qpt_detect(data,tmin,tmax,q)
        speriod.append(smax)
        if smax > spmax:
            mbest.append(MM)
            qbest.append(q)
            spmax = smax

    period = num.array(period)
    speriod = num.array(speriod)
    mbest = num.array(mbest)
    qbest = num.array(qbest)
    snr = speriod/num.sqrt(mbest*qbest)

    return period, snr, nhat

def period_snr_nhat(data,period0,f,nperiod,ds):
    period = []
    speriod = []
    mbest = []
    qbest = []
    


    print period0
    spmax = 0e0
    #period0 = (period0*(1+f/2)/(1-f/2))
    period.append(period0)
    tmin = num.floor(period0*(1e0-f/2e0))
    tmax = num.ceil(period0*(1e0+f/2e0))
    q = num.floor(qats.tdur(ds['rho_s'],ds['b'],period0)/ds['dt'])
    MM, nhat, smax, dc = \
    qats.qpt_detect(data,tmin,tmax,q)
    speriod.append(smax)
    if smax > spmax:
        mbest.append(MM)
        qbest.append(q)
        spmax = smax

    period = num.array(period)
    speriod = num.array(speriod)
    mbest = num.array(mbest)
    qbest = num.array(qbest)
    snr = speriod/num.sqrt(mbest*qbest)

    return period, snr, nhat, tmin, tmax
    
def period_snr_nhatHD(data,period0,f,nperiod,ds,tmin,tmax):
    period = []
    speriod = []
    mbest = []
    qbest = []
    


    print period0
    spmax = 0e0
    #period0 = (period0*(1+f/2)/(1-f/2))
    period.append(period0)
    #tmin = num.floor(period0*(1e0-f/2e0))
    #tmax = num.ceil(period0*(1e0+f/2e0))
    q = num.floor(qats.tdur(ds['rho_s'],ds['b'],period0)/ds['dt'])
    MM, nhat, smax, dc = \
    qats.qpt_detect(data,tmin,tmax,q)
    speriod.append(smax)
    if smax > spmax:
        mbest.append(MM)
        qbest.append(q)
        spmax = smax

    period = num.array(period)
    speriod = num.array(speriod)
    mbest = num.array(mbest)
    qbest = num.array(qbest)
    snr = speriod/num.sqrt(mbest*qbest)

    return period, snr, nhat, tmin, tmax

### PART OF qatspars.py ###
def mkpDict(period0,nperiod,f,N,ds):
    
    periods = [period0]
    for ip in range(nperiod):
        periods.append(periods[-1]*(1+f/2)/(1-f/2))
    
    periods = num.array(periods)
    tdur_vec = num.vectorize(tdur)
    q = num.int_(num.floor(tdur_vec(ds['rho_s'],ds['b'],periods)/ds['dt']))
    tmin = num.floor(periods*(1e0-f/2e0))
    tmax = num.ceil(periods*(1e0+f/2e0))
    mmin= num.int_(num.floor((N+q-1L)/tmax))
    mmax= num.int_(num.floor((N-q)/tmin)+1L)
    
    pDdict = {}
    for ip in range(nperiod):
        pDdict[ip] = {'q':q[ip],
                      'tmin':tmin[ip],
                      'tmax':tmax[ip],
                      'mmin':mmin[ip],
                      'mmax':mmax[ip],
                      'period':periods[ip]
                     }
                     
    return pDdict, periods, q




def gamma_m_list(mlist,tmin,tmax,omDict,q):
    """
    """
    
    outDict = {}
    for m in mlist:
        tempDict = {}
        for MM in omDict[m].keys():
            tempDict[MM] = {}
            for n in range(omDict[m][MM]['k0'],omDict[m][MM]['k1']+1,1):
                tempDict[MM][n] =\
                {'j0':num.int_(max([(m-2L)*tmin,n-tmax])),
                 'j1':num.int_(min([tmax-q+(m-2L)*tmax,n-tmin]))}
        outDict[m] = tempDict
        
    return outDict

def omega_m_list(mlist,MMin,MMax,tmin,tmax,n,q):
    """
    """
    
    outDict = {}
    for m in mlist:
        tempDict = {}
        for MM in range(MMin,MMax+1,1):
            tempDict[MM] = {'k0':num.int_(max([(m-1L)*tmin,n-(MM-m+1)*tmax])),
                            'k1':num.int_(min([-q+m*tmax,n-q-(MM-m)*tmin]))}
                    
        outDict[m] = tempDict
        
    return outDict
    
def qpt_convolve(data,q):
    """
    """
    
    tmp=num.convolve(data,num.ones(q),'valid')
    return tmp
    
def qpt_convolve_multi(data,qlist):
    
    unique_q = list(set(qlist))
    dcDict = {}
    for uq in unique_q:
        dcDict[uq] = qpt_convolve(data,uq)
    
    return dcDict

def fmax_multi(MMin,MMax,Ndata,tmin,tmax,q,dc):
    """
    """
    
    mList = range(1,MMax+1,1)
    
    # Compute omega:
    omDict = omega_m_list(mList,MMin,MMax,tmin,tmax,Ndata,q)
    gmDict = gamma_m_list(mList,tmin,tmax,omDict,q)
    fmaxL = []
    
    for MM in range(MMin,MMax+1,1):
        shape = (MM,Ndata)
        fmnMM = num.zeros(shape)
        for m in mList:
            #om = qats.omega_m(m,MM,tmin,tmax,Ndata,q)
            if m <= MM:
                if m == 1:
                    fmnMM[m-1][omDict[m][MM]['k0']:omDict[m][MM]['k1']+1] =\
                    dc[omDict[m][MM]['k0']:omDict[m][MM]['k1']+1]
                else:
                    maxArr = []
                    #i = 0
                    for k in range(omDict[m][MM]['k0'],omDict[m][MM]['k1']+1,1):
                        #gm = qats.gamma_m(m,tmin,tmax,om[i],q)
                        #i += 1
                        test = fmnMM[m-2][gmDict[m][MM][k]['j0']:gmDict[m][MM][k]['j1']+1]
                        maxArr.append(max(test))
                        #print gm, gmDict[m][MM][k]['j0'],gmDict[m][MM][k]['j1'], len(test)
                    #print len(om), len(dc[omDict[m][MM]['k0']:omDict[m][MM]['k1']+1]), len(maxArr)
                    #print num.shape(dc), omDict[m][MM]['k0'],omDict[m][MM]['k1'], m, MM
                    fmnMM[m-1][omDict[m][MM]['k0']:omDict[m][MM]['k1']+1] =\
                    dc[omDict[m][MM]['k0']:omDict[m][MM]['k1']+1]+\
                    num.array(maxArr)
        
                fmax0 = max(fmnMM[MM-1][omDict[m][MM]['k0']:omDict[m][MM]['k1']+1].ravel())
                #fmax1 = max(fmnMM[MM-1][om].ravel())
                #print om, omDict[m][MM]['k0'],omDict[m][MM]['k1'], num.shape(fmnMM), m, MM
                #print fmax0, fmax1
                fmaxL.append(fmax0)
        
    return fmaxL

def qpt_detect_multi(dcData,pData,Ndata):
    
    spmax = 0
    snr = []
    for ip in pData.keys():
        print ip,pData[ip]['period'],pData[ip]['q']
        fMax =  fmax_multi(pData[ip]['mmin'],\
                           pData[ip]['mmax'],\
                           Ndata,\
                           pData[ip]['tmin'],\
                           pData[ip]['tmax'],\
                           pData[ip]['q'],\
                           dcData[pData[ip]['q']])
        smax = max(fMax)
        imax = smax.argmax()
        MM = pData[ip]['mmin']+imax

        if smax > spmax:
            snr.append(smax/(MM*pData[ip]['q']))
            spmax = smax
            
    return snr

def getListIndicies(Array,ListValues):
    """
    Return a list of indices where a list of values exists
    in a given array
    """ 
    
    lArray = Array.tolist()
    lLV = ListValues.tolist()
    indX = [lArray.index(x) if x in lLV else None for x in lLV]
    
    return indX

def main():
    kid = int(sys.argv[1])
    lc = kep.keplc.keplc(kid)
    kw = kep.quicklc.quickKW()
    lc.runPipeline(kw)
    qlc = qatslc(lc)
    qlc.padLC()
    qlc.addNoise()
    qlc.runNQATS()
    nhat0 = num.array(map(int, list(qlc.nhat0)))
    py.plot(qlc.lcData['x'], qlc.lcData['y'], '.b')
    py.plot(qlc.lcData['x'][nhat0], qlc.lcData['y'][nhat0], 'ro')
    py.show()
    #pdb.set_trace()
    coeff = num.polyfit(num.log10(qlc.periods),num.log10(qlc.snrLC),1)
    outy = scipy.polyval(coeff,num.log10(qlc.periods))
    normalizedPower = 10**(outy)
    
    #py.plot(qlc.periods,qlc.snrFLAT,'ro') 
    #py.plot(qlc.periods,qlc.snrLC,'bo')
    #py.plot(qlc.periods,normalizedPower,'ko')
    #py.setp(py.gca().set_xscale('log'))
    #py.show()
    #py.savefig('sn.'+str(kid)+'.png')
    
    dfile = open('signal.'+str(kid)+'.data','w')
    print >> dfile,'#',kid,'|', qlc.periods[num.argmax(qlc.snrLC/normalizedPower)],\
            qlc.periods[num.argmax(qlc.snrLC)], max(qlc.snrLC/normalizedPower), max(qlc.snrLC)
    for i in range(len(qlc.SignalPower)):
        print >> dfile, qlc.periods[i],'|',qlc.snrLC[i],'|',\
                        qlc.snrFLAT[i],'|',normalizedPower[i],'|',\
                        qlc.nhat0[i],'|',qlc.nhat1[i]
    dfile.close()
    
def mainHD(tmin,tmax, kid):
    lc = kep.keplc.keplc(kid)
    kw = kep.quicklc.quickKW()
    lc.runPipeline(kw)
    qlc = qatslc(lc)
    qlc.padLC()
    qlc.addNoise()
    qlc.runNQATSHD(tmin,tmax)
    nhat0 = num.array(map(int, list(qlc.nhat0)))
    py.clf()
    py.plot(qlc.lcData['x'], qlc.lcData['y'], '.b')
    py.plot(qlc.lcData['x'][nhat0], qlc.lcData['y'][nhat0], 'ro')
    py.savefig('nhat0_' + kid + '.png')
    #py.show()
    #pdb.set_trace()
    coeff = num.polyfit(num.log10(qlc.periods),num.log10(qlc.snrLC),1)
    outy = scipy.polyval(coeff,num.log10(qlc.periods))
    normalizedPower = 10**(outy)
    


    
    return qlc.nhat0
    #print qlc.snr0
#main()