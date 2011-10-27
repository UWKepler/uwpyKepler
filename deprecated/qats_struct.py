
import numpy as num
import qats

###older qatspars
def runQATS(self, **kwargs)
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
        
        ds = {'rho_s':rho_s,'b':b,'dt':self.dt}
        self.ds = ds
        period, snr0 = period_snr(self.lcData['y'],period0,f,self.nperiod,ds)
        period, snr1 = period_snr(self.lcData['flat'],period0,f,self.nperiod,ds)
        
        self.periods = period
        self.snr0 = snr0
        self.snr1 = snr1

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
        q = num.floor(tdur(ds['rho_s'],ds['b'],period0)/ds['dt'])
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

    return period, snr

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