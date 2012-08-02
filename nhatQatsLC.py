import uwpyKepler as kep
import numpy as num
from uwpyKepler import qats

class nhatQatsLC(kep.qats.qatslc):
    def __init__(self, lcData, kid, period):
        self.lcData = lcData
        self.kid = kid
        self.period = period
    
    def runNQATS(self, **kwargs):
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

        #getting resolved period from discoveries.txt or just going with default

        period0 = self.period
        # convert period0 to cadence
        period0 = period0*24.0*60.0/29.425
        
        ds = {'rho_s':rho_s,'b':b,'dt':self.dt}
        self.ds = ds
        period, snr0, nhat0, tmin, tmax = period_snr_nhat(self.lcData['y'],period0,f,self.nperiod,ds)
        period, snr1, nhat1, tmin, tmax = period_snr_nhat(self.lcData['flat'],period0,f,self.nperiod,ds)
        
        self.periods = period
        self.tmin = tmin
        self.tmax = tmax
        self.nhat = nhat0
        self.snrLC = snr0/self.sigma
        self.snrFLAT = snr1/self.sigma
        self.SignalPower = snr0/snr1
    
def period_snr_nhat(data,period0,f,nperiod,ds):
    period = []
    speriod = []
    mbest = []
    qbest = []
    


    spmax = 0e0
    #period0 = (period0*(1+f/2)/(1-f/2))
    period.append(period0)
    tmin = num.floor(period0*(1e0-f/2e0))
    tmax = num.ceil(period0*(1e0+f/2e0))
    q = num.floor(qats.tdur(ds['rho_s'],ds['b'],period0)/ds['dt'])
    MM, nhat, smax, dc = \
    qpt_detect(data,tmin,tmax,q)
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

def gamma_m(m,tmin,tmax,n,q):
    """
    """
    
    j1 = max([(m-2L)*tmin,n-tmax])
    j2 = min([tmax-q+(m-2L)*tmax,n-tmin])
    
    return j1,j2

def omega_m(m,mm,tmin,tmax,n,q):
    """
    """
    
    j1 = max([(m-1L)*tmin,n-(mm-m+1)*tmax])
    j2 = min([-q+m*tmax,n-q-(mm-m)*tmin])
    out = (j1+num.arange(j2-j1+1,dtype=num.int64)).tolist()
    #print j1, j2, out
    return out

def fmax(MM,N,tmin,tmax,q,dc):
    """
    """

    # Computes F_max
    shape = (MM,N)
    fmnMM = num.zeros(shape)
    for m in range(1,MM+1,1):
        # Compute omega:
        om = omega_m(m,MM,tmin,tmax,N,q)
        if (m == 1):
            fmnMM[m-1][om] = dc[om]
        else:
            for i in range(len(om)):
                j = gamma_m(m,tmin,tmax,om[i],q)
                #print num.shape(fmnMM[m-2][j[0]:j[1]+1]), j[0],j[1], om[i]
                fmnMM[m-1][om[i]] = dc[om[i]]\
                +max(fmnMM[m-2][j[0]:j[1]+1])
    
    singlefmnMN = fmnMM[MM-1][om].ravel()
    fmax0 = max(singlefmnMN)
    imax = singlefmnMN.argmax()
    nhat = num.zeros(MM)
    nhat[MM-1] = om[imax]
    for m in range(MM-1,0,-1):
        gam = gamma_m(m+1,tmin,tmax,nhat[m],q)
        tmp = max(fmnMM[m-1][gam[0]:gam[1]+1])
        imax = fmnMM[m-1][gam[0]:gam[1]+1].argmax()
        nhat[m-1] = gam[0]+imax
        
    #print num.shape(nhat), num.shape(fmax0)
    return fmax0,nhat
    
def qpt_convolve(data,q):
    """
    """
    
    tmp=num.convolve(data,num.ones(q),'valid')
    return tmp

def qpt_detect(data,tmin,tmax,q):
    """
    Uses the Kel'Manov algorithm for detection of
    quasi-periodic transits
    """
    
    # Subtract the data from the mean:
    data=num.mean(data)-data
    N = len(data)
    # Compute d(u,n):
    #dc=2d0*qpt_convolve(data,q)-double(q)
    dc=qpt_convolve(data,q)
    # Minimum and maximum number of transits given
    # the duration of the data, q & (tmin,tmax)
    mmin=long(num.floor((N+q-1L)/tmax))
    mmax=long(num.floor((N-q)/tmin)+1L)
    #print 'Range of transit numbers: ',mmin,mmax
    smaxMM = num.zeros(mmax-mmin+1)
    nhatMM = num.zeros( (mmax-mmin+1,mmax) )
    # Loop over the number of transits, MM:
    #print mmin, mmax
    #print num.shape(smaxMM)
    for MM in range(mmin,mmax+1,1):
        # for MM=mmin,mmin do begin
        # Optimize the likelihood for a given number of transits:
        # smaxMM[MM-mmin]=fmax(MM,N,tmin,tmax,q,dc,nhat)^2/double(q*MM)
        smaxMM[MM-mmin],nhat = fmax(MM,N,tmin,tmax,q,dc)
        nhatMM[MM-mmin][0:MM] = nhat

    smax = max(smaxMM)
    imax = smax.argmax()
    MM = mmin+imax
    nhat = (nhatMM[imax][0:MM]).ravel()
    
    return MM,nhat,smax,dc