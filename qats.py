import numpy as num
from lcmod import returnData
from dbinfo import returnLOGG, returnRstar
import qats_cython
import func

def getListIndicies(Array,ListValues):
    """
    Return a list of indices where a list of values exists
    in a given array
    """ 
    
    lArray = Array.tolist()
    lLV = ListValues.tolist()
    
    indX = [lArray.index(x) if x in lLV else None for x in lLV]
    
    return indX

def stellar_dens(KID):
    """
    compute stellar density
    """
    
    Logg = returnLOGG(KID)
    Rstar = returnRstar(KID)
    
    if Logg == -99 or Rstar == -99:
        rho_s = -99
    else:
        Grav = 6.674e-8
        Msol = 1.99e33
        Rsol = 6.96e10
        rho_s = 4e0*num.pi*(10**(Logg))/(Grav*Rstar*Rsol)
        
    return rho_s

def tdur(rho_s,b,period):
    """
    compute transit duration given stellar density,
    impact parameter and period
    """
    
    p_sec = period*86400e0
    Grav = 6.674e-8
    tdur = ((3e0*num.pi/Grav)**(1e0/3e0))*\
           ((1e0/num.pi)*(p_sec/rho_s)**(1e0/3e0))*\
           (num.sqrt(1e0-b**2))
         
    tdur = tdur/86400e0
    
    return tdur

def getPTQM(period0,nperiod,f,N,ds):
    
    periods = [period0]
    for ip in range(nperiod):
        periods.append(periods[-1]*(1+f/2))
    periods = num.array(periods)
    
    tdur_vec = num.vectorize(tdur)
    q = num.int_(num.floor(tdur_vec(ds['rho_s'],ds['b'],periods)/ds['dt']))
    tmin = num.int_(num.floor(periods*(1e0-f/2e0)))
    tmax = num.int_(num.ceil(periods*(1e0+f/2e0)))
    mmin= num.int_(num.floor((N+q-1L)/tmax))
    mmax= num.int_(num.floor((N-q)/tmin)+1L)
    
    return periods, tmin, tmax, mmin, mmax, q

class qatslc:

    def __init__(self,lcData,KID):
        
        self.lcData = lcData
        self.kid = KID
        self.status = 'Initial Input'

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
        xcomplete = num.arange(x0,x1,Tcad)
        cadcomplete = num.arange(c0,c1+1,1)

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
        
    def runQATS(self, **kwargs):
        
        NPoints = len(self.lcData['x'])
        self.dt = num.median(self.lcData['x'][1L:NPoints]\
                            -self.lcData['x'][0L:NPoints-1])

        # default pmin, pmax
        pmin = long(num.floor(1.3e0/self.dt))
        pmax = long(num.ceil(100e0/self.dt))
        # print pmin, pmax
        f = 0.005e0
        b = 0e0
        rho_s = stellar_dens(self.kid)
        if rho_s == -99: rho_s = 1.4  #average solar density in g/cc
        
        #sort through keywords
        for key in kwargs:
            if key.lower() == 'rho_s':
                rho_s = kwargs[key]
            elif key.lower() == 'pmin':
                pmin = long(num.floor(kwargs[key]/self.dt))
            elif key.lower() == 'pmax':
                pmax = long(num.floor(kwargs[key]/self.dt))
            elif key.lower() == 'f':
                f = kwargs[key]
            elif key.lower() == 'b':
                b = kwargs[key]
            else:
                pass

        nperiod = long(num.log(pmax/pmin)/num.log((1+f/2)))

        period0 = pmin/(1+f/2)
        ds = {'rho_s':rho_s,'b':b,'dt':self.dt}
        periods, tmin, tmax, mmin, mmax, q = getPTQM(period0,nperiod,f,NPoints,ds)
        
        snr0 = qats_cython.snr(self.lcData['y'],NPoints,tmin,tmax,q,nperiod+1)
        snr1 = qats_cython.snr(self.lcData['flat'],NPoints,tmin,tmax,q,nperiod+1)
        self.nperiod = nperiod+1
        self.periods = periods*self.dt
        self.Ndata = NPoints 
        self.snrLC = snr0/self.sigma
        self.snrFLAT = snr1/self.sigma
        self.SignalPower = snr0/snr1
