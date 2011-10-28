import numpy as num
import sys

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
    
