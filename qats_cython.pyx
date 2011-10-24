
#from __future__ import division

import numpy as num

#compile time naming info
cimport numpy as num

#cython fixing datatype for arrays 

# faster min max functions
cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

def gamma_m(long int m,
            long int tmin,
            long int tmax,
            long int n,
            long int q):
    
    j0 = int_max((m-2)*tmin,n-tmax)
    j1 = int_min(tmax-q+(m-2)*tmax,n-tmin)
    
    return j0,j1

def omega_m(long int m,
            long int mm, 
            long int tmin,
            long int tmax,
            long int n,
            long int q):

    k0 = int_max((m-1)*tmin,n-(mm-m+1)*tmax)
    k1 = int_min(-q+m*tmax,n-q-(mm-m)*tmin)
    #out = (j1+num.arange(j2-j1+1,dtype=num.int64)).tolist()
    #print j1, j2, out
    return k0,k1

def fmax(long int MM,
         long int N,
         long int tmin,
         long int tmax,
         long int q,
         num.ndarray dc not None):

    fmnMM = num.zeros( (MM,N) ) 
    
    for m in range(1,MM+1,1):
        # Compute omega:
        k0,k1 = omega_m(m,MM,tmin,tmax,N,q)
        if (m == 1):
            fmnMM[m-1][k0:k1+1] = dc[k0:k1+1]
        else:
            for k in range(k0,k1+1,1):
                j0,j1 = gamma_m(m,tmin,tmax,k,q)
                fmnMM[m-1][k] = dc[k]\
                +max(fmnMM[m-2][j0:j1+1])
    
    singlefmnMN = fmnMM[MM-1][k0:k1+1].ravel()
    fmax0 = max(singlefmnMN)
    imax = singlefmnMN.argmax()

    return fmax0
    
def qpt_convolve(num.ndarray data not None,
                 long int q):

    tmp=num.convolve(data,num.ones(q),'valid')
    return tmp

def qpt_detect(num.ndarray data not None,
               long int tmin,
               long int tmax,
               long int q):
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
    # Loop over the number of transits, MM:
    for MM in range(mmin,mmax+1,1):
        # for MM=mmin,mmin do begin
        # Optimize the likelihood for a given number of transits:
        # smaxMM[MM-mmin]=fmax(MM,N,tmin,tmax,q,dc,nhat)^2/double(q*MM)
        smaxMM[MM-mmin] = fmax(MM,N,tmin,tmax,q,dc)
    
    smax = max(smaxMM)
    imax = smax.argmax()
    MM = mmin+imax
    snr = smax/num.sqrt(MM*q)
    return snr
    
def snr(num.ndarray data not None,
        num.ndarray tmin_arr not None,
        num.ndarray tmax_arr not None,
        num.ndarray q_arr not None,
        long int nperiod):
    
    snr = num.zeros( nperiod )
    
    for ip in range(nperiod):
        spmax = 0
        tmin = tmin_arr[ip]
        tmax = tmax_arr[ip]
        q = q_arr[ip]
        snr[ip] = qpt_detect(data,tmin,tmax,q)

    return snr
