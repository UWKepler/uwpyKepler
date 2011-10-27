
from __future__ import division
import numpy as num
cimport numpy as num

#cython fixing datatype for arrays 
DTYPE = num.int
FTYPE = num.float64

#cdeftype
ctypedef num.int_t DTYPE_t
ctypedef num.float64_t FTYPE_t

# faster min max functions
cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

cdef extern from "math.h":
    float sqrt(float x)
    long floor(float x)
    long ceil(float x)

cdef long gamma_m0(long int m, long int tmin, long int tmax, long int n, long int q) except *:
    
    cdef int j0 = int_max((m-2)*tmin,n-tmax)
    return j0

cdef long gamma_m1(long int m, long int tmin, long int tmax, long int n, long int q) except *:
    
    cdef int j1 = int_min(tmax-q+(m-2)*tmax,n-tmin)
    return j1+1
    
cdef long omega_m0(long int m, long int mm, long int tmin, long int tmax, long int n, long int q) except *:
    
    cdef int k0 = int_max((m-1)*tmin,n-(mm-m+1)*tmax)
    return k0

cdef long omega_m1(long int m, long int mm, long int tmin, long int tmax, long int n, long int q) except *:
    
    cdef int k1 = int_min(-q+m*tmax,n-q-(mm-m)*tmin)
    return k1+1

def fmax(long int MM, long int N, long int tmin, long int tmax, long int q, num.ndarray[FTYPE_t, ndim=1] dc):

    cdef num.ndarray[FTYPE_t, ndim=2] fmnMM = num.zeros( [MM,N], dtype=FTYPE )
    cdef long int j0, j1
    cdef float fmax0
    cdef long int imax
    cdef num.ndarray[FTYPE_t, ndim=1] singlefmnMN
    cdef int m
    cdef int k
    
    for m from 1 <= m < MM+1:
        # Compute omega:
        k0 = omega_m0(m,MM,tmin,tmax,N,q)
        k1 = omega_m1(m,MM,tmin,tmax,N,q)
        if (m == 1):
            for k from k0 <= k < k1:
                fmnMM[m-1][k] = dc[k]
        else:
            for k from k0 <= k < k1:
                j0 = gamma_m0(m,tmin,tmax,k,q)
                j1 = gamma_m1(m,tmin,tmax,k,q)
                fmnMM[m-1][k] = dc[k]+max(fmnMM[m-2][j0:j1])
    
    singlefmnMN = fmnMM[MM-1][k0:k1].ravel()
    fmax0 = max(singlefmnMN)
    imax = singlefmnMN.argmax()

    return fmax0
    
def qpt_convolve(num.ndarray[FTYPE_t, ndim=1] data, long int q):

    tmp=num.convolve(data,num.ones(q),'valid')
    return tmp

def qatsconvolve(num.ndarray[FTYPE_t, ndim=1] data,long int Ndata, long int q):
    
    cdef num.ndarray[FTYPE_t, ndim=1] ker = num.ones(q)
    cdef num.ndarray[FTYPE_t, ndim=1] Rt = num.zeros(Ndata-q+1)
    cdef int k
    cdef int t
    cdef int i
    cdef int qi
    
    for qi from 0 <= qi < q:
        ker[qi] = 1

    for t  from 0 <= t < Ndata:
        for k from 0 <= k < q:
            if t >= k:
                Sum = 0
                for i from 0 <= i < k:
                    Sum = Sum + data[t-i]*ker[i]
            else:
                Sum = 0
        if t >= q:
            Rt[t-q+1] = Sum
    
    return Rt

def qpt_detect(num.ndarray[FTYPE_t, ndim=1] data, long int N, long int tmin, long int tmax, long int q):
    """
    Uses the Kel'Manov algorithm for detection of
    quasi-periodic transits
    """
    
    # Subtract the data from the mean:
    data = num.mean(data)-data
    cdef long int mmin = floor((N+q-1)/tmax)
    cdef long int mmax = floor((N-q)/tmin)+1
    cdef num.ndarray[FTYPE_t, ndim=1] smaxMM = num.zeros([mmax-mmin+1], dtype=FTYPE)
    cdef num.ndarray[FTYPE_t, ndim=1] dc = qpt_convolve(data,q)
    cdef int MM
    
    # Loop over the number of transits, MM:
    for MM from mmin <= MM < mmax+1:
        # for MM=mmin,mmin do begin
        # Optimize the likelihood for a given number of transits:
        # smaxMM[MM-mmin]=fmax(MM,N,tmin,tmax,q,dc,nhat)^2/double(q*MM)
        smaxMM[MM-mmin] = fmax(MM,N,tmin,tmax,q,dc)
    
    cdef float smax = max(smaxMM)
    cdef long int imax = smaxMM.argmax()
    MM = mmin+imax
    cdef float snr = smax/sqrt(MM*q)
    
    return snr
    
def snr(num.ndarray[FTYPE_t, ndim=1] data,long int Ndata, num.ndarray[DTYPE_t, ndim=1] tmin_arr, num.ndarray[DTYPE_t, ndim=1] tmax_arr, num.ndarray[DTYPE_t, ndim=1] q_arr, long int nperiod):
    
    #snr = num.zeros( nperiod )
    cdef num.ndarray[FTYPE_t, ndim=1] snr = num.zeros([nperiod], dtype=FTYPE)
    cdef long int tmin
    cdef long int tmax
    cdef long int q
    cdef int ip
    
    for ip from 0 <= ip < nperiod:
        tmin = tmin_arr[ip]
        tmax = tmax_arr[ip]
        q = q_arr[ip]
        snr[ip] = qpt_detect(data,Ndata,tmin,tmax,q)
        #print ip, snr[ip]

    return snr
