import numpy as num

def gamma_m(m,tmin,tmax,n,q):
    """
    """
    
    j1 = num.max([(m-2L)*tmin,n-tmax])
    j2 = num.min([tmax-q+(m-2L)*tmax,n-tmin])
    return j1, j2

def omega_m(m,mm,tmin,tmax,n,q):
    """
    """
    
    j1 = num.max([(m-1L)*tmin,n-(mm-m+1)*tmax])
    j2 = num.min([-q+m*tmax,n-q-(mm-m)*tmin])
    return,j1+lindgen(j2-j1+1)

def fmax(MM,N,tmin,tmax,q,dc,nhat):
    """
    Computes F_max
    """
    
    fmnMM=dblarr(MM,N)
    for m in range(MM)+1:
        # Compute omega
        om = omega_m(m,MM,tmin,tmax,N,q)
        if (m == 1):
            fmnMM[m-1,om]=dc[om]
        else:
            for i in range(len(om)):
                j = gamma_m(m,tmin,tmax,om[i],q)
                fmnMM[m-1,om[i]] = dc[om[i]]+num.max(fmnMM[m-2,j[0]:j[1]])

    fmax0=num.max(fmnMM[MM-1,om])
    imax = index(num.max(fmnMM[MM-1,om]))
    # Now determine optimum values:
    nhat=dblarr(MM)
    nhat[MM-1]=om[imax]
    for m in=MM-1,1,-1 do begin
        gam = gamma_m(m+1,tmin,tmax,nhat[m],q)
        tmp = num.max(fmnMM[m-1,gam[0]:gam[1]])
        imax = index(num.max(fmnMM[m-1,gam[0]:gam[1]]))
        nhat[m-1]=gam[0]+imax
    
    return fmax0

def qpt_convolve(data,q):
    """
    """
    
    tmp=convol(data,replicate(1d0,q),center=0)
    return tmp[(q-1L):*]

def qpt_detect(data,MM,tmin,tmax,q,dc,nhat,smax):
    """ Uses the Kel'Manov algorithm for detection of
    quasi-periodic transits
    """
        
    # Subtract the data from the mean:
    data=num.mean(data)-data
    N = n_elements(data)
    # Compute d(u,n):
    #dc=2d0*qpt_convolve(data,q)-double(q)
    dc=qpt_convolve(data,q)
    # Minimum and maximum number of transits given
    # the duration of the data, q & (tmin,tmax)
    mmin=floor((N+q-1L)/tmax) & mmax=floor((N-q)/tmin)+1L
    #print,'Range of transit numbers: ',mmin,mmax
    smaxMM=dblarr(mmax-mmin+1)
    nhatMM=dblarr(mmax-mmin+1,mmax)
    # Loop over the number of transits, MM:
    for MM=mmin,mmax do begin
        # Optimize the likelihood for a given number of transits:
        #  smaxMM[MM-mmin]=fmax(MM,N,tmin,tmax,q,dc,nhat)^2/double(q*MM)
        smaxMM[MM-mmin]=fmax(MM,N,tmin,tmax,q,dc,nhat)
        nhatMM[MM-mmin,0:MM-1]=nhat
        #  print,MM,smaxMM[MM-mmin]
    smax=max(smaxMM,imax)
    MM=mmin+imax
    nhat=reform(nhatMM[imax,0:MM-1])
    return


