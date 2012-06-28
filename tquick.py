""" Contains the functions that transitquick depends on. Some routines here are the translations of routines in E. Agol's transitquick.pro into python, written by J. Eastman. The remaining functions here are used by MultiTransitQuick. """

import numpy as np
import os
import sys
np.seterr(divide='ignore')

def ekepler(m,e):
    """ 
        This routine solves Kepler's equation for E as a function of (e,M)
        using the procedure outlined in Murray & Dermott:
        derived from E. Agol's implementation
    """

    eps = 1.e-10
    pi2 = 2.e0*np.arccos(-1e0)
    ms = (m % pi2)
    d3 = 1e10
    e0 = ms+e*0.85e0*np.sin(ms)/np.abs(np.sin(ms))
     
    while max(abs(d3)) > eps:
        f3=e*np.cos(e0)
        f2=e*np.sin(e0)
        f1=1.-f3
        f0=e0-ms-f2
        d1=-f0/f1
        d2=-f0/(f1+0.5*d1*f2)
        d3=-f0/(f1+d2*0.5*(f2+d2*f3/3.))
        e0=e0+d3

    ekep=e0+m-ms
    return ekep

def kepler(m,e):
    """ Translation of kepler.pro """
    
    nm = len(m)
    ekep = np.zeros(nm)
    i = long(0)
    
    if e != 0.e0:
        ekep = ekepler(m,e)
        f = 2e0*np.arctan(np.sqrt((1.e0+e)/(1.e0-e))*np.tan(0.5e0*ekep))
    else:
        f = m
    
    nm0 = []
    for idx in range(len(m)):
        if m[idx] == 0.e0:
            nm0.append(idx)

    if len(nm0) > 0:
        f[nm0] = 0.e0

    return f

def ellke(k):
    """ Computes Hasting's polynomial approximation for the complete
        elliptic integral of the first (ek) and second (kk) kind
    """

    m1=1.-k**2
    logm1 = np.log(m1)

    a1=0.44325141463
    a2=0.06260601220
    a3=0.04757383546
    a4=0.01736506451
    b1=0.24998368310
    b2=0.09200180037
    b3=0.04069697526
    b4=0.00526449639
    ee1=1.+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
    ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*(-logm1)
    ek = ee1+ee2
        
    a0=1.38629436112
    a1=0.09666344259
    a2=0.03590092383
    a3=0.03742563713
    a4=0.01451196212
    b0=0.5          
    b1=0.12498593597
    b2=0.06880248576
    b3=0.03328355346
    b4=0.00441787012
    ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
    ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*logm1
    kk = ek1-ek2
    
    return [ek,kk]

def ellpic_bulirsch(n,k):
    """Computes the complete elliptical integral of the third kind using
    the algorithm of Bulirsch (1965):
    """

    kc=np.sqrt(1.-k**2); p=n+1.
    if(np.min(p) < 0.):
        print 'Negative p'
    m0=1.; c=1.; p=np.sqrt(p); d=1./p; e=kc
    while 1:
        f = c; c = d/p+c; g = e/p; d = 2.*(f*g+d)
        p = g + p; g = m0; m0 = kc + m0
        if np.max(np.abs(1.-kc/g)) > 1.e-8:
            kc = 2*np.sqrt(e); e=kc*m0
        else:
            return 0.5*np.pi*(c*m0+d)/(m0*(m0+p))

def occultquad(z,u1,u2,p0):
    """ Python translation of IDL code.
    This routine computes the lightcurve for occultation of a
    quadratically limb-darkened source without microlensing.  Please
    cite Mandel & Agol (2002) and Eastman & Agol (2008) if you make use
    of this routine in your research.  Please report errors or bugs to
    jdeast@astronomy.ohio-state.edu
    """

    nz = np.size(z)
    lambdad = np.zeros(nz)
    etad = np.zeros(nz)
    lambdae = np.zeros(nz)
    omega=1.-u1/3.-u2/6.

    ## tolerance for double precision equalities
    ## special case integrations
    tol = 1e-14

    p = np.abs(p0)
    
    z = np.where(np.abs(p-z) < tol,p,z)
    z = np.where(np.abs((p-1)-z) < tol,p-1.,z)
    z = np.where(np.abs((1-p)-z) < tol,1.-p,z)
    z = np.where(z < tol,0.,z)
               
    x1=(p-z)**2.
    x2=(p+z)**2.
    x3=p**2.-z**2.
    
    ## trivial case of no planet
    if p <= 0.:
        muo1 = np.zeros(nz) + 1. 
        mu0  = np.zeros(nz) + 1.
        return [muo1,mu0]

    ## Case 1 - the star is unocculted:
    ## only consider points with z lt 1+p
    notusedyet = np.where( z < (1. + p) )
    notusedyet = notusedyet[0]
    if np.size(notusedyet) == 0:
        muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+ u2*etad)/omega
        mu0=1.-lambdae
        return [muo1,mu0]

    # Case 11 - the  source is completely occulted:
    if p >= 1.:
        occulted = np.where(z[notusedyet] <= p-1.)#,complement=notused2)
        if np.size(occulted) != 0:
            ndxuse = notusedyet[occulted]
            etad[ndxuse] = 0.5 # corrected typo in paper
            lambdae[ndxuse] = 1.
            # lambdad = 0 already
            notused2 = np.where(z[notusedyet] > p-1)
            if np.size(notused2) == 0:
                muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+u2*etad)/omega
                mu0=1.-lambdae
                return [muo1,mu0]
            notusedyet = notusedyet[notused2]
                
    # Case 2, 7, 8 - ingress/egress (uniform disk only)
    inegressuni = np.where((z[notusedyet] >= np.abs(1.-p)) & (z[notusedyet] < 1.+p))
    if np.size(inegressuni) != 0:
        ndxuse = notusedyet[inegressuni]
        tmp = (1.-p**2.+z[ndxuse]**2.)/2./z[ndxuse]
        tmp = np.where(tmp > 1.,1.,tmp)
        tmp = np.where(tmp < -1.,-1.,tmp)
        kap1 = np.arccos(tmp)
        tmp = (p**2.+z[ndxuse]**2-1.)/2./p/z[ndxuse]
        tmp = np.where(tmp > 1.,1.,tmp)
        tmp = np.where(tmp < -1.,-1.,tmp)
        kap0 = np.arccos(tmp)
        tmp = 4.*z[ndxuse]**2-(1.+z[ndxuse]**2-p**2)**2
        tmp = np.where(tmp < 0,0,tmp)
        lambdae[ndxuse] = (p**2*kap0+kap1 - 0.5*np.sqrt(tmp))/np.pi
        # eta_1
        etad[ndxuse] = 1./2./np.pi*(kap1+p**2*(p**2+2.*z[ndxuse]**2)*kap0-(1.+5.*p**2+z[ndxuse]**2)/4.*np.sqrt((1.-x1[ndxuse])*(x2[ndxuse]-1.)))
    
    # Case 5, 6, 7 - the edge of planet lies at origin of star
    ocltor = np.where(z[notusedyet] == p)#, complement=notused3)
    t = np.where(z[notusedyet] == p)
    if np.size(ocltor) != 0:
        ndxuse = notusedyet[ocltor] 
        if p < 0.5:
            # Case 5
            q=2.*p  # corrected typo in paper (2k -> 2p)
            Ek,Kk = ellke(q)
            # lambda_4
            lambdad[ndxuse] = 1./3.+2./9./np.pi*(4.*(2.*p**2-1.)*Ek+(1.-4.*p**2)*Kk)
            # eta_2
            etad[ndxuse] = p**2/2.*(p**2+2.*z[ndxuse]**2)
            lambdae[ndxuse] = p**2 # uniform disk
        elif p > 0.5:
            # Case 7
            q=0.5/p # corrected typo in paper (1/2k -> 1/2p)
            Ek,Kk = ellke(q)
            # lambda_3
            lambdad[ndxuse] = 1./3.+16.*p/9./np.pi*(2.*p**2-1.)*Ek-(32.*p**4-20.*p**2+3.)/9./np.pi/p*Kk
            # etad = eta_1 already
        else:
            # Case 6
            lambdad[ndxuse] = 1./3.-4./np.pi/9.
            etad[ndxuse] = 3./32.
        notused3 = np.where(z[notusedyet] != p)
        if np.size(notused3) == 0:
            muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+u2*etad)/omega
            mu0=1.-lambdae
            return [muo1,mu0]
        notusedyet = notusedyet[notused3]

    # Case 2, Case 8 - ingress/egress (with limb darkening)
    inegress = np.where( ((z[notusedyet] > 0.5+np.abs(p-0.5)) & (z[notusedyet] < 1.+p))  | ( (p > 0.5) & (z[notusedyet] > abs(1.-p)) & (z[notusedyet] < p)) )#, complement=notused4)
    if np.size(inegress) != 0:

        ndxuse = notusedyet[inegress]
        q=np.sqrt((1.-x1[ndxuse])/(x2[ndxuse]-x1[ndxuse]))
        Ek,Kk = ellke(q)
        n=1./x1[ndxuse]-1.

        # lambda_1:
        lambdad[ndxuse]=2./9./np.pi/np.sqrt(x2[ndxuse]-x1[ndxuse])*(((1.-x2[ndxuse])*(2.*x2[ndxuse]+x1[ndxuse]-3.)-3.*x3[ndxuse]*(x2[ndxuse]-2.))*Kk+(x2[ndxuse]-x1[ndxuse])*(z[ndxuse]**2+7.*p**2-4.)*Ek-3.*x3[ndxuse]/x1[ndxuse]*ellpic_bulirsch(n,q))

        notused4 = np.where( ( (z[notusedyet] <= 0.5+np.abs(p-0.5)) | (z[notusedyet] >= 1.+p) ) & ( (p <= 0.5) | (z[notusedyet] <= np.abs(1.-p)) | (z[notusedyet] >= p) ))
        if np.size(notused4) == 0:
            muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+u2*etad)/omega
            mu0=1.-lambdae
            return [muo1,mu0]
        notusedyet = notusedyet[notused4]

    # Case 3, 4, 9, 10 - planet completely inside star
    if p < 1.:
        inside = np.where(z[notusedyet] <= (1.-p))#, complement=notused5)
        if np.size(inside) != 0:
            ndxuse = notusedyet[inside]

            ## eta_2
            etad[ndxuse] = p**2/2.*(p**2+2.*z[ndxuse]**2)

            ## uniform disk
            lambdae[ndxuse] = p**2

            ## Case 4 - edge of planet hits edge of star
            edge = np.where(z[ndxuse] == 1.-p)#, complement=notused6)
            if np.size(edge[0]) != 0:
                ## lambda_5
                lambdad[ndxuse[edge]] = 2./3./np.pi*np.arccos(1.-2.*p)-4./9./np.pi*np.sqrt(p*(1.-p))*(3.+2.*p-8.*p**2)
                if p > 0.5:
                    lambdad[ndxuse[edge]] -= 2./3.
                notused6 = np.where(z[ndxuse] != 1.-p)
                if np.size(notused6) == 0:
                    muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+u2*etad)/omega
                    mu0=1.-lambdae
                    return [muo1,mu0]
                ndxuse = ndxuse[notused6[0]]

            ## Case 10 - origin of planet hits origin of star
            origin = np.where(z[ndxuse] == 0)#, complement=notused7)
            if np.size(origin) != 0:
                ## lambda_6
                lambdad[ndxuse[origin]] = -2./3.*(1.-p**2)**1.5
                notused7 = np.where(z[ndxuse] != 0)
                if np.size(notused7) == 0:
                    muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+u2*etad)/omega
                    mu0=1.-lambdae
                    return [muo1,mu0]
                ndxuse = ndxuse[notused7[0]]

            q=np.sqrt((x2[ndxuse]-x1[ndxuse])/(1.-x1[ndxuse]))
            n=x2[ndxuse]/x1[ndxuse]-1.
            Ek,Kk = ellke(q)

            ## Case 3, Case 9 - anywhere in between
            ## lambda_2
            lambdad[ndxuse] = 2./9./np.pi/np.sqrt(1.-x1[ndxuse])*((1.-5.*z[ndxuse]**2+p**2+x3[ndxuse]**2)*Kk+(1.-x1[ndxuse])*(z[ndxuse]**2+7.*p**2-4.)*Ek-3.*x3[ndxuse]/x1[ndxuse]*ellpic_bulirsch(n,q))

        ## if there are still unused elements, there's a bug in the code
        ## (please report it)
        notused5 = np.where(z[notusedyet] > (1.-p))
        if notused5[0] != 0:
            print "ERROR: the following values of z didn't fit into a case:"
            return [-1,-1]

        muo1 =1.-((1.-u1-2.*u2)*lambdae+(u1+2.*u2)*(lambdad+2./3.*(p > z))+u2*etad)/omega
        mu0=1.-lambdae
        return [muo1,mu0]

def TransitLC(timeIn,F0,inc,aRs,Period,RpRs,u1,u2,T0):
    """
        Computes transit lightcurve.
        timeIn - array of timestamps
        F0 - flat flux level
        inc - orbital inclination in radians
        aRs - a/Rstar
        Period in days
        RpRs - Rp/Rs
        u1,u2 - Limb-Darkening Coefficient
        T0 - Mid Transit Time
        
    """

    # RpRs = tqe.MTQ_getRpRs(u1,u2,tT,tG,D)
    th = (2e0*np.pi/Period)*(timeIn-T0)
    z0 = aRs*np.sqrt(1e0 - (np.cos(th)*np.sin(inc))**2e0)
    out = occultquad(z0,u1,u2,RpRs)
    flux= out[:][0]*F0
    return flux
