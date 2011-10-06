
import numpy as num
import pylab
import uwpyKepler as kep
import sys

file = open('QATSTEST.txt','r')
file = file.readlines()
tflat = []
fflat = []
flagflat = []

for line in file:
    split = map(str,line.split('|'))
    fflat.append(float(split[0]))
    flagflat.append(float(split[1]))
    tflat.append(float(split[2]))

fflat = num.array(fflat)
flagflat = num.array(flagflat)
tflat = num.array(tflat)
extranoise = 0
MEDGAP = 0.020434495
N0 = -1
NTRANS = 41
PDOT = 0.0000000
PERIOD = 3.5222970
T0 = 54954.118
TD = 0.13691667

#pylab.plot(tflat,fflat,'b.')
#pylab.show()
#pylab.plot(tflat,flagflat,'b.')
#pylab.show()

#sys.exit()

nt = len(tflat)
ndata = 3
data = num.zeros( (ndata,nt) )
data[0][:] = 1e0+fflat
data[1][:] = 1e0+fflat*(flagflat == 0)
sigma = num.std(data[1][num.where(data[1][:] != 0e0)[0]])
data[2][:] = 1e0+sigma*num.random.randn(nt)*(num.sum(flagflat) == 0)
# Add additional noise where there are no gaps:
if (extranoise != 0e0):
    iseg = num.where(fflat != 0e0)[0]
    nadd = len(iseg)
    addnoise = num.zeros(nt)
    addnoise[iseg] = extranoise*num.random.randn(nadd)
    for idata in range(ndata):
        inon=(num.where(data[idata][:] != 0e0)[0]).ravel()
        data[idata,inon]=data[idata,inon]+addnoise[inon]

gap = num.median(tflat[1L:nt]-tflat[0L:nt-1])

# Now call detection routine:
# p_i  = p_i-1 *(1+f/2)/(1-f/2)  from 100 -10,000# f=0.2# round tmin & tmax
# Search from 5 days to 150 days:
pmin = long(num.floor(1.3e0/gap))
#pmin = long(num.floor(60e0/gap))
pmax = long(num.ceil(100e0/gap))
print 'Range of periods: ',pmin,pmax
#pmin=num.floor(5d0/gap)
#pmax=num.ceil(150d0/gap)
#pmin=646.88611d0
#pmax=646.88611d0
#pmin=1647.7602d0
#pmax=1647.7602d0
f = 0.005e0
nperiod = long(num.log(pmax/pmin)/num.log((1+f/2)/(1-f/2)))
#print long(nperiod), pmax,pmin
#nperiod=2
#nperiod=1
speriod = num.zeros((ndata,nperiod))
#print num.shape(speriod)
#sys.exit()

period0 = pmin/(1+f/2)*(1-f/2)
period = num.zeros(nperiod)
#period=[pmin,pmax]
#period=[pmin]
nhatbest = num.zeros( (ndata,nperiod,500),dtype='int64')
mbest = num.zeros( (ndata,nperiod) )
qbest = num.zeros( (ndata,nperiod) )

#for ip in range(nperiod):
for ip in [0]:
    spmax = num.zeros(ndata)
    period0 = (period0*(1+f/2)/(1-f/2))
    #  period0=period[ip]
    period[ip] = period0
    tmin = num.floor(period0*(1-f/2))
    tmax = num.ceil(period0*(1+f/2))
    q = num.floor(8e0*(float(period0)/600e0)**(1./3.))
    #print ip,period0
    #for idata in range(ndata):
    for idata in [0]:
        datatmp=(data[idata][:].ravel())
        MM,nhat,smax,dc = kep.qats.qpt_detect(datatmp,tmin,tmax,q)
        print MM,num.shape(nhat),smax,num.shape(dc)
        speriod[idata,ip]=smax
        if(smax > spmax[idata]):
            mbest[idata][ip]=MM
            qbest[idata][ip]=q
            nhatbest[idata][ip][:]=0
            nhatbest[idata][ip][0:MM]=nhat
            spmax[idata]=smax


#ax = pylab.subplot(1,1,1)
#snr=speriod/num.sqrt(mbest*qbest)
#print num.shape(snr)
#pylab.plot(period,snr[0][:],'b-')
#ax.set_xscale('log')
#pylab.show()

#for idata in range(ndata):
    #ibest = num.where(snr[idata][:] == max(snr[idata][:]))[0]
    #pbest = period[ibest]
    #mm = mbest[idata][ibest]
    #coeff = num.polyfit(num.arange(mm),tflat[nhatbest[idata][ibest][0:mm]],1)
    #pylab.plot(tflat[nhatbest[idata][ibest][1:mm]]-coeff[0]-coeff[1]*arange(mm),'b.')
    ##print idata,pbest,mm

#print 'Finished'

#import uwpyKepler as kep

#MM = 10
#N = 10
#tmin = 1
#tmax = 50
#q = 0.4
#dc = 10
#nhat = 10

#print kep.qats.fmax(MM,N,tmin,tmax,q,dc,nhat)