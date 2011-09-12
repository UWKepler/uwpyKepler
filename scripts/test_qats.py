pro test_qpt_koi010,extranoise
restore,'koi010_detrend.sav'
nt=n_elements(tflat)
ndata=3
data=dblarr(ndata,nt)
data[0,*]=1d0+fflat
data[1,*]=1d0+fflat*(flagflat[0,*] eq 0)
sigma=stddev(data[1,where(data[1,*] ne 0d0)])
data[2,*]=1d0+sigma*randomn(seed,nt)*(total(flagflat,1) eq 0)
; Add additional noise where there are no gaps:
if(extranoise ne 0d0) then begin
  iseg=where(fflat ne 0d0)
  nadd=n_elements(iseg)
  addnoise=dblarr(nt)
  addnoise[iseg]=extranoise*randomn(seed,nadd)
  for idata=0,ndata-1 do begin
    inon=reform(where(data[idata,*] ne 0d0))
    data[idata,inon]=data[idata,inon]+addnoise[inon]
  endfor
endif

gap=median(tflat[1L:nt-1L]-tflat[0L:nt-2L])

; Now call detection routine:
; p_i  = p_i-1 *(1+f/2)/(1-f/2)  from 100 -10,000; f=0.2; round tmin & tmax
; Search from 5 days to 150 days:
pmin=floor(1.3d0/gap)
pmax=ceil(50d0/gap*2d0)
print,'Range of periods: ',pmin,pmax
;pmin=floor(5d0/gap)
;pmax=ceil(150d0/gap)
;pmin=646.88611d0
;pmax=646.88611d0
;pmin=1647.7602d0
;pmax=1647.7602d0
f=.005d0
nperiod=alog(pmax/pmin)/alog((1+f/2)/(1-f/2))
;nperiod=2
;nperiod=1
speriod=dblarr(ndata,nperiod)
period0=pmin/(1+f/2)*(1-f/2)
period=dblarr(nperiod)
;period=[pmin,pmax]
;period=[pmin]
nhatbest=intarr(ndata,nperiod,500)
mbest=dblarr(ndata,nperiod)
qbest=dblarr(ndata,nperiod)
for ip=0,nperiod-1 do begin
  spmax=dblarr(ndata)
  period0=period0*(1+f/2)/(1-f/2)
;  period0=period[ip]
  period[ip]=period0
  tmin=floor(period0*(1-f/2))
  tmax=ceil(period0*(1+f/2))
  q=floor(8d0*(double(period0)/600d0)^(1./3.))
;  q=floor(4d0*(double(period0)/600d0)^(1./3.))
  print,ip,period0
  for idata=0,ndata-1 do begin
    datatmp=reform(data[idata,*])
    qpt_detect,datatmp,MM,tmin,tmax,q,dc,nhat,smax
    speriod[idata,ip]=smax
    if(smax gt spmax[idata]) then begin
       mbest[idata,ip]=MM
       qbest[idata,ip]=q
       nhatbest[idata,ip,*]=0
       nhatbest[idata,ip,0:MM-1]=nhat
       spmax[idata]=smax
    endif
  endfor
endfor
snr=speriod/sqrt(mbest*qbest)
plot,period,snr[0,*],/xl
oplot,period,snr[1,*],col=255L
oplot,period,snr[2,*],col=255L*256L
;oplot,period,snr[3,*],col=255L*256L+255L
;2oplot,period,snr[4,*],col=255L*256L+255L*256L^2
;p1=13.197120
;alias=[2./5.,1./2.,2./3.,3./4.,4./5.,5./6.,6./7.,7./8.,1.,7./6.,6./5.,5./4.,4./3.,3./2.,2.,3] & nalias=n_elements(alias)
;for ialias=0,nalias-1 do oplot,p1/gap*[1.,1.]*alias[ialias],[.1,100],col=255
;p2=33.653d0
;p2=50.45d0
;alias=[.5,1.,2.] & nalias=n_elements(alias)
;alias=[1./3.,.5,1.,2.,3.] & nalias=n_elements(alias)
;for ialias=0,nalias-1 do oplot,p2/gap*[1.,1.]*alias[ialias],[.1,100],col=255L*256L
c=get_kbrd(1)
for idata=0,ndata-1 do begin &$
  ibest=where(snr[idata,*] eq max(snr[idata,*])) &$
  pbest=period[ibest] &$
  mm=mbest[idata,ibest] &$
  coeff=poly_fit(dindgen(mm),tflat[nhatbest[idata,ibest,0:mm-1]],1,/double) &$
;  plot,nhatbest[idata,ibest,1:mm-1]-nhatbest[idata,ibest,0:mm-2],ys=1 &$
  plot,tflat[nhatbest[idata,ibest,1:mm-1]]-coeff[0]-coeff[1]*dindgen(mm),ys=1,psym=8,syms=3 &$
  print,idata,pbest,mm &$
  c=get_kbrd(1)&$
endfor
print,'Finished'
c=get_kbrd(1)
return
e8nd