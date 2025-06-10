#define IOSGB iosgb
#define JCAP jcap
#define JOSGB josgb
#define LEVS levs
#define NCLDG ncldg
#define NTOTAL ntotal
#define NTRAC ntrac
#include <define.h>
      subroutine sgb1(idrt,fhour,idate,nsig,si,sl,
     $                ntrace,ncloud,
     &                nsgb,ncpus,ids,icen,icen2,igen,lhgrad)
#include <paramodel.h>
c-paramodel
c
      dimension idate(4),si(*),sl(*),ids(255),ipdsx(100)
c
      parameter(io=io_, jo=jo_, mxbit=16)
      parameter(jcap=jcap_,levs=levs_)
      parameter(nc=(jcap_+1)*(jcap_+2)+1,nctop=(jcap_+1)*2)
c
c      parameter(nflds=18*levs_+12)	
      parameter(nflds=18*levs_+(ntotal_-1)*levs_+12)	
c      parameter(mflds=13*levs_+(ntotal_-1)*levs_+6)	
c
      dimension slat((jo+1)/2),clat((jo+1)/2),wlat((jo+1)/2)
      dimension ifax(20),trig(io*2)
      dimension eps(nc/2),epstop(nctop/2)
      dimension s(nc,nflds),sstop(nctop,nflds)
      dimension fxs(2*io+6,nflds),fxys(2*io,(jo+1)/2,nflds)
c
      parameter(iput=11,ipuu=33,ipuv=34,ipuvor=43,ipudiv=44,ipuq=51)
      parameter(ipuo3=154,ipupr=152,ipucl=153) !! added for ntrace/ncloud
      parameter(ipup=1,ipuz=7,ipupx=181,ipupy=182,ipuzx=183,ipuzy=184)
c
c  added for ceop
c
      parameter(ipugz=70,ipugzi=142,ipugzx=195,ipugzy=196,ipugzd= 29,
     1          ipugzdi=110)
      parameter(ipuct=60,ipucti=187,iputvx=191,iputvy=192,ipuctd= 30,
     1          ipuctdi=189)
      parameter(         ipuqi = 54,ipuqx =193,ipuqy =194,ipuqd =200,
     1          ipuqdi =108)
      parameter(ipuomeg=39)
c
      parameter(isfc=1,isglev=107,iscol=10)
c
      dimension ipu(nflds),itl(nflds),il2(nflds)
      character grib(200+io*jo*(mxbit+1)/8)
c
      logical lhgrad,lngrad
      logical  lppr,lpcl,lpo3 !! added for ntrace/ncloud
c
      dimension mpf(255)
c
      data sigmn,smgmx/0.,1./
      data ipdsx/100*0/
c
      character*3 fni
      parameter(fni='sgb',nchi=3)
      character*80 fno
#ifdef ASSIGN
      character*80 asgnstr
#endif
c
      call fnam(fni,nchi,fhour,fno,ncho)
#ifdef ASSIGN
      write(asgnstr,'(23hassign -s unblocked  u:,I2,)') npgb
      call assign('assign -R')
      call assign(asgnstr)
#endif
      open(unit=nsgb,file=fno(1:ncho),form='unformatted',err=900)
      go to 901
  900 continue
      write(6,*) ' error in opening file ',fno(1:ncho)
      call abort
  901 continue
      write(6,*) ' file ',fno(1:ncho),' opened. unit=',nsgb 
c
      call rdss(nsig,si,sl,idrt,
     $          nflds,ntrace,ncloud,
     $          clat,slat,trig,ifax,eps,epstop,s,sstop)
c
#ifndef ISOTOPE      
      lpcl=ncloud.ge.1
      lpo3=ntrace.ge.2
      lppr=ncloud.ge.2
#else
      lpcl=.true.
      lpo3=.true.
      lppr=.true.
#endif
      jfhour=nint(fhour)
      len=io*jo
c
      ksz   = 1
      ksd   = 1+   levs
      kst   = 1+ 2*levs
      ksq   = 1+ 3*levs
      ksu   = 1+ 4*levs
      ksv   = 1+ 5*levs
      ksps  = 1+ 6*levs
      kszs  = 2+ 6*levs
c
      kspsx = 3+ 6*levs
      kspsy = 4+ 6*levs
      kszsx = 5+ 6*levs
      kszsy = 6+ 6*levs
c
      ksgz  = 7+ 6*levs
      kstvx = 7+ 7*levs
      kstvy = 7+ 8*levs
      ksqx  = 7+ 9*levs
      ksqy  = 7+10*levs
      ksgzx = 7+11*levs
      ksgzy = 7+12*levs
c
      ksct  = 7+13*levs
      kscti = 7+14*levs
      ksctd = 8+14*levs
      ksctdi= 8+15*levs
      ksqi  = 9+15*levs
      ksqd  =10+15*levs
      ksqdi =10+16*levs
      ksgzi =11+16*levs
      ksgzd =12+16*levs
      ksgzdi=12+17*levs
c
      ksomeg=13+17*levs
c
      if (lpo3) then
         kso3=ksomeg+levs
      else
         kso3=ksomeg
      endif
      if (lpcl) then
         kscl=kso3+levs
      else
         kscl=kso3
      endif
      if (lppr) then
         kspr=kscl+levs
      else
         kspr=kscl
      endif
c
c  1.  definition of grib header rec no. 5
c
c  1.1 dependent variables
c
      do i=ksu,ksu+levs-1
        ipu(i)=ipuu
      enddo
      do i=ksv,ksv+levs-1
        ipu(i)=ipuv
      enddo
      do i=ksz,ksz+levs-1
        ipu(i)=ipuvor
      enddo
      do i=ksd,ksd+levs-1
        ipu(i)=ipudiv
      enddo
      do i=kst,kst+levs-1
        ipu(i)=iput
      enddo
      do i=ksq,ksq+levs-1
        ipu(i)=ipuq
      enddo
      ipu(ksps )=ipup
      ipu(kszs )=ipuz
c
c  1.2 horizontal derivatives of dependent variables
c
      ipu(kspsx)=ipupx
      ipu(kspsy)=ipupy
      ipu(kszsx)=ipuzx
      ipu(kszsy)=ipuzy
c
c  1.3 enthalpy and potential energy
c
      do i=ksct,ksct+levs-1
        ipu(i)=ipuct
      enddo    
      do i=ksgz,ksgz+levs-1
        ipu(i)=ipugz
      enddo    
c
c  1.4 divergence of enthalpy, potential energy and moisture
c
      do i=ksctd,ksctd+levs-1
        ipu(i)=ipuctd
      enddo    
      do i=ksqd,ksqd+levs-1
        ipu(i)=ipuqd
      enddo    
      do i=ksgzd,ksgzd+levs-1
        ipu(i)=ipugzd
      enddo    
c
c  1.5 horiontal derivatives of enthalpy, potential energy and moisture
c
      do i=kstvx,kstvx+levs-1
        ipu(i)=iputvx
      enddo    
      do i=kstvy,kstvy+levs-1
        ipu(i)=iputvy
      enddo    
      do i=ksqx,ksqx+levs-1
        ipu(i)=ipuqx
      enddo    
      do i=ksqy,ksqy+levs-1
        ipu(i)=ipuqy
      enddo    
      do i=ksgzx,ksgzx+levs-1
        ipu(i)=ipugzx
      enddo    
      do i=ksgzy,ksgzy+levs-1
        ipu(i)=ipugzy
      enddo    
c
c  1.6 vertical integral of enthalpy, potential energy and moisture
c
      ipu(kscti)=ipucti
      ipu(ksqi )=ipuqi
      ipu(ksgzi)=ipugzi
c
c  1.7 vertical integral of the divergence of enthalpy, potential energy and moisture
c
      ipu(ksctdi)=ipuctdi
      ipu(ksqdi )=ipuqdi
      ipu(ksgzdi)=ipugzdi
c
c  1.8 omega
c
      do i=ksomeg,ksomeg+levs-1
        ipu(i)=ipuomeg
      enddo    
c
c  1.9 ntrace/ncloud
c
      if (lpo3) then
      do i=kso3,kso3+levs-1
         ipu(i)=ipuo3
      enddo
      endif
      if (lpcl) then
      do i=kscl,kscl+levs-1
         ipu(i)=ipucl
      enddo
      endif
      if (lppr) then
      do i=kspr,kspr+levs-1
         ipu(i)=ipupr
      enddo
      endif
c      write(6,*) 'check 1.9'
c
c  2.  definition of the number of levels
c
c  2.1 dependent variables
c
      do i=ksu,ksu+levs-1
        itl(i)=isglev
      enddo
      do i=ksv,ksv+levs-1
        itl(i)=isglev
      enddo
      do i=ksz,ksz+levs-1
        itl(i)=isglev
      enddo
      do i=ksd,ksd+levs-1
        itl(i)=isglev
      enddo
      do i=kst,kst+levs-1
        itl(i)=isglev
      enddo
      do i=ksq,ksq+levs-1
        itl(i)=isglev
      enddo
      itl(ksps)=isfc
      itl(kszs)=isfc
c
c  2.2 horizontal derivatives of dependent variables
c
      itl(kspsx)=isfc
      itl(kspsy)=isfc
      itl(kszsx)=isfc
      itl(kszsy)=isfc
c
c  2.3 enthalpy and potential energy
c
      do i=ksct,ksct+levs-1
        itl(i)=isglev
      enddo
      do i=ksgz,ksgz+levs-1
        itl(i)=isglev
      enddo
c
c  2.4 divergence of enthalpy, potential energy and moisture
c
      do i=ksctd,ksctd+levs-1
        itl(i)=isglev
      enddo
      do i=ksqd,ksqd+levs-1
        itl(i)=isglev
      enddo
      do i=ksgzd,ksgzd+levs-1
        itl(i)=isglev
      enddo
c
c  2.5 horiontal derivatives of enthalpy, potential energy and moisture
c
      do i=kstvx,kstvx+levs-1
        itl(i)=isglev
      enddo
      do i=kstvy,kstvy+levs-1
        itl(i)=isglev
      enddo
      do i=ksqx,ksqx+levs-1
        itl(i)=isglev
      enddo
      do i=ksqy,ksqy+levs-1
        itl(i)=isglev
      enddo
      do i=ksgzx,ksgzx+levs-1
        itl(i)=isglev
      enddo
      do i=ksgzy,ksgzy+levs-1
        itl(i)=isglev
      enddo
c
c  2.6 vertical integral of enthalpy, potential energy and moisture
c
      itl(kscti)=iscol
      itl(ksqi )=iscol
      itl(ksgzi)=iscol
c
c  2.7 vertical integral of the divergence of enthalpy, potential energy and moisture
c
      itl(ksctdi)=iscol
      itl(ksqdi )=iscol
      itl(ksgzdi)=iscol
c
c  2.8 omega
c
      do i=ksomeg,ksomeg+levs-1
        itl(i)=isglev
      enddo
c
c  2.9 ntrace/ncloud
c
      if (lpo3) then
      do i=kso3,kso3+levs-1
        itl(i)=isglev
      enddo
      endif
      if (lpcl) then
      do i=kscl,kscl+levs-1
        itl(i)=isglev
      enddo
      endif
      if (lppr) then
      do i=kspr,kspr+levs-1
        itl(i)=isglev
      enddo
      endif
c      write(6,*) 'check 2.9'
c
c  3. definition of sigma levels
c
c  3.1 dependent variables
c
      do i=ksu,ksu+levs-1
        il2(i)=nint(sl(i-ksu+1)*1.e4)
      enddo
      do i=ksv,ksv+levs-1
        il2(i)=nint(sl(i-ksv+1)*1.e4)
      enddo
      do i=ksz,ksz+levs-1
        il2(i)=nint(sl(i-ksz+1)*1.e4)
      enddo
      do i=ksd,ksd+levs-1
        il2(i)=nint(sl(i-ksd+1)*1.e4)
      enddo
      do i=kst,kst+levs-1
        il2(i)=nint(sl(i-kst+1)*1.e4)
      enddo
      do i=ksq,ksq+levs-1
        il2(i)=nint(sl(i-ksq+1)*1.e4)
      enddo
      il2(ksps)=0
      il2(kszs)=0
c
c  3.2 horizontal derivatives of dependent variables
c
      il2(kspsx)=0
      il2(kspsy)=0
      il2(kszsx)=0
      il2(kszsy)=0
c
c  3.3 enthalpy and potential energy
c
      do i=ksct,ksct+levs-1
        il2(i)=nint(sl(i-ksct+1)*1.e4)
      enddo
      do i=ksgz,ksgz+levs-1
        il2(i)=nint(sl(i-ksgz+1)*1.e4)
      enddo
c
c  3.4 divergence of enthalpy, potential energy and moisture
c
      do i=ksctd,ksctd+levs-1
        il2(i)=nint(sl(i-ksctd+1)*1.e4)
      enddo
      do i=ksqd,ksqd+levs-1
        il2(i)=nint(sl(i-ksqd+1)*1.e4)
      enddo
      do i=ksgzd,ksgzd+levs-1
        il2(i)=nint(sl(i-ksgzd+1)*1.e4)
      enddo
c
c  3.5 horiontal derivatives of enthalpy, potential energy and moisture
c
      do i=kstvx,kstvx+levs-1
        il2(i)=nint(sl(i-kstvx)*1.e4)
      enddo
      do i=kstvy,kstvy+levs-1
        il2(i)=nint(sl(i-kstvy+1)*1.e4)
      enddo
      do i=ksqx,ksqx+levs-1
        il2(i)=nint(sl(i-ksqx+1)*1.e4)
      enddo
      do i=ksqy,ksqy+levs-1
        il2(i)=nint(sl(i-ksqy+1)*1.e4)
      enddo
      do i=ksgzx,ksgzx+levs-1
        il2(i)=nint(sl(i-ksgzx+1)*1.e4)
      enddo
      do i=ksgzy,ksgzy+levs-1
        il2(i)=nint(sl(i-ksgzy+1)*1.e4)
      enddo
c
c  3.6 vertical integral of enthalpy, potential energy and moisture
c
      il2(kscti)=0
      il2(ksqi )=0
      il2(ksgzi)=0
c
c  3.7 vertical integral of the divergence of enthalpy, potential energy and moisture
c
      il2(ksctdi)=0
      il2(ksqdi )=0
      il2(ksgzdi)=0
c
c  3.8 omega
c
      do i=ksomeg,ksomeg+levs-1
        il2(i)=nint(sl(i-ksomeg+1)*1.e4)
      enddo
c
c  3.9 ntrace/ncloud
c
      if (lpo3) then
      do i=kso3,kso3+levs-1
        il2(i)=nint(sl(i-kso3+1)*1.e4)
      enddo
      endif
      if (lpcl) then
      do i=kscl,kscl+levs-1
        il2(i)=nint(sl(i-kscl+1)*1.e4)
      enddo
      endif
      if (lppr) then
      do i=kspr,kspr+levs-1
        il2(i)=nint(sl(i-kspr+1)*1.e4)
      enddo
      endif
c
      il2mn=nint(sigmn*1.e4)
      il2mx=nint(sigmx*1.e4)
c
c      write(6,*) 'check 3.9'
      call mpfdef(2,mpf)
c
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(jo,jcap,nc,nctop,levs,trig,ifax,eps,epstop,s,sstop)
cmic$1 shared(io,clat,slat,ksps,fxys,nflds)
cmic$1 private(i,j,k,fxs)
#endif
      do j=1,(jo+1)/2
        call trss(trig,ifax,eps,epstop,s,sstop,clat(j),slat(j),
     $            nflds,
     1            fxs,j,si,sl)
c
        do i=1,io*2
          fxs(i,ksps)=1.e3*fxs(i,ksps)
          do k=1,nflds
            fxys(i,j,k)=fxs(i,k)
          enddo
        enddo
      enddo
c
      do k=1,nflds
c
c  skip gribbing of x- and y- derivatives if desired
c
        lngrad=.true.
        if(ipu(k).ge.181.and.ipu(k).le.184.or.
     1     ipu(k).ge.191.and.ipu(k).le.196) lngrad=.false.
        if(lhgrad) lngrad=.true.
        if((.not.lpo3).and.k.ge.kso3.and.k.lt.kso3+levs) lngrad=.false.
        if((.not.lpcl).and.k.ge.kscl.and.k.lt.kscl+levs) lngrad=.false.
        if((.not.lppr).and.k.ge.kspr.and.k.lt.kspr+levs) lngrad=.false.
c
        if(lngrad) then
          call polext(mpf(ipu(k)),io,fxys(1,2,k),fxys(1+io,2,k),
     &                fxys(1,1,k),fxys(1+io,1,k))
          call rowsep(fxys(1,1,k))
c
          fmax=fxys(1,1,k)
          fmin=fxys(1,1,k)
          do j=1,(jo+1)/2
            do i=1,io*2
              if(fmax.lt.fxys(i,j,k)) fmax=fxys(i,j,k)
              if(fmin.gt.fxys(i,j,k)) fmin=fxys(i,j,k)
            enddo
          enddo
          print 100,k,fmax,fmin
  100     format('k=',i5,' max=',e14.5,' min=',e14.5)
c
          call gribit(fxys(1,1,k),ibmap,idrt,io,jo,mxbit,acos(slat(1)),
     &                28,2,icen,igen,0,ipu(k),itl(k),0,il2(k),
     &                idate(4),idate(2),idate(3),idate(1),1,jfhour,0,10,
     &                0,0,icen2,ids(ipu(k)),ipdsx,
     &                0.,0.,0.,0.,0.,0.,0.,0.,
     &                grib,lgrib,ierr)
          if(ierr.eq.0) call wryte(51,lgrib,grib)
        endif
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
