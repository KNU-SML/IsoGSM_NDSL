#define G g
#define HVAP hvap
#define IJCAP ijcap
#define ILATG ilatg
#define ILEVS ilevs
#define ILONF ilonf
#define JCAP jcap
#define LATG latg
#define LDIAG ldiag
#define LEVS levs
#define LGRID lgrid
#define LONF lonf
#define LPNT lpnt
#define LSOIL lsoil
#define LTSTP ltstp
#define LVAR lvar
#define PI pi
#define PSAT psat
#define RD rd
#define RV rv
#define TTP ttp
#define CP cp
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine extrap(lnew)
#include <paramodel.h>
      save
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    extrap      vertically interpolate spectral data.
c   prgmmr: mark iredell     org: w/nmc23    date: 91-03-15
c
c abstract: this routine transforms spectral data to gridspace,
c           fixes moisture in the stratosphere if necessary,
c           interpolates surface pressure to the new orography,
c           vertically interpolates wind, temperature and moisture,
c           and transforms the data back to spectral space.
c
c program history log:
c   91-03-15  mark iredell  docblock written (prehistorical program)
c
c usage:    call extrap(lnew)
c   input argument list:
c     lnew     - logical .true. to use newer moisture climatology
c
c   subprograms called:
c     dztouv   - compute winds from vorticity and divergence
c     pln2t    - compute legendre polynomials
c     sums2a   - transform scalar spectral to fourier
c     sumv2a   - transform vector spectral to fourier
c     fftidim  - fourier transform
c     exmois   - older moisture extrapolation
c     exmoisn  - newer moisture extrapolation
c     newps    - interpolate surface pressure
c     sg2sg    - interpolate upper air quantities
c     symasa   - combine fourier coefficients
c     fl22a    - transform scalar fourier to spectral
c     gozrmt   - compute legendre derivatives
c     psum2a   - transform vector fourier to spectral vorticity
c     msum2a   - transform vector fourier to spectral divergence
c
c attributes:
c   language: fortran 77.
c   machine:
c
c$$$
c
#include <paramter.h>
#include <parmchgr.h>
c
      parameter(idimt=idim*2,jdimhf=jdim/2,
     1          mdim=(mwave+1)*(mwave+2),
     3          mdimv=(mwave+1)*(mwave+4),
     4          mwvp2=(mwave+1)*2,
     5          mdimhf=mdim/2,mdimvh=mdimv/2,
     6          kdimp=kdim+1,kdimm=kdim-1,
     7          kdimt=kdim*2)
c
c
c hmhj 2        u v t q ps
      parameter(mkdi=(kdimi*3+kdimqi+1)*2)
      parameter(mkd =(kdim *3+kdimq +1)*2)
c
      parameter(idimti=idimi*2,jdimhfi=jdimi/2,ijdimi=idimi*jdimi,
     1          mdimi=(mwavei+1)*(mwavei+2),
     3          mdimvi=(mwavei+1)*(mwavei+4),
     4          mwvp2i=(mwavei+1)*2,
     5          mdimhfi=mdimi/2,mdimvhi=mdimvi/2,
     6          kdimpi=kdimi+1,kdimmi=kdimi-1)
c
#include <comchgr.h>
c
c hmhj 3
      dimension fj(idimt,kdimi),gj(idimt,kdimi),
     1          work1(idim,mkdi),work2(idim,mkdi),
     2          work3(idim,mkd)
c hmhj 2
      common /inpver/ siin(kdimpi),slin(kdimi),delin(kdimi),
     1                ciin(kdimpi),clin(kdimi),rpiin(kdimmi)
     2               ,ak(kdimpi),bk(kdimpi)
c
c hmhj 5
      common pss(idimt),tts(idimt,kdim),qqs(idimt,kdimq),
     1                  uus(idimt,kdim),vvs(idimt,kdim)
      common psa(idimt),tta(idimt,kdim),qqa(idimt,kdimq),
     1                  uua(idimt,kdim),vva(idimt,kdim)
c
      common/plncom/eps   (mdimv ),
     1              colrad(jdimhf),wgt(jdimhf ),
     2              wgtcs (jdimhf),rcs2(jdimhf)
c hmhj1
      common/comver/ si  (kdimp ),sl  (kdim ),del  (kdim ),
     1               ci  (kdimp ),cl  (kdim ),rpi  (kdimm),
     $               ak5 (kdimp ),bk5 (kdimp)
c
      dimension qlnt  (mdim),qlnv(mdimv),
     1          qdert (mdim),qlnwct(mdim),
     2          dplnxn(mdim),dplnyn(mdim),
     3          plnw(mdim),plnd2w(mdim),dplnxw(mdim),dplnyw(mdim)
c
      dimension uln(mdimv,kdimi),vln(mdimv,kdimi)
      dimension zs(idimt),zso(idimt)
c
      dimension expps(idimt),rqfx(idimt,kdimi),rqfy(idimt,kdim)
      dimension exppso(idimt)
c
      dimension spin(idimt,kdimi+1),spln(idimt,kdimi)
      dimension spino(idimt,kdim+1),splno(idimt,kdim)
c
      logical lnew
#include <funvap.h>
c
c hmhj 2
      in=0
      tensn=10. 0
c      tensn=50. 0 !! linear
c transfer from slin to sl
*     print *,' ENTER EXTRAP transfer from slin to sl '
*     print *,' slin ',slin
*     print *,' sl   ',sl
c
      do 10 k=1,kdimq
      do 10 i = 1, mdim
      rqo(i,k)=0.0
   10 continue
      do 11 i = 1, mdim
      qo(i)=0.0
 11   continue
      do 12 k=1,kdim
      do 12 i = 1, mdim
      teo(i,k)=0.0
      dio(i,k)=0.0
      zeo(i,k)=0.0
   12 continue
      do 121 k=1,kdimi
      call dztouv(di(1,k),ze(1,k),uln(1,k),vln(1,k),eps)
121   continue
c
c  lat loop
c
      do 1000 lat = 1,jdimhf
c
      call pln2t(qlnt,qlnv,colrad,lat)
c
      xlat=pi_/2.0-colrad(lat)
c hmhj
      call sums2a(gz,zs,qlnt,     1)
      call fftidim(zs,zs,2,1)
c
      call sums2a(gzo,zso,qlnt,     1)
      call fftidim(zso,zso,2,1)
c
      call sums2a(q ,ps ,qlnt,     1)
      call sums2a(te,tf ,qlnt, kdimi)
      call sums2a(rq,rqf,qlnt,kdimqi)
c
      call sumv2a(uln,uf,qlnv,kdimi)
      call sumv2a(vln,vf,qlnv,kdimi)
c
      call fftidim(ps,ps,mkdi,1)
c
c  at this point grid, values are available
c
c       tf ... temperature
c      rqf ... specific humidity
c       ps ... pai
c       uf ... wind at x direction
c       vf ... wind at y direction
c
 9999 format(2x,10e12.5)
c
      sinlat=cos(colrad(lat))
      rcl=rcs2(lat)
c
      do 8745 i = 1, idimt
      expps(i)=exp(ps(i))
 8745 continue
c
      rmaxsiin=siin(1)
      do k=1,kdimpi
        rmaxsiin=max(rmaxsiin,siin(k))
      enddo
c
      if(rmaxsiin.gt.1. or. rmaxsiin.eq.0. )then
        do k=1,kdimi+1   
          do i=1,idimt  
            spin(i,k) = ak(k)/expps(i) + bk(k) 
          enddo
        enddo
        rk  = rd_/cp_
        rk1 = rk + 1.
        rkinv=1./rk
        do i=1,idimt
          do k=1,kdimi
            dif = spin(i,k)**rk1 - spin(i,k+1)**rk1
            dif = dif / (rk1*(spin(i,k)-spin(i,k+1)))
            spln(i,k) = dif**rkinv
c            spln(i,k) = (spin(i,k)+spin(i,k+1))*0.5
          enddo
        enddo
      else
        do k=1,kdimi+1
          do i=1,idimt
            spin(i,k) = siin(k)
          enddo
        enddo
        do k=1,kdimi
          do i=1,idimt
            spln(i,k) = slin(k)
          enddo
        enddo
      endif
c
c  extrapolate humidity into stratosphere
c kdimqi --> kdimi for rqf --> rqfx
c
       if(lnew) then
        call exmoisn(expps,rqf,tf,rqfx,lat,xlat)
       else
        call exmois (expps,rqf,tf,rqfx)
       endif
c
c ---   end of grid point and extrap q to stratosphere at old lev
c
c     if (lat .eq. 1 ) then
c      print *,' check grid point values before sg2sg.  lat=', lat
c      call maxmn(uf,idimt,1,kdimi,' uf ')
c      call maxmn(vf,idimt,1,kdimi,' vf ')
c      call maxmn(tf,idimt,1,kdimi,' tf ')
c      call maxmn(rqf,idimt,1,kdimqi,' rqf')
c      call maxmn(rqfx,idimt,1,kdimi,' rqfx')
c      call maxmn(ps,idimt,1,1,' ps ')
c     endif
c
c hmhj nov. 90
      call newps(tf,ps,zs,zso,pso,spin)
c     do 300 i = 1, idimt
c     pso(i)=ps(i)
c300  continue
c
      do i = 1, idimt
      exppso(i)=exp(pso(i))
      enddo
c
#ifndef HYBRID
      do k=1,kdimp
        do i=1,idimt
          spino(i,k)=si(k)
        enddo
      enddo
      do k=1,kdim
        do i=1,idimt
          splno(i,k)=sl(k)
        enddo
      enddo
#else
      rk  = rd_/cp_
      rk1 = rk + 1.
      rkinv=1./rk
c
      do k=1,kdimp
        do i=1,idimt
          spino(i,k) = ak5(k)/exppso(i) + bk5(k) 
        enddo
      enddo
      do k=1,kdim
        do i=1,idimt
          dif = spino(i,k)**rk1 - spino(i,k+1)**rk1
          dif = dif / (rk1*(spino(i,k)-spino(i,k+1)))
          splno(i,k) = dif**rkinv
c          splno(i,k) = (spino(i,k)+spino(i,k+1))*0.5
        enddo
      enddo
#endif
c
c start sg2sg by hmhj
c
c
      call sg2sg(ps,spln,uf,pso,splno,ufo,idimt,1,kdimi,kdim,in,tensn,1)
      call sg2sg(ps,spln,vf,pso,splno,vfo,idimt,1,kdimi,kdim,in,tensn,1)
      call sg2sg(ps,spln,tf,pso,splno,tfo,idimt,1,kdimi,kdim,in,tensn,1)
      call sg2sg(ps,spln,rqfx,pso,splno,rqfy,idimt,1,
     $     kdimi,kdim,in,tensn,1)
      if (ntotal_.ge.2) then
        do n=2,ntotal_
          call sg2sg(ps,spln,rqf(1,kdimi*(n-1)+1),pso,splno,
     $          rqfo(1,kdim*(n-1)+1),idimt,1,kdimi,kdim,in,tensn,1)
        enddo
      endif
c
c sg2sg sets values to constant outside of input domain.
c fix lapse rate and relative humidity below input surface.
c
      gamma=6.5e-3
      do 299 i=1,idimt
      ps1=log(spln(i,1))+ps(i)
      rh1=rqfx(i,1)/fqs(tf(i,1),exp(ps1))
      do 298 k=1,kdim
      psk=log(splno(i,k))+pso(i)
      if(psk.lt.ps1) goto 299
      tfo(i,k)=tf(i,1)*exp(gamma*rd_/g_*(psk-ps1))
      rqfy(i,k)=rh1*fqs(tfo(i,k),exp(psk))
c      print*,i,k,exp(psk),exp(ps1),tfo(i,k),tf(i,1)
298   continue
299   continue
c
c kdim --> kdimq
c
c      do 301 kq=1,kdimq
      do 301 kq=1,kdim
      do 301 i = 1, idimt
      rqfo(i,kq)=rqfy(i,kq)
 301  continue
c
        print *,'zs:(zs(1))=',zs(1)
        print *,'zso:(zso(1))=',zso(1)
        print *,'ps:(ps(1))=',expps(1)*10.
        print *,'ps:(pso(1))=',exppso(1)*10.
        print *,'spin:(spin(1,k)=',
     $       (spin(1,k),k=1,kdimi)
        print *,'spino:(spino(1,k)=',
     $       (spino(1,k),k=1,kdim)
        print *,'tf:(tf(1,k),k=1,kdimi)=',
     $       (tf(1,k),k=1,kdimi)
        print *,'tfo:(tfo(1,k),k=1,kdim)=',
     $       (tfo(1,k),k=1,kdim)
        print *,'rqf:(rqf(1,k),k=1,kdimi)=',
     $       (rqf(1,k),k=1,kdimi)
        print *,'rqfo:(rqfo(1,k),k=1,kdim)=',
     $       (rqfo(1,k),k=1,kdim)
        print *,'uf:(uf(1,k),k=1,kdimi)=',
     $       (uf(1,k),k=1,kdimi)
        print *,'ufo:(ufo(1,k),k=1,kdim)=',
     $       (ufo(1,k),k=1,kdim)
        print *,'vf:(vf(1,k),k=1,kdimi)=',
     $       (vf(1,k),k=1,kdimi)
        print *,'vfo:(vfo(1,k),k=1,kdim)=',
     $       (vfo(1,k),k=1,kdim)
        
        do i=1,idimt
          tmpa=0.
          do k=1,kdim
            tmpa=tmpa+tf(i,k)-tfo(i,k)
          enddo
          if (abs(tmpa).gt.300.) then
        print *,'Abnormal!!',i,lat,tmpa
        print *,'zs:(zs(i))=',zs(i)
        print *,'zso:(zso(i))=',zso(i)
        print *,'ps:(ps(i))=',expps(i)*10.
        print *,'ps:(pso(i))=',exppso(i)*10.
        print *,'spin:(spin(i,k)=',
     $       (spin(i,k),k=1,kdimi)
        print *,'spino:(spino(i,k)=',
     $       (spino(i,k),k=1,kdim)
        print *,'tf:(tf(i,k),k=1,kdimi)=',
     $       (tf(i,k),k=1,kdimi)
        print *,'tfo:(tfo(i,k),k=1,kdim)=',
     $       (tfo(i,k),k=1,kdim)
        print *,'rqf:(rqf(i,k),k=1,kdimi)=',
     $       (rqf(i,k),k=1,kdimi)
        print *,'rqfo:(rqfo(i,k),k=1,kdim)=',
     $       (rqfo(i,k),k=1,kdim)
        print *,'uf:(uf(i,k),k=1,kdimi)=',
     $       (uf(i,k),k=1,kdimi)
        print *,'ufo:(ufo(i,k),k=1,kdim)=',
     $       (ufo(i,k),k=1,kdim)
        print *,'vf:(vf(i,k),k=1,kdimi)=',
     $       (vf(i,k),k=1,kdimi)
        print *,'vfo:(vfo(i,k),k=1,kdim)=',
     $       (vfo(i,k),k=1,kdim)
          endif
        enddo
c
c ======= check values ========
c
c     if(lat .eq. 1 ) then
c      print *,' check grid point values after sg2sg.   lat=', lat
c      call maxmn(ufo,idimt,1,kdim,' ufo')
c      call maxmn(vfo,idimt,1,kdim,' vfo')
c      call maxmn(tfo,idimt,1,kdim,' tfo')
c      call maxmn(rqfo,idimt,1,kdimq,'rqfo')
c      call maxmn(pso,idimt,1,1,' pso')
c     endif
c
c  spherical transform
c
      call fftidim(pso,pss,mkd,-1)
c
      do 310 i = 1, mdim
      qdert(i) = qlnt(i) * wgt(lat)
 310  continue
c
      call symasa(pss,psa,1)
      call fl22a(pss,psa,qo,qdert,1)
c
      call symasa(tts,tta,kdim)
      call fl22a(tts,tta,teo,qdert,kdim)
c
      call symasa(qqs,qqa,kdimq)
      call fl22a(qqs,qqa,rqo,qdert,kdimq)
c
      call gozrmt(qlnt,qlnv,qdert,eps,lat,qlnwct,rcs2,wgt(lat))
c
      call symasa(uus,uua,kdim)
      call symasa(vvs,vva,kdim)
      call psum2a(uua,uus,vva,vvs,zeo,qlnwct,qdert,kdim)
      call msum2a(uua,uus,vva,vvs,dio,qlnwct,qdert,kdim)
c
1000   continue
c
      do 110 k=1,kdim
      do 110 m=1,mdim
      dio(m,k)=-dio(m,k)
  110 continue
c
c ======= check values ========
c
c      print *,' check coef(before transp) after sg2sg'
c      call maxmn(zeo,mdim,1,kdim,' zeo')
c      call maxmn(dio,mdim,1,kdim,' dio')
c      call maxmn(teo,mdim,1,kdim,' teo')
c      call maxmn(rqo,mdim,1,kdimq,' rqo')
c      call maxmn(qo,mdim,1,1,' qo ')
c
c
      return
      end
