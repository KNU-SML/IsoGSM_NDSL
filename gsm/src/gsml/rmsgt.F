#define JCAP jcap
#define LATG latg
#define LEVH levh
#define LEVS levs
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define NCLDG ncldg
#define NPES npes
#define NTOTAL ntotal
#define NTRAC ntrac
      subroutine rmsgt(q,x,y,w,del,r)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    rmsgt       computes root mean square.
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-25
c
c abstract: computes the root mean square in each level and of
c   the vertical integral (when appropriate) given the spectral
c   coefficients of the model variables or the tendencies.
c   the results are printed.  rmsgt is strictly diagnostic.
c
c program history log:
c   88-04-25  joseph sela
c   02-09-05  yifeng cui, if rsm, undef mp
c
c usage:    call rmsgt (q, x, y, w, del, r)
c   input argument list:
c     q        - spectral coefs of ln(psfc)     or its tendency.
c     x        - spectral coefs of divergence   or its tendency.
c     y        - spectral coefs of temperature  or its tendency.
c     w        - spectral coefs of vorticity    or its tendency.
c     del      - sigma spacing at each layer.
c     r        - spectral coefs of mixing ratio or its tendency.
c
c   output files:
c     output   - print file.
c
c attributes:
c   language: fortran 200.
c   machine:  cyber 205.
c
c$$$
#include <paramodel.h>
#ifdef RSM
#undef MP
#endif
#ifdef MP
#include <commpi.h>
      real,allocatable::qf(:),xf(:,:),yf(:,:),wf(:,:),rf(:,:)
#endif
#ifdef MP
#define LNT22S lnt22p_
#else
#define LNT22S lnt22_
#endif
      dimension q(LNT22S),x(LNT22S,levs_),y(LNT22S,levs_)
      dimension w(LNT22S,levs_),r(LNT22S,levh_)
      dimension rx(levs_),ry(levs_),rw(levs_),del(levs_),rr(levh_)
#ifdef MP
      allocate(qf(lnt22_))
      allocate(xf(lnt22_,levs_))
      allocate(yf(lnt22_,levs_))
      allocate(wf(lnt22_,levs_))
      allocate(rf(lnt22_,levh_))
c
      call mpsp2f(q,lnt22p_,qf,lnt22_,1)
      call mpsp2f(x,lnt22p_,xf,lnt22_,levs_)
      call mpsp2f(y,lnt22p_,yf,lnt22_,levs_)
      call mpsp2f(w,lnt22p_,wf,lnt22_,levs_)
      call mpsp2f(r,lnt22p_,rf,lnt22_,levh_)
      if( mype.eq.master ) then
#ifdef DBG
       print *,' In rmsgt PE0'
       call maxmin(qf,lnt22_,1,1,1,'q')
       call maxmin(yf,lnt22_,levs_,1,levs_,'y')
       call maxmin(xf,lnt22_,levs_,1,levs_,'x')
       call maxmin(wf,lnt22_,levs_,1,levs_,'w')
       call maxmin(rf,lnt22_,levh_,1,levh_,'r')
#endif
#define QS qf
#define XS xf
#define YS yf
#define WS wf
#define RS rf
#else
#define QS q
#define XS x
#define YS y
#define WS w
#define RS r
#endif
      vr=0.e0
      vx=0.e0
      vy=0.e0
      vw=0.e0
      do 10 k=1,levs_
      call bartri(XS(1,k),XS(1,k),rx(k))
      call bartri(YS(1,k),YS(1,k),ry(k))
      call bartri(WS(1,k),WS(1,k),rw(k))
      vx=vx+rx(k)*del(k)
      vy=vy+ry(k)*del(k)
      vw=vw+rw(k)*del(k)
10    continue
      do 20 n=1,ntotal_
      do 20 k=1,levs_
      nk = (n-1)*levs_+k
      call bartri(RS(1,nk),RS(1,nk),rr(nk))
      vr=vr+rr(nk)*del(k)
20    continue
      call bartri(QS,QS,rq)
c
      print 100,vx,vw,vy,vr,rq
      do k=1,levs_
        print 200,rx(k),rw(k),ry(k),rr(k)
      enddo
      if(ntotal_.ge.2) then
        print 300,ntotal_ - 1
        do k = 1,levs_
          ki = levs_+k
          ke = (ntotal_-1)*levs_  + k
            if(ntotal_.eq.2) print 402, (rr(kk),kk=ki,ke,levs_)
            if(ntotal_.eq.3) print 403, (rr(kk),kk=ki,ke,levs_)
            if(ntotal_.eq.4) print 404, (rr(kk),kk=ki,ke,levs_)
            if(ntotal_.eq.5) print 405, (rr(kk),kk=ki,ke,levs_)
            if(ntotal_.eq.6) print 406, (rr(kk),kk=ki,ke,levs_)
        enddo
      endif
100   format(1h0,'div vort temp mixratio ln(ps) ',5(e15.8,1x))
200   format(1h ,4(2x,e17.10))
300   format(1h0,'tracers and clouds = ',i5)
402   format(1x,1(2x,e17.10))
403   format(1x,2(2x,e17.10))
404   format(1x,3(2x,e17.10))
405   format(1x,4(2x,e17.10))
406   format(1x,5(2x,e17.10))
#undef QS
#undef XS
#undef YS
#undef WS
#undef RS
#ifdef MP
      endif
      deallocate(qf)
      deallocate(xf)
      deallocate(yf)
      deallocate(wf)
      deallocate(rf)
#endif
      return
      end
