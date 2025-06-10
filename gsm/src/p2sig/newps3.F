#define G g
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
#define RD rd
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine newps3(t,q,kdimprs,ps,hold,hnew,psnew,si,sl)
#include <paramodel.h>
      save
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    newps       interpolate surface pressure.
c   prgmmr: mark iredell     org: w/nmc23    date: 91-03-15
c
c abstract: using the hydrostatic equation, the surface pressure
c           is interpolated from input orography to output orography.
c           below the input surface, the temperature lapse rate
c           is fixed at -6.5k/km.  above the input surface,
c           surface pressure is interpolated from height profiles
c           computed as in the mrf model.
c
c program history log:
c   91-03-15  mark iredell
c
c usage:    call newps(t,kdimprs,q0,h0,h1,q1,si)
c   input argument list:
c     t        - temperature
c     q        - specific humidity
c     q0       - old ln(psfc)
c     h0       - old orography
c     h1       - new orography
c     si       - edge sigma values
c
c   output argument list:
c     q1       - new ln(psfc)
c
c   subprograms called:
c     bnmc     - compute mrf temperature to height matrix
c
c attributes:
c   language: fortran 77.
c   machine:  cray ymp.
c
c$$$
c
#include <paramter.h>
#include <parmchgr.h>
c
      parameter(idimt=idim*2,jdimhf=jdim/2)
c
         dimension t(idimt,kdimprs),ps(idimt),psnew(idimt)
         dimension q(idimt,kdimprs)
         dimension hold(idimt,kdimprs),hnew(idimt),si(idimt,kdimprs+1)
         dimension sl(idimt,kdimprs)
         dimension tau(idimt,kdimprs)
         dimension rlsig(idimt,kdimprs+1)
         dimension rllsig(idimt,kdimprs)
         dimension h(idimt,kdimprs+1)
         dimension a(idimt)
         dimension absa(idimt)
c--------
c-------- compute interface heights.
c--------
c        print *,' ------ into newpsfc ---'
         rog = rd_/g_
         do 5 k=1,kdimprs
         do 5 i=1,idimt
         rlsig(i,k) = - alog(si(i,k))
         rllsig(i,k) = - alog(sl(i,k))
         tau(i,k) = (1.+0.61*q(i,k))*t(i,k) * rog
5        continue
         do 10 i=1,idimt
         rlsig(i,kdimprs+1) = - alog ( .005 )
         psnew(i)= - ps(i)
         h(i,1)=hold(i,1)
10       continue
         do 100 k=1,kdimprs
         do 100 i=1,idimt
           h(i,k+1)=hold(i,k)+tau(i,k)*(rlsig(i,k+1)-rllsig(i,k))
100      continue
c--------
c-------- loop over layers, testing to see if new sfc pressure is
c-------- in layer, and obtaining it, if so.
c--------
         eps=1.e-1
         do 200 k=1,kdimprs
c----------
c---------- compute lapse rate
c----------
           kp=min(k+1,kdimprs)
           km=max(1,k-1)
           kppz=min(k+2,kdimprs)
           kpz=kppz-1
           kmmz=max(1,k-2)
           kmz=kmmz+1
           do 2020 iq=1,idimt
              if ((rlsig(iq,kpz)+rlsig(iq,kppz)-
     $           rlsig(iq,kmz)-rlsig(iq,kmmz)).eq.0.) then
                 a(iq)=0.
              else
                 a(iq)=2.0e0*(tau(iq,kp)-tau(iq,km))/
     1          (rlsig(iq,kpz)+rlsig(iq,kppz)-
     $           rlsig(iq,kmz)-rlsig(iq,kmmz))
              endif
           absa(iq)=abs(a(iq))
2020       continue
           do 2030 i=1,idimt
           if(hnew(i).ge.h(i,k).and.hnew(i).le.h(i,k+1).and.
     1        absa(i).gt.eps) psnew(i)=rlsig(i,k)
     1                                -ps(i)+(sqrt(tau(i,k)**2+
     2        2.0e0*a(i)*(hnew(i)-h(i,k))) - tau(i,k) )/a(i)
           if(hnew(i).ge.h(i,k).and.hnew(i).le.h(i,k+1).and.
     1        absa(i).le.eps) psnew(i)=rlsig(i,k)
     2                                -ps(i)+(hnew(i)-h(i,k))/tau(i,k)
2030       continue
200      continue
c--------
c-------- do points which fall below first layer (use fixed lapse rate
c--------  of 6.5 deg per km.
c--------
         gamma=6.5e-3
         gascon=rd_
         g=g_
         c=gascon*gamma/g
         do 8010 iq=1,idimt
         a(iq)=tau(iq,1)*
     1       (1.e0-exp(c*(rlsig(iq,1)-rlsig(iq,2))))/
     $           (rlsig(iq,1)-rlsig(iq,2))
8010     continue
         do 8050 i=1,idimt
         if( hnew(i).le.h(i,1) ) psnew(i)=rlsig(i,1)
     1                                   -ps(i)+(sqrt(tau(i,1)**2+
     1     2.e0*a(i)*(hnew(i)-h(i,1))) - tau(i,1))/a(i)
8050     continue
         icount = 0
         do 8100 i=1,idimt
         psnew(i) = - psnew(i)
ccc      if( psnew(i).eq.ps(i) ) icount = icount + 1
         if( psnew(i).eq.ps(i) ) then
              icount = icount + 1
         else
c         print *,'  gz gzold ps psnew ',hold(i),hnew(i),ps(i),psnew(i)
         endif
8100     continue
c
        write(6,300)icount
300      format(' number of surface pressure points not updated=',i6)
c
         do i=1,idimt
c            do k=1,kdimprs
c         print'(2i5,4f15.7)',i,k,h(i,k),tau(i,k),rlsig(i,k),si(i,k)
c            enddo
         if (abs(ps(i)-psnew(i)).gt.
     $           abs(hold(i,1)-hnew(i))*180.E-6) then
           print*,'PS Difference seems too big...'
         print'(2i5,4f15.7)',i,0,hold(i,1),hnew(i),exp(ps(i))*10,
     $           exp(psnew(i))*10
         print'(2i5,2f15.7)',i,0,(sqrt(tau(i,1)**2+
     1        2.e0*a(i)*(hnew(i)-h(i,1)))- tau(i,1))/a(i),
     $        a(i)
         endif
         enddo
       return
       end
