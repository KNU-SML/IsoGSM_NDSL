      subroutine greduceg(qtt,lnt2,jcap,qttcut,lcapd,lonfd)
c
c subroutine greduceg		programmer: hann-ming juang
c
c purpose: global reduce grid initial routine to compute lcapd and lonfd
c          the maximal value of qtt and the accuray of the digit (ndigit)
c          are used to determine the resolution lcapd and grid point
c          lonfd for each latitude. The lonfd is determined by the
c          checks of factors only 2, 3 and 5 with at least one 2, and
c          it is larger than lcapd*3+1.
c
c input:
c	qtt	coefficient for spectral transform at given lat
c	lnt2	dimension of qtt at given lat
c	jcap	wave dimension
c	qttcut	minimax value of accuracy for qtt
c output:
c	lcapd	wave resolution for reduced grid
c	lonfd	number of reduced grid point for latitude
c
#include <paramodel.h>
      dimension qtt(lnt2)
c
      ind=0
      mwave=0
      do m=0,jcap
        need=0
        do n=m,jcap
          ind=ind+1
          ireal=2*ind-1
          imagi=2*ind
          if(abs(qtt(ireal)).ge.qttcut.and.
     1       abs(qtt(imagi)).ge.qttcut) then
             need=1
          endif
        enddo
        mwave=mwave+need
      enddo
      lcapd=mwave
c
      lonfi=3*(lcapd-1)+1
      lonff=lonfi+mod(lonfi,2)
      do 100 ii=0,100,2
        lonf=lonff+ii
        lonfo=lonf
        jtime=nint(log(float(lonf))/log(2.))
#ifdef DCRFT
        jtime=min(jtime,25)
#endif
        do j=1,jtime
          if( mod(lonf,2).eq.0 ) then
            lonf=lonf/2
            if( lonf.eq.1 ) go to 200
          endif
        enddo
        ktime=nint(log(float(lonf))/log(3.))
#ifdef DCRFT
        ktime=min(ktime,2)
#endif
        do k=1,ktime
          if( mod(lonf,3).eq.0 ) then
            lonf=lonf/3
            if( lonf.eq.1 ) go to 200
          endif
        enddo
        ltime=nint(log(float(lonf))/log(5.))
#ifdef DCRFT
        ltime=min(ltime,1)
#endif
        do l=1,ltime
          if( mod(lonf,5).eq.0 ) then
            lonf=lonf/5
            if( lonf.eq.1 ) go to 200
          endif
        enddo
#ifdef DCRFT
        if( mod(lonf,7).eq.0 ) then
          lonf=lonf/7
          if( lonf.eq.1 ) go to 200
        endif
        if( mod(lonf,11).eq.0 ) then
          lonf=lonf/11
          if( lonf.eq.1 ) go to 200
        endif
#endif
 100  continue
 200  lonfd=lonfo
c
      return
      end   
