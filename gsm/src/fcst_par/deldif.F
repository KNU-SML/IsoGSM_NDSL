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
#define NTRAC ntrac
#define RERTH rerth
      subroutine deldif(rt,w,deltim,qm,sl,x,y,lnt2s,lnoffset)
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    deldif      horizontal diffusion.
c   prgmmr: mark iredell     org: w/nmc23    date: 91-03-15
c
c abstract: horizontal diffusion of temperature, moisture,
c           vorticity and divergence. the implicit linear equation
c           is solved using the lapse rates of globally averaged
c           temperature and moisture to transform the laplacian
c           from constant pressure to constant sigma surfaces.
c           for the t126 operational model, deldif invokes
c           second order leith diffusion only above wavenumber 69
c           with a time scale of 13080 seconds at wavenumber 126.
c
c program history log:
c   91-03-15  mark iredell
c
c usage:    call deldif (rt,w,deltim,qm,sl,x,y)
c   input argument list:
c     rt       - specific humidity
c     w        - vorticity
c     deltim   - timestep
c     qm       - ln(psfc)
c     sl       - sigma layer values
c     x        - divergence
c     y        - temperature
c
c   output argument list:
c     rt       - specific humidity
c     w        - vorticity
c     x        - divergence
c     y        - temperature
c
c remarks: local variables that can be modified to change diffusion are:
c          rtnp- reciprocal of time scale of diffusion at wavenumber np
c          np  - wavenumber at which rtnp diffusion is applied
c          n0  - maximum wavenumber for zero diffusion
c          jdel- order of diffusion
c
c attributes:
c   language: fortran 77.
c   machine:  cray ymp.
c
c$$$
#include <paramodel.h>
#include <comio.h>
#ifdef MP
#include <commpi.h>
#define LNT22S lnt22p_
#else
#define LNT22S lnt22_
#endif
      dimension sl(levs_),qm(LNT22S)
      dimension w(LNT22S,levs_),x(LNT22S,levs_),
     &          y(LNT22S,levs_),rt(LNT22S,levh_)
      dimension dn(lnt2_)
      data ifirst/0/
      save ifirst,dn
c
      parameter(difcof=3.e15,lefres=80,jdel=2)
c
      if(ifirst.eq.0) then
        ifirst=1
c
c  reciprocal of time scale of diffusion at reference wavenumber np
c
        rtnp=difcof/(rerth_**4)*float(lefres*(lefres+1))**2
        np=jcap_
c
c  maximum wavenumber for zero diffusion
c
        n0=0.55*jcap_
c
#ifndef NOPRINT
        if( mype.eq.master ) print 6,rtnp,np,n0,jdel
6       format(' horizontal diffusion parameters'/
     &  '   effective ',6pf10.3,' microhertz at wavenumber ',i4/
     &  '   maximum wavenumber for zero diffusion ',i4/
     &  '   order of diffusion ',i2)
#endif
c
        jdelh=jdel/2
        npd=max(np-n0,0)
        dn1=2.*rtnp/float(npd*(npd+1))**jdelh
        i=0
        do 10 nm=0,jcap_
        do 10 m=0,jcap_-nm
          nd=max(nm+m-n0,0)
          dn(i+1)=dn1*float(nd*(nd+1))**jdelh
          dn(i+2)=dn(i+1)
          i=i+2
10      continue
#ifdef MP
        call spcshfli(dn,lnt2_,1,jcap_,lwvdef)
#endif
#ifdef DBG
        print *,'dn computation completed'
#endif
      endif
c
      istr=1
      if(lnoffset.eq.0) istr=3
#ifdef CRAY_THREAD
cfpp$ nodepchk
#endif
      do k=1,levs_
        do i=istr,lnt2s
          ii=i+lnoffset
          w(i,k)= w(i,k)/(1.+deltim*dn(ii))
          x(i,k)= x(i,k)/(1.+deltim*dn(ii))
          y(i,k)= y(i,k)/(1.+deltim*dn(ii))
        enddo
      enddo
#ifndef NISLQ
      do k=1,levh_
        do i=istr,lnt2s
          ii=i+lnoffset
          rt(i,k)=rt(i,k)/(1.+deltim*dn(ii))
        enddo
      enddo
#endif
#ifdef DBG
        print *,'diffusion computation completed'
#endif
      return
      end
