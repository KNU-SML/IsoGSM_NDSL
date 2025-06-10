      subroutine cldpar
!----------------------------------------------------------------------
!
!
#include <paramodel.h>
!
#undef DEBUG
!
      real n0r,n0s,n0smax,n0g
      real lamdarmax,lamdasmax,lamdagmax
      common /concld35/den0,denr,cpd,cpv,xlv0,xlv1,xlf0,xls,
     1             sigtop,pidn0r,pidn0s,alpha,n0r,n0s,
     2             bvtr,bvtr2,bvtr3,pvtr,pacrr,precr1,precr2,
     3  betai,xn0,xm0,xmmax,bvts,bvts2,bvts3,pvts,pacrs,precs1,precs2,
     4  t40c,pacrc,pfrz1,pfrz2,qmin,qcrmin,
     5  n0smax,dimax,roqimax,
     6       rslopermax,rslopesmax,rslopegmax,       
     7       rsloperbmax,rslopesbmax,rslopegbmax,   
     8       rsloper2max,rslopes2max,rslopeg2max,  
     9       rsloper3max,rslopes3max,rslopeg3max
      common /concld6/pidn0g,bvtg,precg1,precg2,pvtg,qs0,n0g,
     1                dens,pacrg,avtr,g6pbr,cl

!
!.... constants which may not be tunable
      pi = 4.*atan(1.)
      t0c = 273.15
      den0 = 1.28
      denr = 1000.
      dens = 100.
      cpd = 1005.7
      cpv = 1870.0
      cl  = 4190.0
      xlv0 = 2.501e6   ! latent heat of vaporization at 0c
      xlv1 = cl-cpv
      xlf0 = .3337e6
      xls = xlv0 + xlf0
!
!.... sigma values for possible cloud top
      sigtop = .08
!
!.... tunable parameters for warm rain process
      n0r = 8.e6
#ifdef DEN_SNOW_OPT
      alpha = .1  
#else
      alpha = .12        ! .122 exponen factor for n0s ! for cloud3.2 only
#endif
!
      avtr = 841.9
      bvtr = 0.8
      avts = 11.72
      bvts = .41
      avtg = 330.
      bvtg = .8
      deng = 500.
      n0g = 4.e6
!
      bvtr1 = 1.+bvtr
      bvtr2 = 2.5+.5*bvtr
      bvtr3 = 3.+bvtr
      bvtr4 = 4.+bvtr
      bvtr6 = 6.+bvtr
      g1pbr = rgmma(bvtr1)
      g3pbr = rgmma(bvtr3)
      g4pbr = rgmma(bvtr4)            ! 17.837825
      g6pbr = rgmma(bvtr6)
      g5pbro2 = rgmma(bvtr2)          ! 1.8273
      pvtr = avtr*g4pbr/6.
      eacrr = 1.0
      pacrr = pi*n0r*avtr*g3pbr*.25*eacrr
      precr1 = 2.*pi*n0r*.78
      precr2 = 2.*pi*n0r*.31*avtr**.5*g5pbro2
!
!.... tunable parameters for cold rain process
      betai = .6
      xn0 = 1.e-2
#ifdef NEW_ICE_GEN
      dicon = 11.9
#else
      dicon = 16.3
#endif
      di0 = 12.9e-6
      dimax = 500.e-6
#ifdef N0S_CONST
      n0s = 2.e7
#else
#ifdef DEN_SNOW_OPT
      n0s = 1.e7             ! temperature dependent n0s
#else
      n0s = 2.e6             ! temperature dependent n0s
#endif
#endif
!
      xm0  = (di0/dicon)**2
      xmmax = (dimax/dicon)**2
      bvts1 = 1.+bvts
      bvts2 = 2.5+.5*bvts
      bvts3 = 3.+bvts
      bvts4 = 4.+bvts
      g1pbs = rgmma(bvts1)    !.8875
      g3pbs = rgmma(bvts3)
      g4pbs = rgmma(bvts4)    ! 12.0786
      g5pbso2 = rgmma(bvts2)
      pvts = avts*g4pbs/6.
      pacrs = pi*n0s*avts*g3pbs*.25
      precs1 = 4.*n0s*.65
      precs2 = 4.*n0s*.44*avts**.5*g5pbso2
!
      pidn0r =  pi*denr*n0r
#ifdef DEN_SNOW_WATER
      pidn0s =  pi*denr*n0s
#else
      pidn0s =  pi*dens*n0s
#endif
!
      lamdarmax = 8.e4
      lamdasmax = 1.e5
      lamdagmax = 6.e4
!
      rslopermax = 1./lamdarmax
      rslopesmax = 1./lamdasmax
      rslopegmax = 1./lamdagmax
      rsloperbmax = rslopermax ** bvtr
      rslopesbmax = rslopesmax ** bvts
      rslopegbmax = rslopegmax ** bvtg
      rsloper2max = rslopermax * rslopermax
      rslopes2max = rslopesmax * rslopesmax
      rslopeg2max = rslopegmax * rslopegmax
      rsloper3max = rsloper2max * rslopermax
      rslopes3max = rslopes2max * rslopesmax
      rslopeg3max = rslopeg2max * rslopegmax
!
!     additional parameters for mixed phase clouds
!
      t40c = t0c - 40.
      eacrc = 1.0
      pacrc = pi*n0s*avts*g3pbs*.25*eacrc
      pfrz1 = 100.
      pfrz2 = 0.66
!
!     additional parameters
!
      qmin=1.e-15
      qcrmin=1.e-9
      roqimax = 2.08e22*dimax**8
#ifndef N0S_CONST
      n0smax=1.e11
#endif
!
      qs0 = 6.e-4
      bvtg1 = 1.+bvtg
      bvtg2 = 2.5+.5*bvtg
      bvtg3 = 3.+bvtg
      bvtg4 = 4.+bvtg
      g1pbg = rgmma(bvtg1)
      g3pbg = rgmma(bvtg3)
      g4pbg = rgmma(bvtg4)
      pacrg = pi*n0g*avtg*g3pbg*.25
      g5pbgo2 = rgmma(bvtg2)
      pvtg = avtg*g4pbg/6.
      precg1 = 2.*pi*n0g*.78
      precg2 = 2.*pi*n0g*.31*avtg**.5*g5pbgo2
      pidn0g =  pi*deng*n0g
!
      return
      end
! ...................................................................
      real function rgmma(x)
!     rgmma function:  use infinite product form
      parameter (euler=0.577215664901532)
      if(x.eq.1.)then
        rgmma=0.
          else
        rgmma=x*exp(euler*x)
        do i=1,10000
          y=float(i)
          rgmma=rgmma*(1.000+x/y)*exp(-x/y)
        enddo
        rgmma=1./rgmma
      endif
      return
      end
