         subroutine albedo(vegtyp,rlat,ijmax,mon,alb,frac)
#include <paramodel.h>
c
c***************************************************************
c   this program  modified from yu tai's albedo program 
c
c   Input data:
c     vegtyp: real vegitation index array of ijmax
c     rlat:   real latitude in degree. array of ijmax 
c     mon :   month
c
c   Outputs:
c     alb(*,1) vis spectrum albedo (strong cosz dep.) array of ijmax
c     alb(*,2) vis spectrum albedo (weak cosz dep.) array of ijmax
c     alb(*,3) nir spectrum albedo (strong cosz dep.) array of ijmax
c     alb(*,4) nir spectrum albedo (weak cosz dep.), array of ijmax
c     frac(*,1) fraction of coverage for strong cosz dep. array of ijmax
c     frac(*,2) fraction of coverage for weaj   cosz dep. array of ijmax
c
c     note that frac(*,*)=1.0 or 0.0 for dominant vegetation type
c
c****************************************************************
c
      dimension vegtyp(ijmax),rlat(ijmax)
c
      dimension albv1(13,5),iswx(13),albn1(13,5)
c
      dimension alb(ijmax,4),frac(ijmax,2)
c
c-for original zcos form with error
c
      data  albv1
     1  / .06,.09,.11,.14,.15,.36,.25,.07,.07,.06,.06,.06,.90
     2,   .06,.10,.12,.15,.13,.36,.25,.07,.07,.07,.07,.07,.90
     3,   .06,.11,.13,.16,.12,.36,.25,.07,.07,.08,.08,.08,.90
     4,   .06,.10,.12,.15,.13,.36,.25,.07,.07,.07,.07,.07,.90
     5,   .06,.10,.12,.15,.13,.36,.25,.07,.07,.07,.07,.07,.90 /
c
c-for original zcos form with error
c
      data  albn1
     1  / .06,.28,.32,.33,.30,.51,.40,.24,.24,.28,.31,.27,.65
     2,   .06,.30,.35,.36,.28,.51,.40,.24,.24,.30,.33,.29,.65
     3,   .06,.32,.38,.39,.26,.51,.40,.24,.24,.32,.35,.31,.65
     4,   .06,.30,.35,.36,.28,.51,.40,.24,.24,.30,.33,.29,.65
     5,   .06,.30,.35,.36,.28,.51,.40,.24,.24,.30,.33,.29,.65 /
c
      data  iswx / 3,1,1,1,1,1,1,2,2,2,2,2,2 /
c
c
c --- preset default seasonal index = 1 (winter nh)
c
c
c
       if (mon.eq.12.or.mon.eq.1.or.mon.eq.2) then
          kssn=1
       elseif (mon.eq.3.or.mon.eq.4.or.mon.eq.5) then
          kssn=2
       elseif (mon.eq.6.or.mon.eq.7.or.mon.eq.8) then
          kssn=3
       elseif (mon.eq.9.or.mon.eq.10.or.mon.eq.11) then
          kssn=4
       endif
c
      do ij=1,ijmax
          if(vegtyp(ij).gt.13) then
            print *,'vegatation index is too large'
#ifdef MP
#ifdef RSM
            call rmpabort
#else
            call mpabort
#endif
#else
            call abort
#endif
          endif
c
          idx = vegtyp(ij)+1
          if(kssn.eq.0.or.kssn.gt.4) then
            print *,'illigal season kssn=',kssn
#ifdef MP
#ifdef RSM
            call rmpabort
#else
            call mpabort
#endif
#else
            call abort
#endif
          endif
c
          do ka=1,4
            alb(ij,ka) = 0.0
          enddo
c
c  frac is taken from grib file.  Possible error here
c
          do ka=1,2
            frac(ij,ka)=0.
          enddo
c
          issx=kssn
c
c --if the point is in the sh, reverse  the season
c
          if(rlat(ij).lt.0.) then
            issx=mod(issx+2,4)
            if (issx.eq.0) issx=4
          endif
c
          indx = idx
          isw  = iswx(indx)
          if(isw.le.2) then
            alb(ij,isw)   = albv1(indx,issx) 
            alb(ij,isw+2) = albn1(indx,issx)
            frac(ij,isw)  = 1.
          endif
c
      enddo
c
      return
      end
