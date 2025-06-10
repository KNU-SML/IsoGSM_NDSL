#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine buoyan(q       ,t       ,p       ,z       ,pf      ,
     $                  tp      ,qstp    ,tl      ,rl      ,cape    ,
     $                  pblt    ,lcl     ,lel     ,lon     ,mx      ,
     $                  rd      ,grav    ,cp      ,msg     ,nstep   ,
     $                  lat     ,tpert   ,qpert   )
#include <paramodel.h>
#ifdef MP
#ifdef RSM
#ifdef RSMVECTORIZE
#define ILOTS igrd12p_*jgrd12p_
#else
#define ILOTS igrd12p_
#endif
#else
#define ILOTS lonf2p_
#endif
#else
#undef RSMVECTORIZE
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
C-----------------------------------------------------------------------
C This is contributed code not fully standardized by the CCM core group.
C
c the documentation has been enhanced to the degree that we are able
c
C Original version:  G. Zhang and collaborators
C Standardized:      Core group staff, 1994 and 195
C Reviewed:          P. Rasch, April 1996
c-----------------------------------------------------------------------
c  $Id: buoyan.F,v 1.2 2009/10/02 00:35:46 kana Exp $
c-----------------------------------------------------------------------
#include <implicit.h>
C-----------------------------------------------------------------------
#include <guang.h>
C-----------------------------------------------------------------------
C
C jul 14/92 - guang jun zhang, m.lazare, n.mcfarlane.  as in
C             previous version buoyan except remove pathalogical
C             cases of "zig-zags" in profiles where lel defined
C             far too high (by use of lelten array, which assumes
C             a maximum of five such crossing points).
C feb 18/91 - guang jun zhang, m.lazare, n.mcfarlane.  previous
C             version buoyan.
C
C input arguments
C
      real q(ILOTS,levs_)        ! spec. humidity
      real t(ILOTS,levs_)        ! temperature
      real p(ILOTS,levs_)        ! pressure
      real z(ILOTS,levs_)        ! height
      real pf(ILOTS,levs_+1)     ! pressure at interfaces
      real pblt(ILOTS)          ! index of pbl depth
      real tpert(ILOTS)         ! perturbation temperature by pbl processes
      real qpert(ILOTS)         ! perturbation moisture by pbl processes
c
c output arguments
c
      real tp(ILOTS,levs_)       ! parcel temperature
      real qstp(ILOTS,levs_)     ! saturation mixing ratio of parcel
      real tl(ILOTS)            ! parcel temperature at lcl
      real cape(ILOTS)          ! convective aval. pot. energy.
      integer lcl(ILOTS)        ! 
      integer lel(ILOTS)        ! 
      integer lon(ILOTS)        ! level of onset of deep convection
      integer mx(ILOTS)         ! level of max moist static energy
C
C--------------------------Local Variables------------------------------
C
      real capeten(ILOTS,5)     ! provisional value of cape
      real tv(ILOTS,levs_)       ! 
      real tpv(ILOTS,levs_)      ! 
      real buoy(ILOTS,levs_)

      real a1(ILOTS) 
      real a2(ILOTS) 
      real estp(ILOTS) 
      real pl(ILOTS) 
      real plexp(ILOTS) 
      real hmax(ILOTS) 
      real hmn(ILOTS) 
      real y(ILOTS)

      logical plge600(ILOTS) 
      integer knt(ILOTS) 
      integer lelten(ILOTS,5)

      real cp
      real e
      real grav

      integer i
      integer k
      integer lat
      integer msg
      integer n
      integer nstep

      real rd
      real rl
C
C-----------------------------------------------------------------------
C
      do n = 1,5
        do i = 1,ILOTS
          lelten(i,n) = levs_
          capeten(i,n) = 0.
        end do
      end do
C
      do i = 1,ILOTS
        lon(i) = levs_
        knt(i) = 0
        lel(i) = levs_
        mx(i) = lon(i)
        cape(i) = 0.
        hmax(i) = 0.
      end do
C
C set "launching" level(mx) to be at maximum moist static energy.
C search for this level stops at planetary boundary layer top.
C
      do k = levs_,msg + 1,-1
        do i = 1,ILOTS
          hmn(i) = cp*t(i,k) + grav*z(i,k) + rl*q(i,k)
          if (k.ge.nint(pblt(i)) .and. k.le.lon(i) .and.
     $         hmn(i).gt.hmax(i)) then
            hmax(i) = hmn(i)
            mx(i) = k
          end if
        end do
      end do
C
      do i = 1,ILOTS
        lcl(i) = mx(i)
        e = p(i,mx(i))*q(i,mx(i))/ (eps1+q(i,mx(i)))
        tl(i) = 2840./ (3.5*log(t(i,mx(i)))-log(e)-4.805) + 55.
        if (tl(i).lt.t(i,mx(i))) then
          plexp(i) = (1./ (0.2854* (1.-0.28*q(i,mx(i)))))
          pl(i) = p(i,mx(i))* (tl(i)/t(i,mx(i)))**plexp(i)

        else
          tl(i) = t(i,mx(i))
          pl(i) = p(i,mx(i))
        end if
      end do
C
C calculate lifting condensation level (lcl).
C
      do k = levs_,msg + 2,-1
        do i = 1,ILOTS
          if (k.le.mx(i) .and. (p(i,k).gt.pl(i).and.
     $         p(i,k-1).le.pl(i))) then
            lcl(i) = k - 1
          end if
        end do
      end do
C
C if lcl is above the nominal level of non-divergence (600 mbs),
C no deep convection is permitted (ensuing calculations
C skipped and cape retains initialized value of zero).
C
      do i = 1,ILOTS
        plge600(i) = pl(i).ge.600.
      end do
C
C initialize parcel properties in sub-cloud layer below lcl.
C
      do k = levs_,msg + 1,-1
        do i=1,ILOTS
          if (k.gt.lcl(i) .and. k.le.mx(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
            qstp(i,k) = q(i,mx(i))
            tp(i,k) = t(i,mx(i))* (p(i,k)/p(i,mx(i)))**
     $                (0.2854* (1.-0.28*q(i,mx(i))))
C
C buoyancy is increased by 0.5 k as in tiedtke
C
c-jjh          tpv (i,k)=tp(i,k)*(1.+1.608*q(i,mx(i)))/
c-jjh     1                     (1.+q(i,mx(i)))
            tpv(i,k) = (tp(i,k)+tpert(i))*
     $                 (1.+1.608*q(i,mx(i)))/ (1.+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5
          end if
        end do
      end do
C
C define parcel properties at lcl (i.e. level immediately above pl).
C
      do k = levs_,msg + 1,-1
        do i=1,ILOTS
          if (k.eq.lcl(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
            qstp(i,k) = q(i,mx(i))
            tp(i,k) = tl(i)* (p(i,k)/pl(i))**
     $                (0.2854* (1.-0.28*qstp(i,k)))
C              estp(i)  =exp(a-b/tp(i,k))
C use of different formulas for est has about 1 g/kg difference
C in qs at t= 300k, and 0.02 g/kg at t=263k, with the formula
C above giving larger qs.
C
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/
     $                ((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
            a1(i) = cp/rl + qstp(i,k)* (1.+qstp(i,k)/eps1)*rl*
     $              eps1/ (rd*tp(i,k)**2)
            a2(i) = .5* (qstp(i,k)* (1.+2./eps1*qstp(i,k))*
     $              (1.+qstp(i,k)/eps1)*eps1**2*rl*rl/
     $              (rd**2*tp(i,k)**4)-qstp(i,k)*
     $              (1.+qstp(i,k)/eps1)*2.*eps1*rl/
     $              (rd*tp(i,k)**3))
            a1(i) = 1./a1(i)
            a2(i) = -a2(i)*a1(i)**3
            y(i) = q(i,mx(i)) - qstp(i,k)
            tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
c          estp(i)  =exp(a-b/tp(i,k))
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/
     $              ((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
C
C buoyancy is increased by 0.5 k in cape calculation.
C dec. 9, 1994
c-jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/(1.+q(i,mx(i)))
C
            tpv(i,k) = (tp(i,k)+tpert(i))* (1.+1.608*qstp(i,k))/
     $                 (1.+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5
          end if
        end do
      end do
C
C main buoyancy calculation.
C
      do k = levs_ - 1,msg + 1,-1
        do i=1,ILOTS
          if (k.lt.lcl(i) .and. plge600(i)) then
            tv(i,k) = t(i,k)* (1.+1.608*q(i,k))/ (1.+q(i,k))
            qstp(i,k) = qstp(i,k+1)
            tp(i,k) = tp(i,k+1)* (p(i,k)/p(i,k+1))**
     $                (0.2854* (1.-0.28*qstp(i,k)))
c          estp(i) = exp(a-b/tp(i,k))
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/
     $                ((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
            a1(i) = cp/rl + qstp(i,k)* (1.+qstp(i,k)/eps1)*rl*
     $              eps1/ (rd*tp(i,k)**2)
            a2(i) = .5* (qstp(i,k)* (1.+2./eps1*qstp(i,k))*
     $              (1.+qstp(i,k)/eps1)*eps1**2*rl*rl/
     $              (rd**2*tp(i,k)**4)-qstp(i,k)*
     $              (1.+qstp(i,k)/eps1)*2.*eps1*rl/
     $              (rd*tp(i,k)**3))
            a1(i) = 1./a1(i)
            a2(i) = -a2(i)*a1(i)**3
            y(i) = qstp(i,k+1) - qstp(i,k)
            tp(i,k) = tp(i,k) + a1(i)*y(i) + a2(i)*y(i)**2
c          estp(i)  =exp(a-b/tp(i,k))
            estp(i) = c1*exp((c2* (tp(i,k)-tfreez))/
     $                ((tp(i,k)-tfreez)+c3))

            qstp(i,k) = eps1*estp(i)/ (p(i,k)-estp(i))
c-jjh          tpv(i,k) =tp(i,k)*(1.+1.608*qstp(i,k))/
cjt            (1.+q(i,mx(i)))
            tpv(i,k) = (tp(i,k)+tpert(i))* (1.+1.608*qstp(i,k))/
     $                 (1.+q(i,mx(i)))
            buoy(i,k) = tpv(i,k) - tv(i,k) + 0.5
          end if
        end do
      end do
C
      do k = msg + 2,levs_
        do i = 1,ILOTS
          if (k.lt.lcl(i) .and. plge600(i)) then
            if (buoy(i,k+1).gt.0. .and. buoy(i,k).le.0.) then
              knt(i) = min(5,knt(i) + 1)
              lelten(i,knt(i)) = k
            end if
          end if
        end do
      end do
C
C calculate convective available potential energy (cape).
C
      do n = 1,5
        do k = msg + 1,levs_
          do i = 1,ILOTS
            if (plge600(i) .and.k.le.mx(i) .and.k.gt.lelten(i,n)) then
              capeten(i,n) = capeten(i,n) +
     $                       rd*buoy(i,k)*log(pf(i,k+1)/pf(i,k))
            end if
          end do
        end do
      end do
C
C find maximum cape from all possible tentative capes from
C one sounding,
C and use it as the final cape, april 26, 1995
C
      do n = 1,5
        do i = 1,ILOTS
          if (capeten(i,n).gt.cape(i)) then
            cape(i) = capeten(i,n)
            lel(i) = lelten(i,n)
          end if
        end do
      end do
C
C put lower bound on cape for diagnostic purposes.
C
      do i = 1,ILOTS
        cape(i) = max(cape(i), 0.)
      end do
C
      return
      end
 
