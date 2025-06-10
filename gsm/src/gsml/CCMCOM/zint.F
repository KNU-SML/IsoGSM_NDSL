#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine zint(piln    ,pmln    ,rair    ,gravit  ,tv      ,
     $                zm      ,zi      )
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
C
C Compute the geopotential height at the interface points from *HEIGHT
C ABOVE THE SURFACE* at midlayer midpoints, using the supplied virtual
C temperatures and pressures.
C
C---------------------------Code history--------------------------------
C
C Original version:  D. Williamson, J. Hack
C Standardized:      L. Buja, Jun 1992, Feb 1996
C Reviewed:          J. Hack, B. Boville, Aug 1992, Apr 1996
C
C-----------------------------------------------------------------------
c
c $Id: zint.F,v 1.2 2009/10/02 00:35:47 kana Exp $
c $Author: kana $
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Commons----------------------------------
C
C Input arguments
C
      real piln(ILOTS,levs_+1)    ! Log interface pressures
      real pmln(ILOTS,levs_)     ! Log midpoint pressures
      real rair                 ! Gas constant for dry air
      real gravit               ! Acceleration of gravity
      real tv(ILOTS,levs_)       ! Virtual temperature
      real zm(ILOTS,levs_)       ! Height above surface at midpoints
C
C Output arguments
C
      real zi(ILOTS,levs_+1)      ! Height above surface at interfaces
C
C---------------------------Local variables-----------------------------
C
      integer i,k               ! Lon, level indices
      real rog                  ! Rair / gravit
C
C-----------------------------------------------------------------------
C
C Add increment to midlayer height
C
      rog = rair/gravit
      do i=1,ILOTS
        zi(i,1) = zm(i,1) + rog*(pmln(i,1) - piln(i,1))*tv(i,1)
      end do
C
      do k=2,levs_
        do i=1,ILOTS
          zi(i,k) = zm(i,k) + rog*(pmln(i,k) - piln(i,k))*0.5*
     $              (2.0*tv(i,k) - (tv(i,k) - tv(i,k-1))/
     $               (pmln(i,k) - pmln(i,k-1))*
     $               (pmln(i,k) - piln(i,k)))
        end do
      end do
C
C The surface height is zero by definition.
C
      do i=1,ILOTS
        zi(i,levs_+1) = 0.0
      end do
C
      return
      end
 
