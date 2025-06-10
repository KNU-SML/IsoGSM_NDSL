#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine zmid(pstarln ,pmln    ,rair    ,gravit  ,tv      ,
     $                z       )
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
C Compute the geopotential height *ABOVE THE SURFACE* at layer
C midpoints from the virtual temperatures and pressures.
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Boville, D. Williamson, Jan 1990
C Standardized:      L. Buja, Jun 1992, Feb 1996
C Reviewed:          J. Hack, B. Boville, Aug 1992, Apr 1996
C
C-----------------------------------------------------------------------
c
c $Id: zmid.F,v 1.2 2009/10/02 00:35:47 kana Exp $
c
C-----------------------------------------------------------------------
#include <implicit.h>
C------------------------------Parameters-------------------------------
C
C Input arguments
C
      real pstarln(ILOTS)        ! Log surface pressures
      real pmln(ILOTS,levs_)      ! Log midpoint pressures
      real rair                  ! Gas constant for dry air
      real gravit                ! Acceleration of gravity
      real tv(ILOTS,levs_)        ! Virtual temperature
C
C Output arguments
C
      real z(ILOTS,levs_)         ! Height above surface at midpoints
C
C---------------------------Local variables-----------------------------
C
      integer i,k,l              ! Lon, level, level indices
      real rog                   ! Rair / gravit
C
C-----------------------------------------------------------------------
C
C Diagonal term of hydrostatic equation
C
      rog = rair/gravit
      do k=1,levs_-1
        do i=1,ILOTS
          z(i,k) = rog*tv(i,k)*0.5*(pmln(i,k+1) - pmln(i,k))
        end do
      end do
      do i=1,ILOTS
        z(i,levs_) = rog*tv(i,levs_)*(pstarln(i) - pmln(i,levs_))
      end do
C
C Bottom level term of hydrostatic equation
C
      do  k=1,levs_-1
        do i=1,ILOTS
          z(i,k) = z(i,k) + rog*tv(i,levs_)*(pstarln(i) -
     $                      0.5*(pmln(i,levs_-1) + pmln(i,levs_)))
        end do
      end do
C
C Interior terms of hydrostatic equation
C
      do k=1,levs_-2
        do l=k+1,levs_-1
          do i=1,ILOTS
            z(i,k) = z(i,k) + rog*(tv(i,l)) *
     $                        0.5*(pmln(i,l+1) - pmln(i,l-1))
          end do
        end do
      end do
C
      return
      end
 
