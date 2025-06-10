#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine virtem(t       ,q       ,zvir    ,tv      )
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
C Compute the virtual temperature.
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Boville
C Standardized:      J. Rosinski, June 1992
C Reviewed:          D. Williamson, J. Hack, August 1992
C Reviewed:          D. Williamson, March 1996
C
C-----------------------------------------------------------------------
c
c $Id: virtem.F,v 1.2 2009/10/02 00:35:47 kana Exp $
c
#include <implicit.h>
C------------------------------Parameters-------------------------------
C
C Input arguments
C
      real t(ILOTS,levs_)       ! temperature
      real q(ILOTS,levs_)       ! specific humidity
      real zvir                ! virtual temperature constant
C
C Output arguments
C
      real tv(ILOTS,levs_)      ! virtual temperature
C
C---------------------------Local storage-------------------------------
C
      integer i,k              ! longitude and level indexes
C
      do k=1,levs_
        do i=1,ILOTS
          tv(i,k) = t(i,k)*(1.0 + zvir*q(i,k))
        end do
      end do
C
      return
      end
 
