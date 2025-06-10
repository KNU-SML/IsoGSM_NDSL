#include <define.h>
#ifndef RSM
#ifdef MP
#define LONF2S lonf2p_
#else
#define LONF2S lonf2_
#endif
#else
#ifdef MP
#ifdef RSMVECTORIZE
#define LONF2S igrd12p_*jgrd12p_
#else
#define LONF2S igrd12p_
#endif
#else
#define LONF2S igrd12_
#endif
#endif

#define IGRD1 igrd1
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD1 jgrd1
#define JGRD1P jgrd1p
#define JGRD jgrd
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define LONF lonf
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONFP lonfp
#define LEVS levs

      module ppgrid
!-------------------------------------------------------
!
!-------------------------------------------------------
#include <paramodel.h>

!        integer,parameter :: pcols = PCOLS   ! number of columns (max)
        integer,parameter :: pcols = LONF2S    ! number of columns (max)
        integer,parameter :: pver  = levs_     ! number of vertical levels
        integer,parameter :: pverp = pver + 1  ! pver + 1

      end module ppgrid
