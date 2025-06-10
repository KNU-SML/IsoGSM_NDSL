#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine inddia
c  set levels and indices for diagnostics
c
c  use kd as defined below when using subroutine accdia
c
c  kd     cnmgda
c   1     dtlarg
c   2     dtconv
c   3     dqconv
c   4     dtshal
c   5     dqshal
c   6     dtvrdf
c   7     duvrdf
c   8     dvvrdf
c   9     dqvrdf
c   10    dthsw
c   11    dthlw
c   12    cloud
c   13    cvcld
c
#include <paramodel.h>
#include <comgda.h>
c
c  set identification index
c
c  large-scale heating
#ifdef DG3
      cnmgda( 1)='dtlarg  '
      ipugda( 1)=241
      ibmgda( 1)=0
#endif
c  deep convective heating
#ifdef DG3
      cnmgda( 2)='dtconv  '
      ipugda( 2)=242
      ibmgda( 2)=0
#endif
c  deep convective moistening
#ifdef DG3
      cnmgda( 3)='dqconv  '
      ipugda( 3)=243
      ibmgda( 3)=0
#endif
c  shallow convective heating
#ifdef DG3
      cnmgda( 4)='dtshal  '
      ipugda( 4)=244
      ibmgda( 4)=0
#endif
c  shallow convective moistening
#ifdef DG3
      cnmgda( 5)='dqshal  '
      ipugda( 5)=245
      ibmgda( 5)=0
#endif
c  vertical diffusion of temperature
#ifdef DG3
      cnmgda( 6)='dtvrdf  '
      ipugda( 6)=246
      ibmgda( 6)=0
#endif
c  vertical diffusion of zonal wind
#ifdef DG3
      cnmgda( 7)='duvrdf  '
      ipugda( 7)=247
      ibmgda( 7)=0
#endif
c  vertical diffusion of meridional wind
#ifdef DG3
      cnmgda( 8)='dvvrdf  '
      ipugda( 8)=248
      ibmgda( 8)=0
#endif
c  vertical diffusion of moisture
#ifdef DG3
      cnmgda( 9)='dqvrdf  '
      ipugda( 9)=249
      ibmgda( 9)=0
#endif
c  short wave radiation heating
#ifdef DG3
      cnmgda(10)='dthsw   '
      ipugda(10)=250
      ibmgda(10)=0
#endif
c  long wave radiation heating
#ifdef DG3
      cnmgda(11)='dthlw   '
      ipugda(11)=251
      ibmgda(11)=0
#endif
c  cloud amount
#ifdef DG3
      cnmgda(12)='cloud   '
      ipugda(12)=71
      ibmgda(12)=0
#endif
c  convective cloud amount
#ifdef DG3
      cnmgda(13)='cvcld   '
      ipugda(13)=72
      ibmgda(13)=0
#endif
c
      return
      end
