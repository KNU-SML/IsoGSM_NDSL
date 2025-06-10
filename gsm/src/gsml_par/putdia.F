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
      subroutine putdia(j,l,a)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c  put diagnostics onto disk or memory
c  call another getdia or putdia to finish i/o
c
c  j is record number
c  l is record length
c  a is array of length l to get
c
#include <paramodel.h>
#ifdef DG3
      dimension a(l)
#include <comgda.h>
c
      do n=j*l-l+1,j*l
        gdd(n)=a(n-(j*l-l+1)+1)
      enddo
#endif
c
      return
c
      end
