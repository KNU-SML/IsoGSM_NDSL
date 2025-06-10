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
      subroutine zerdia
c
#include <paramodel.h>
#include <comgda.h>
#ifdef DG3
      dimension gda(nwgda*kdgda)
#endif
#ifdef DG3
      do n=1,nwgda*kdgda
        gda(n)=0.
      enddo
      do 30 jr=1,nrgda
      call putdia(jr,nwgda*kdgda,gda)
#endif
30    continue
c
      return
      end
