#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JCAPP jcapp
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LCAPP lcapp
#define LEVH levh
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LNUV lnuv
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NTRAC ntrac
#define TWOJ1 twoj1
#include <define.h>
      subroutine wriznl(nznl)
#include <paramodel.h>
#include <comfibm.h>
      real dswsfc(lonf2_,latg2_),uswsfc(lonf2_,latg2_)
      real dumtn(lonf2_,latg2_)
c.................................................................
#ifdef DGZ
      secphy=shour
      secswr=max(shour,3600.*dtswav)
      seclwr=max(shour,3600.*dtlwav)
      do j=1,latg2_
        do i=1,lonf2_
          dswsfc(i,j)=fluxr(i,j,4)
          uswsfc(i,j)=fluxr(i,j,3)
        enddo
      enddo
c
      if(nznl.gt.0) then
        call mtntrq(shour,snnp1,colrad,gz,z00,psmean,dumtn)
        call znlavz(latg2_,lonf2_,secphy,secswr,
     &              dumtn,slmsk,sheleg,dswsfc,uswsfc)
        call znldia(nznl,thour,idate,kdt,latg2_,levs_,wgb,colrab,del)
      endif
#endif
c
      return
      end
