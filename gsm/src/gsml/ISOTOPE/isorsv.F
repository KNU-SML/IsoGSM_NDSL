#define CAL cal
#define CP cp
#define G g
#define HFUS hfus
#define HVAP hvap
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NSTYPE nstype
#define NVTYPE nvtype
#define RD rd
#define RVRDM1 rvrdm1
#define RV rv
#define SBC sbc
#define T0C t0c
#define LSOIL lsoil
#include <define.h>
      subroutine isorsv(ims2,imx2,dtf,
     &                  rain,raini,
     $                  smsoil,smsoil0,canopy,sheleg,
     $                  drain,
     &                  rresv1,rresv2)
#include <paramodel.h>
c
* temporarily very simple reservoir is considered
c
      dimension rain(imx2),raini(imx2,2)
      dimension resv(imx2,lsoil_),rresv(imx2,2,lsoil_)
      dimension resv0(imx2,lsoil_),rresv0(imx2,2,lsoil_)
      dimension rresv1(imx2,lsoil_),rresv2(imx2,lsoil_)
      dimension smsoil(imx2,lsoil_),canopy(imx2),
     $     sheleg(imx2),zsoil(imx2,lsoil_)
      dimension runof(imx2),drain(imx2)
      dimension smsoil0(imx2,lsoil_)
      integer   k,i,im
      real      rsvlim
      data      rsvlim /50./
      dimension sldpth(4)
      data sldpth /0.1, 0.3, 0.6, 1.0/
      integer   srclv1,srclv2
      real      iflt(imx2,lsoil_)
c
      im = ims2
c
      do i = 1, im
        do k=1,lsoil_
           resv(i,k)=smsoil(i,k)*sldpth(k)*1000.D0
           resv0(i,k)=smsoil0(i,k)*sldpth(k)*1000.D0
        enddo
        resv(i,1)=resv(i,1)+sheleg(i)+canopy(i)
c
        iflt(i,4)=0.
        iflt(i,3)=resv(i,4)-resv0(i,4)+drain(i)*dtf
        iflt(i,2)=resv(i,3)-resv0(i,3)+iflt(i,3)
        iflt(i,1)=resv(i,2)-resv0(i,2)+iflt(i,2)
c
        do k=1,lsoil_
           rresv(i,1,k)=rresv1(i,k)
           rresv(i,2,k)=rresv2(i,k)
           rresv0(i,1,k)=rresv1(i,k)
           rresv0(i,2,k)=rresv2(i,k)
        enddo
c
        do iso=1,2 !! 1 18O, 2 HDO
c          rresv(i,iso)=rresv(i,iso)*rsvlim+raini(i,iso)*1000.D0
c          rresv(i,iso)=rresv(i,iso)/resv(i)
c assuming no fractionation for evaporation and runoff
          do k=lsoil_,2,-1
            if (iflt(i,k-1).ge.0.) then
              srclv1=k-1
            else
              srclv1=k
            endif
            if (iflt(i,k).ge.0.) then
              srclv2=k
            else
              srclv2=k+1
            endif
c            
            rresv(i,iso,k)=rresv0(i,iso,k)
     $           +(iflt(i,k-1)*(rresv0(i,iso,srclv1)-rresv0(i,iso,k))
     $           - iflt(i,k  )*(rresv0(i,iso,srclv2)-rresv0(i,iso,k)))
     $           /resv(i,k)
          enddo
          if (iflt(i,1).ge.0.) then
             srclv2=1
          else
             srclv2=2
          endif
          rresv(i,iso,1)=rresv0(i,iso,1)
     $         +((raini(i,iso)*1000.-rain(i)*1000.D0*rresv0(i,iso,1))
     $         -iflt(i,k)*(rresv0(i,iso,srclv2)-rresv0(i,iso,k)))
     $         /resv(i,1)
c
          do k=1,lsoil_
          if ((rresv(i,iso,k).gt.100).or.rresv(i,iso,k).lt.-100) then
             print*,'isorsv',i,iso,k,
     $            rresv(i,iso,k),resv(i,k),rain(i),raini(i,iso)
          endif
          rresv1(i,k)=rresv(i,1,k)
          rresv2(i,k)=rresv(i,2,k)
          enddo
        enddo
      enddo
c
      return
      end
