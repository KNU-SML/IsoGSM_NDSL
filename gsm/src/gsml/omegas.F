#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
      subroutine omegas(ims2,imx2,kmx,
     1      dphi,dlam,cg,ug,vg,dg,del,rcl,vvel,ps,sl)
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c....   code lifted from post (mcp1840) jun 88--computes vvel (cb/sec)
c....    input ps in cb,output vvel in cb/sec
c....   do loops altered for better vectorization possibilities..k.a.c.
#include <paramodel.h>
      dimension dphi(imx2),dlam(imx2),
     1 cg(imx2,kmx),ug(imx2,kmx),vg(imx2,kmx),
     2 dg(imx2,kmx),del(imx2,kmx),ps(imx2),sl(imx2,kmx)
c...   vvel contains omega in layers on return from subroutine...
      dimension vvel(imx2,kmx)
c local array
c
c prgmmr:  yifeng cui added rsm mp    org:sdsc     date:02-09-05
c
#ifdef MP
#ifdef RSM
#define ILOTS igrd12p_
#else
#define ILOTS lonf2p_
#endif
#else
#ifdef RSM
#define ILOTS igrd12_
#else
#define ILOTS lonf2_
#endif
#endif
      dimension db(ILOTS,levs_),
     1          cb(ILOTS,levs_),
     2          dot(ILOTS,levs_+1)
      nx =ims2
      nz =kmx
      do 1 k=1,nz+1
        do 1 lo=1,nx
          dot(lo,k) = 0.e0
    1 continue
c...  compute c=v(true)*del(ln(ps)).divide by cos for del cos for v
      do 3 lo=1,nx
        dphi(lo)=dphi(lo)*rcl
        dlam(lo)=dlam(lo)*rcl
    3 continue
      do 5 le=1,nz
        do 4 lo=1,nx
          cg(lo,le)=ug(lo,le)*dlam(lo)+vg(lo,le)*dphi(lo)
    4   continue
    5 continue
      do 10 lo=1,nx
        db(lo,1)=del(lo,1)*dg(lo,1)
        cb(lo,1)=del(lo,1)*cg(lo,1)
   10 continue
      do 6 le=1,nz-1
        do 6 lo=1,nx
          db(lo,le+1)=db(lo,le)+del(lo,le+1)*dg(lo,le+1)
          cb(lo,le+1)=cb(lo,le)+del(lo,le+1)*cg(lo,le+1)
    6 continue
c...    sigma dot computed only at interior interfaces
      do 7 k=1,nz-1
        do 7 lo=1,nx
          dot(lo,k+1)=dot(lo,k)+del(lo,k)
     1               *(db(lo,nz)+cb(lo,nz)-dg(lo,k)-cg(lo,k))
    7 continue
      do 8 k=1,nz
        do 8 lo=1,nx
          vvel(lo,k)=  sl(lo,k)*(cg(lo,k)-cb(lo,nz)-db(lo,nz))-
     1                0.5*(dot(lo,k+1)+dot(lo,k))
          vvel(lo,k)=vvel(lo,k)*ps(lo)
ccc       vvel(lo,k)=vvel(lo,k)*ps(lo)*10.
    8 continue
      return
      end
