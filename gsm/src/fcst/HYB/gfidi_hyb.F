#define CP cp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LONF22 lonf22
#define LONF22P lonf22p
#define LONF lonf
#define LONFP lonfp
#define NCLDG ncldg
#define NTOTAL ntotal
#define NTRAC ntrac
#define OMEGA omega
#define RD rd
#define RERTH rerth
#define CVAP cvap
#define NPES npes
#define JCAP jcap
#define LATG latg

#include <define.h>
      subroutine gfidi_hyb(lons_lat,lat,
     $  ak5,bk5,
     &  dg,tg,zg,ug,vg,rqg,dphi,dlam,qg,
     &  rcl,spdmax,deltim,nvcn,xvcn,
     &  dtdf,dtdl,drdf,drdl,dudl,dvdl,dudf,dvdf,
     &  dqdt,dtdt,drdt,dudt,dvdt)
      implicit none
! sela  add calculation of spdlat as in gfidiu
 
c
#include <paramodel.h>
#ifdef MP
#define LONF22S lonf22p_
#include <commpi.h>
#else
#define LONF22S lonf22_
#endif
c
      integer lons_lat
      integer j,k,n,nvcn,ifirst,lat
      real coriol,rcl,rk,sinra,deltim,xvcn,sinlat
      real
     1    dg(LONF22S,levs_), tg(LONF22S,levs_),  zg(LONF22S,levs_),
     2    ug(LONF22S,levs_), vg(LONF22S,levs_),
     2   rqg(LONF22S,levs_,ntotal_),
     3  dphi(LONF22S), dlam(LONF22S), qg(LONF22S)
      real
     1  dtdf(LONF22S,levs_),       dtdl(LONF22S,levs_),
     1  dudf(LONF22S,levs_),       dudl(LONF22S,levs_),
     1  dvdf(LONF22S,levs_),       dvdl(LONF22S,levs_),
     1  drdf(LONF22S,levs_,ntotal_), drdl(LONF22S,levs_,ntotal_)
      real
     1  dudt(LONF22S,levs_),       dvdt(LONF22S,levs_),
     1  dqdt(LONF22S),
     1  dtdt(LONF22S,levs_),
     1  drdt(LONF22S,levs_,ntotal_), spdmax(levs_)
 
       real pk5(lons_lat,levp1_), dpk(lons_lat,levs_)
      real zadv(lons_lat,levs_,3+ntotal_)
 
      real
     &     dot(lons_lat,levp1_), dotinv(lons_lat,levp1_),
     &      ek(lons_lat,levs_),      cg(lons_lat,levs_),
     &      cb(lons_lat,levs_),      db(lons_lat,levs_),
     &    zlam(lons_lat,levs_),    zphi(lons_lat,levs_),
     &   worka(lons_lat,levs_),   workb(lons_lat,levs_),
     &   workc(lons_lat,levs_),
     &    phiu(lons_lat,levs_),    phiv(lons_lat,levs_),
     &    uprs(lons_lat,levs_),    vprs(lons_lat,levs_),
     &    cofa(lons_lat,levs_),    cofb(lons_lat,levs_),
     &    alfa(lons_lat,levs_),    rlnp(lons_lat,levs_),
     &    px1u(lons_lat,levs_),    px1v(lons_lat,levs_),
     &    px2u(lons_lat,levs_),    px2v(lons_lat,levs_),
     &     px2(lons_lat,levs_),
     &    px3u(lons_lat,levs_),    px3v(lons_lat,levs_),
     &    px4u(lons_lat,levs_),    px4v(lons_lat,levs_),
     &    px5u(lons_lat,levs_),    px5v(lons_lat,levs_),
     &    uphi(lons_lat,levs_),    vphi(lons_lat,levs_),
     &    expq(lons_lat),
     &    rdel(lons_lat,levs_),   rdel2(lons_lat,levs_),
     &  sumdel(lons_lat),          del(lons_lat,levs_),
     &      si(lons_lat,levp1_),     sl(lons_lat,levs_),
     &    rk1,rkr,                 dif(lons_lat)
 
      real cons0,cons0p5,cons1,cons2,clog2   !constant
      real rmin,rmax,delta,delta1
      save clog2,ifirst,delta,delta1
      data ifirst /1/
 
      real ak5(levp1_),bk5(levp1_),dbk(levs_),bkl(levs_),ck(levs_),diff
      integer lons,i

!     print *,' enter gfidi_hyb_fd ' 		! hmhj

      cons0   = 0.d0      !constant
      cons0p5 = 0.5d0     !constant
      cons1   = 1.d0      !constant
      cons2   = 2.d0      !constant
      rk= rd_ /cp_
 
!-------------------------------------------------------
      sinra=sqrt(cons1-cons1/rcl)     !constant
      coriol=cons2*omega_*sinra          !constant
      sinra=sinra/rerth_
 
!-------------------------------------------------------
 
!!sela     rcl = cons1/(cons1-sinlat*sinlat)  !constant
 
      if(ifirst.eq.1)then
 
      clog2=log(cons2)     ! constant
      delta=cvap_/cp_  ! check these cpv cpd (at const p for vapor and dry
      delta1=delta-cons1
      rk1 = rk + 1.e0
      rkr = 1.0/rk
 
      ifirst=0
      endif
c     
      do k=1,levs_
         dbk(k) = bk5(k+1)-bk5(k)
         bkl(k) = (bk5(k+1)+bk5(k))*0.5
         ck(k)  = ak5(k+1)*bk5(k)-ak5(k)*bk5(k+1)
      enddo
c
      do j=1,lons_lat
       expq(j)=exp(qg(j))
      enddo
!.............................................................
! get vertical coordinate for vcnhyb going bot. to top.
      do k=1,levp1_
      do j=1,lons_lat
       si(j,levs_+2-k)= ak5(k)+bk5(k)*expq(j) !ak(k) bk(k) go top to bottom
      enddo
      enddo
 
!sela sumdel=0.
!sela do k=1,levs_
!sela   del(k)= si(k)-si(k+1)
!sela   sumdel=sumdel+del(k)
!sela enddo
 
      do  k=1,levs_
      do j=1,lons_lat
       diff = si(j,k)**rk1 - si(j,k+1)**rk1
       diff = diff / (rk1*(si(j,k)-si(j,k+1)))
       sl(j,k) = diff**rkr
c       sl(j,k) = cons0p5*(si(j,k)+si(j,k+1))
      enddo
      enddo
!.............................................................
 
      do k=1,levp1_
      do j=1,lons_lat
        pk5(j,k)=ak5(k) + bk5(k)*expq(j)
!sela  print*,'sela pk5=',pk5(k),'k=',k
c          if (mype.eq.master) then
c          print'(a4,5f12.4,2i5)','hyb',
c     $         pk5(j,k),expq(j),ak5(k),bk5(k),si(j,k),j,k
c          endif
      enddo
      enddo
 
      do k=1,levs_
      do j=1,lons_lat
         dpk(j,k)=    pk5(j,k+1) - pk5(j,k)
        rdel(j,k)=    cons1/dpk(j,k)            ! constant
       rdel2(j,k)=  cons0p5/dpk(j,k)            ! constant
      enddo
      enddo
 
      k=1
      do j=1,lons_lat
       alfa(j,1)=clog2                          ! constant
!sela  alfa(j,1)=cons1                          ! constant
      enddo
 
      do j=1,lons_lat
        rlnp(j,1)= 99999.99
      enddo
      do  k=2,levs_
      do j=1,lons_lat
        rlnp(j,k)= log( pk5(j,k+1)/pk5(j,k) )
        alfa(j,k)= cons1-( pk5(j,k)/dpk(j,k) )*rlnp(j,k)
      enddo
      enddo
 
      spdmax=0.
      do  k=1,levs_
      do j=1,lons_lat
        ek(j,k)=(ug(j,levp1_-k)*ug(j,levp1_-k)+
     1           vg(j,levp1_-k)*vg(j,levp1_-k))*rcl
      if (ek(j,k) .gt. spdmax(levp1_-k))  spdmax(levp1_-k)=ek(j,k)
      enddo
      enddo
      do k=1,levs_
      do j=1,lons_lat
       cg(j,k)=(ug(j,levp1_-k)*dlam(j)+vg(j,levp1_-k)*dphi(j))*rcl
      enddo
      enddo
      k=1
      do j=1,lons_lat
       db(j,1)=dg(j,levs_)*dpk(j,1)
       cb(j,1)=cg(j,1)*dbk(1)
      enddo
 
      do k=1,levm1_
      do j=1,lons_lat
       db(j,k+1)=db(j,k)+dg(j,levs_-k)*dpk(j,k+1)
       cb(j,k+1)=cb(j,k)+cg(j,k+1)*dbk(k+1)
      enddo
      enddo
 
      do j=1,lons_lat
       dqdt(j)= -db(j,levs_)/expq(j)-cb(j,levs_)
c       dqdt(j)= -cb(j,levs_)
       dot(j,    1)=cons0                    !constant
       dot(j,levp1_)=cons0                    !constant
      enddo
      do k=1,levm1_
      do j=1,lons_lat
        dot(j,k+1)=-expq(j)*(bk5(k+1)*dqdt(j)+cb(j,k)) -db(j,k)
c        dot(j,k+1)=-expq(j)*(bk5(k+1)*(dqdt(j)-db(j,levs_)/expq(j))
c     $        +cb(j,k)) -db(j,k)
      enddo
      enddo
      do k=1,levp1_
      do j=1,lons_lat
       dotinv(j,k)=dot(j,levp1_+1-k)
      enddo
      enddo
 
 
! variables are in bottom to top order  !!!!!!!!!!!!!!!!!
! do horizontal advection.
 
      k=1
      do j=1,lons_lat
        dudt(j,levp1_-k)=-ug(j,levp1_-k)*dudl(j,levp1_-k)
     1                   -vg(j,levp1_-k)*dudf(j,levp1_-k)
        dvdt(j,levp1_-k)=-ug(j,levp1_-k)*dvdl(j,levp1_-k)
     1                   -vg(j,levp1_-k)*dvdf(j,levp1_-k)
        dtdt(j,levp1_-k)=-ug(j,levp1_-k)*dtdl(j,levp1_-k)
     1                   -vg(j,levp1_-k)*dtdf(j,levp1_-k)
      enddo
 
      k=levs_
      do j=1,lons_lat
        dudt(j,levp1_-k)=-ug(j,levp1_-k)*dudl(j,levp1_-k)
     1                   -vg(j,levp1_-k)*dudf(j,levp1_-k)
        dvdt(j,levp1_-k)=-ug(j,levp1_-k)*dvdl(j,levp1_-k)
     1                   -vg(j,levp1_-k)*dvdf(j,levp1_-k)
        dtdt(j,levp1_-k)=-ug(j,levp1_-k)*dtdl(j,levp1_-k)
     1                   -vg(j,levp1_-k)*dtdf(j,levp1_-k)
      enddo
 
      do k=2,levm1_
      do j=1,lons_lat
        dudt(j,levp1_-k)=-ug(j,levp1_-k)*dudl(j,levp1_-k)
     1                   -vg(j,levp1_-k)*dudf(j,levp1_-k)
        dvdt(j,levp1_-k)=-ug(j,levp1_-k)*dvdl(j,levp1_-k)
     1                   -vg(j,levp1_-k)*dvdf(j,levp1_-k)
        dtdt(j,levp1_-k)=-ug(j,levp1_-k)*dtdl(j,levp1_-k)
     1                   -vg(j,levp1_-k)*dtdf(j,levp1_-k)
      enddo
      enddo
c      if (mype.eq.master) then
c      do k=1,1
c        do i=1,5
c        print'(a4,5e12.4,2i5)','dt1',
c     $       dudt(i,k),dvdt(i,k),dtdt(i,k),drdt(i,k,1),dqdt(i),i,lat
c        enddo
c      enddo
c      endif
! add coriolis,deformation note  coriolis sign for s.hemi
      lons=lons_lat/2
      do k=1,levs_
      do j=1,lons
        dudt(j,levp1_-k)=dudt(j,levp1_-k)
     $       +vg(j,levp1_-k)*coriol
 
        dudt(j+lons,levp1_-k)=dudt(j+lons,levp1_-k)
     $       -vg(j+lons,levp1_-k)*coriol
 
        dvdt(j,levp1_-k)=dvdt(j,levp1_-k)
     $       -ug(j,levp1_-k)*coriol
     &       -sinra*ek(j,k)

        dvdt(j+lons,levp1_-k)=dvdt(j+lons,levp1_-k)
     $       +ug(j+lons,levp1_-k)*coriol
     &       +sinra*ek(j+lons,k)
      enddo
      enddo
c      if (mype.eq.master) then
c      do k=1,1
c        do i=1,5
c        print'(a4,5e12.4,2i5)','dt2',
c     $       dudt(i,k),dvdt(i,k),dtdt(i,k),drdt(i,k,1),dqdt(i),i,lat
c        enddo
c      enddo
c      endif
!.................................................................
! calculate pressure force:
       k=1
      do j=1,lons_lat
        cofb(j,k)=-rdel(j,k)*(                 alfa(j,k)*dbk(k))
      enddo
 
      do k=2,levs_
      do j=1,lons_lat
        cofb(j,k)=-rdel(j,k)*(bk5(k)*rlnp(j,k)+alfa(j,k)*dbk(k))
      enddo
      enddo
 
      do k=1,levs_
      do j=1,lons_lat
        uprs(j,k)=cofb(j,k)*rd_*tg(j,levp1_-k)*expq(j)*dlam(j)
        vprs(j,k)=cofb(j,k)*rd_*tg(j,levp1_-k)*expq(j)*dphi(j)
      enddo
      enddo
      do k=1,levs_
      do j=1,lons_lat
        cofa(j,k)=-rdel(j,k)*(
     &   bk5(k+1)*pk5(j,k)/pk5(j,k+1) - bk5(k)
     &  +rlnp(j,k)*( bk5(k)-pk5(j,k)*dbk(k)*rdel(j,k) )  )
      enddo
      enddo
!.................................................................
      do k=1,levs_
      do j=1,lons_lat
        px1u(j,k)=cons0              ! grid topography =0 for testing
        px1v(j,k)=cons0              ! grid tpopgraphy =0 for testing
      enddo
      enddo
! see programming notes for looping in  calculating  px2 and px3
      do j=1,lons_lat
        px2(j,levs_)=cons0                             ! constant
        px2(j,levm1_)=
     &  -rd_*( bk5(levp1_)/pk5(j,levp1_)-bk5(levs_)/pk5(j,levs_) )
     $       *tg(j,1)
      enddo
 
 
      do k=2,levm1_
      do j=1,lons_lat
      px2(j,levs_-k)=px2(j,levp1_-k)
     & -rd_*(bk5(levs_+2-k)/
     $        pk5(j,levs_+2-k)-bk5(levp1_-k)/pk5(j,levp1_-k))*
     &                                                  tg(j,k)
      enddo
      enddo
 
      do k=1,levs_
      do j=1,lons_lat
        px2u(j,k)=px2(j,k)*expq(j)*dlam(j)
        px2v(j,k)=px2(j,k)*expq(j)*dphi(j)
      enddo
      enddo
      do j=1,lons_lat
       px3u(j,levs_)=cons0 ! constant
       px3v(j,levs_)=cons0 ! constant
      enddo
 
      do j=1,lons_lat
       px3u(j,levm1_)=-rd_*rlnp(j,levs_)*dtdl(j,1)
       px3v(j,levm1_)=-rd_*rlnp(j,levs_)*dtdf(j,1)
      enddo
 
      do k=2,levm1_
      do j=1,lons_lat
       px3u(j,levs_-k)=px3u(j,levp1_-k)-rd_*rlnp(j,levp1_-k)*dtdl(j,k)
       px3v(j,levs_-k)=px3v(j,levp1_-k)-rd_*rlnp(j,levp1_-k)*dtdf(j,k)
      enddo
      enddo
      do k=1,levs_
      do j=1,lons_lat
       px3u(j,k)=px3u(j,k)/rcl
       px3v(j,k)=px3v(j,k)/rcl
      enddo
      enddo
      do k=1,levs_
      do j=1,lons_lat
        px4u(j,k)=-rd_*alfa(j,k)*dtdl(j,levp1_-k)/rcl
        px4v(j,k)=-rd_*alfa(j,k)*dtdf(j,levp1_-k)/rcl
      enddo
      enddo
 
      do k=1,levs_
      do j=1,lons_lat
        px5u(j,k)=-cofa(j,k)*rd_*tg(j,levp1_-k)*expq(j)*dlam(j)
        px5v(j,k)=-cofa(j,k)*rd_*tg(j,levp1_-k)*expq(j)*dphi(j)
      enddo
      enddo
 
      do k=1,levs_
      do j=1,lons_lat
        uphi(j,k)=px1u(j,k)+px2u(j,k)+px3u(j,k)+px4u(j,k)+px5u(j,k)
        vphi(j,k)=px1v(j,k)+px2v(j,k)+px3v(j,k)+px4v(j,k)+px5v(j,k)
      enddo
      enddo
      do k=1,levs_
      do j=1,lons_lat
        dudt(j,levp1_-k)=dudt(j,levp1_-k)+uphi(j,k)+uprs(j,k)
        dvdt(j,levp1_-k)=dvdt(j,levp1_-k)+vphi(j,k)+vprs(j,k)
      enddo
      enddo
!..............................................................

c      if (mype.eq.master) then
c      do k=1,1
c        do i=1,5
c        print'(a4,5e12.4,2i5)','dt3',
cc     $       dudt(i,k),dvdt(i,k),dtdt(i,k),drdt(i,k,1),dqdt(i),i,k
cc     $       uphi(i,k),uprs(i,k),vphi(i,k),vprs(i,k),dqdt(i),i,k
c     $       px2(i,k),dlam(i),rlnp(i,k),dtdl(i,k),tg(i,k),i,lat
c        enddo
c      enddo
c      endif
      do k=1,levs_
      do j=1,lons_lat
        worka(j,k)=rk*tg(j,levp1_-k)/(cons1+(delta1)*
     &   rqg(j,levp1_-k,1)) * rdel(j,k)
      enddo
      enddo
 
      k=1
      do j=1,lons_lat
       workb(j,1)=
     &  alfa(j,1)*( dg(j,levs_)*dpk(j,1)+expq(j)*cb(j,1)*dbk(1) )
      enddo
 
      do k=2,levs_
      do j=1,lons_lat
        workb(j,k)=rlnp(j,k)*( db(j,k-1)+expq(j)*cb(j,k-1) )
     &  +alfa(j,k)*( dg(j,levp1_-k)*dpk(j,k)+expq(j)*cg(j,k)*dbk(k) )
      enddo
      enddo
 
      k=1
      do j=1,lons_lat
       workc(j,1)=expq(j)*cg(j,1)*dbk(1)
      enddo
 
      do k=2,levs_
      do j=1,lons_lat
        workc(j,k)=expq(j)*cg(j,k)*( dbk(k)+ck(k)*rlnp(j,k)*rdel(j,k) )
      enddo
      enddo
 
      do k=1,levs_
      do j=1,lons_lat
       dtdt(j,levp1_-k)=
     & dtdt(j,levp1_-k)+worka(j,k)*( -workb(j,k) + workc(j,k))
      enddo
      enddo
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do 330 n=1,ntotal_
      k=1
      do j=1,lons_lat
        drdt(j,levp1_-k,n)=-ug(j,levp1_-k)*drdl(j,levp1_-k,n)
     1                     -vg(j,levp1_-k)*drdf(j,levp1_-k,n)
      enddo
      k=levs_
      do j=1,lons_lat
        drdt(j,levp1_-k,n)=-ug(j,levp1_-k)*drdl(j,levp1_-k,n)
     1                     -vg(j,levp1_-k)*drdf(j,levp1_-k,n)
      enddo
      do k=2,levm1_
      do j=1,lons_lat
        drdt(j,levp1_-k,n)=-ug(j,levp1_-k)*drdl(j,levp1_-k,n)
     1                     -vg(j,levp1_-k)*drdf(j,levp1_-k,n)
      enddo
      enddo
330   continue
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
! do vertical advection
 
      k=1
      do j=1,lons_lat
        zadv(j,levp1_-k,1)=
     1 -rdel2(j,k)*dot(j,k+1)*( ug(j,levs_-k)-ug(j,levp1_-k))
 
        zadv(j,levp1_-k,2)=
     1 -rdel2(j,k)*dot(j,k+1)*( vg(j,levs_-k)-vg(j,levp1_-k))
 
        zadv(j,levp1_-k,3)=
     1 -rdel2(j,k)*dot(j,k+1)*( tg(j,levs_-k)-tg(j,levp1_-k))
 
      enddo
 
      k=levs_
      do j=1,lons_lat
        zadv(j,levp1_-k,1)=
     1 -rdel2(j,k)*dot(j,k)*( ug(j,levp1_-k)-ug(j,levs_+2-k) )
 
        zadv(j,levp1_-k,2)=
     1 -rdel2(j,k)*dot(j,k)*( vg(j,levp1_-k)-vg(j,levs_+2-k) )
 
        zadv(j,levp1_-k,3)=
     1 -rdel2(j,k)*dot(j,k)*( tg(j,levp1_-k)-tg(j,levs_+2-k) )
 
      enddo
 
      do k=2,levm1_
      do j=1,lons_lat
        zadv(j,levp1_-k,1)=
     1 -rdel2(j,k)*( dot(j,k+1)*( ug(j,levs_  -k)-ug(j,levp1_-k) ) +
     2               dot(j,k  )*( ug(j,levp1_-k)-ug(j,levs_+2-k) ) )
 
        zadv(j,levp1_-k,2)=
     1 -rdel2(j,k)*( dot(j,k+1)*( vg(j,levs_  -k)-vg(j,levp1_-k) ) +
     2               dot(j,k  )*( vg(j,levp1_-k)-vg(j,levs_+2-k) ) )
 
        zadv(j,levp1_-k,3)=
     1 -rdel2(j,k)*( dot(j,k+1)*( tg(j,levs_  -k)-tg(j,levp1_-k) ) +
     2               dot(j,k  )*( tg(j,levp1_-k)-tg(j,levs_+2-k) ) )
 
      enddo
      enddo
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      do 340 n=1,ntotal_
      k=1
      do j=1,lons_lat
        zadv(j,levp1_-k,3+n)=
     1 -rdel2(j,k)*dot(j,k+1)*( rqg(j,levs_-k,n)-rqg(j,levp1_-k,n) )
      enddo
 
      k=levs_
      do j=1,lons_lat
        zadv(j,levp1_-k,3+n)=
     1 -rdel2(j,k)*dot(j,k)*( rqg(j,levp1_-k,n)-rqg(j,levs_+2-k,n) )
      enddo
 
      do k=2,levm1_
      do j=1,lons_lat
        zadv(j,levp1_-k,3+n)=-rdel2(j,k)*
     1 ( dot(j,k+1)*( rqg(j,levs_  -k,n)-rqg(j,levp1_-k,n) ) +
     2   dot(j,k  )*( rqg(j,levp1_-k,n)-rqg(j,levs_+2-k,n) ) )
      enddo
      enddo
340   continue
c      if (mype.eq.master) then
c      do k=1,1
c      do i=1,5
c        print'(a5,4e12.4,2i5)','zad1',
c     $       zadv(i,k,1),zadv(i,k,2),zadv(i,k,3),zadv(i,k,4),i,lat
c      enddo
c      enddo
c      endif
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      call vcnhyb(lons_lat,levs_,3+ntotal_,deltim,
     &            si,sl,dotinv,zadv,nvcn,xvcn)
!sela if(xvcn.ne.0.) print*,'xvcn=',xvcn,' nvcn=',nvcn
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 
! add vertical filterd advection
c      if (mype.eq.master) then
c      do k=1,1
c      do i=1,5
c        print'(a5,4e12.4,2i5)','zad2',
c     $       zadv(i,k,1),zadv(i,k,2),zadv(i,k,3),zadv(i,k,4),i,lat
c      enddo
c      enddo
c      endif
      do k=1,levs_
      do j=1,lons_lat
       dudt(j,k)=dudt(j,k)+zadv(j,k,1)
       dvdt(j,k)=dvdt(j,k)+zadv(j,k,2)
       dtdt(j,k)=dtdt(j,k)+zadv(j,k,3)
      enddo
      enddo
      do  n=1,ntotal_
       do k=1,levs_
       do j=1,lons_lat
        drdt(j,k,n)=drdt(j,k,n)+zadv(j,k,3+n)
       enddo
       enddo
      enddo
 
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! this multiplication must be on  completed tendencies.
      do k=1,levs_
      do j=1,lons_lat
        dudt(j,levp1_-k)=dudt(j,levp1_-k)*rcl
        dvdt(j,levp1_-k)=dvdt(j,levp1_-k)*rcl
      enddo
      enddo

c      if (mype.eq.master) then
c      do k=1,1
c        do i=1,5
c        print'(a4,5e12.4,2i5)','dt4',
c     $       dudt(i,k),dvdt(i,k),dtdt(i,k),drdt(i,k,1),
c     $          dqdt(i)+db(i,levs_)/expq(i),i,lat
c        enddo
c      enddo
c      endif


!     print *,' leave gfidi_hyb_fd ' 		! hmhj

      return
      end
 
      subroutine vcnhyb(im,km,nm,dt,zint,zmid,zdot,zadv,nvcn,xvcn)
c                .      .    .                                       .
c subprogram:    vcnhyb      vertical advection instability filter
c   prgmmr: iredell          org: w/nmc23    date: 91-05-07
c
c abstract: filters vertical advection tendencies
c   in the dynamics tendency equation in order to ensure stability
c   when the vertical velocity exceeds the cfl criterion.
c   the vertical velocity in this case is sigmadot.
c   for simple second-order centered eulerian advection,
c   filtering is needed when vcn=zdot*dt/dz>1.
c   the maximum eigenvalue of the linear advection equation
c   with second-order implicit filtering on the tendencies
c   is less than one for all resolvable wavenumbers (i.e. stable)
c   if the nondimensional filter parameter is nu=(vcn**2-1)/4.
c
c program history log:
c   97-07-30  iredell
c
c usage:    call vcnhyb(im,km,nm,dt,zint,zmid,zdot,zadv,nvcn,xvcn)
c
c   input argument list:
c     im       - integer number of gridpoints to filter
c     km       - integer number of vertical levels
c     nm       - integer number of fields
c     dt       - real timestep in seconds
c     zint     - real (im,km+1) interface vertical coordinate values
c     zmid     - real (im,km) midlayer vertical coordinate values
c     zdot     - real (im,km+1) vertical coordinate velocity
c     zadv     - real (im,km,nm) vertical advection tendencies
c
c   output argument list:
c     zadv     - real (im,km,nm) vertical advection tendencies
c     nvcn     - integer number of points requiring filtering
c     xvcn     - real maximum vertical courant number
c
c   subprograms called:
c     tridim_hyb   - tridiagonal matrix solver
c
      implicit none
      integer,intent(in):: im,km,nm
      real,intent(in):: dt,zint(im,km+1),zmid(im,km),zdot(im,km+1)
      real,intent(inout):: zadv(im,km,nm)
      integer,intent(out):: nvcn
      real,intent(out):: xvcn
      integer i,j,k,n,ivcn(im)
      logical lvcn(im)
      real zdm,zda,zdb,vcn(im,km-1)
      real rnu,cm(im,km),cu(im,km-1),cl(im,km-1)
      real rr(im,km,nm)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  compute vertical courant number
c  increase by 10% for safety
      nvcn=0
      xvcn=0.
      lvcn=.false.
      do k=1,km-1
        do i=1,im
          zdm=zmid(i,k)-zmid(i,k+1)
          vcn(i,k)=abs(zdot(i,k+1)*dt/zdm)*1.1
          lvcn(i)=lvcn(i).or.vcn(i,k).gt.1
          xvcn=max(xvcn,vcn(i,k))
        enddo
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  determine points requiring filtering
      if(xvcn.gt.1) then
        do i=1,im
          if(lvcn(i)) then
            ivcn(nvcn+1)=i
            nvcn=nvcn+1
          endif
        enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  compute tridiagonal matrim
        do j=1,nvcn
          cm(j,1)=1
        enddo
        do k=1,km-1
          do j=1,nvcn
            i=ivcn(j)
            if(vcn(i,k).gt.1) then
             zdm=zmid(i,k)-zmid(i,k+1)
             zda=zint(i,k+1)-zint(i,k+2)
             zdb=zint(i,k)-zint(i,k+1)
              rnu=(vcn(i,k)**2-1)/4
              cu(j,k)=-rnu*zdm/zdb
              cl(j,k)=-rnu*zdm/zda
              cm(j,k)=cm(j,k)-cu(j,k)
              cm(j,k+1)=1-cl(j,k)
            else
              cu(j,k)=0
              cl(j,k)=0
              cm(j,k+1)=1
            endif
          enddo
        enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  fill fields to be filtered
        do n=1,nm
          do k=1,km
            do j=1,nvcn
              i=ivcn(j)
              rr(j,k,n)=zadv(i,k,n)
            enddo
          enddo
        enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  solve tridiagonal system
        call tridim_hyb(nvcn,im,km,km,nm,cl,cm,cu,rr,cu,rr)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  replace filtered fields
        do n=1,nm
          do k=1,km
            do j=1,nvcn
              i=ivcn(j)
              zadv(i,k,n)=rr(j,k,n)
            enddo
          enddo
        enddo
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
c-----------------------------------------------------------------------
      subroutine tridim_hyb(l,lx,n,nx,m,cl,cm,cu,r,au,a)
c                .      .    .                                       .
c subprogram:    tridim_hyb      solves tridiagonal matrix problems.
c   prgmmr: iredell          org: w/nmc23    date: 91-05-07
c
c abstract: this routine solves multiple tridiagonal matrix problems
c   with multiple right-hand-side and solution vectors for every matrix.
c   the solutions are found by eliminating off-diagonal coefficients,
c   marching first foreward then backward along the matrix diagonal.
c   the computations are vectorized around the number of matrices.
c   no checks are made for zeroes on the diagonal or singularity.
c
c program history log:
c   97-07-30  iredell
c
c usage:    call tridim_hyb(l,lx,n,nx,m,cl,cm,cu,r,au,a)
c
c   input argument list:
c     l        - integer number of tridiagonal matrices
c     lx       - integer first dimension (lx>=l)
c     n        - integer order of the matrices
c     nx       - integer second dimension (nx>=n)
c     m        - integer number of vectors for every matrix
c     cl       - real (lx,2:n) lower diagonal matrix elements
c     cm       - real (lx,n) main diagonal matrix elements
c     cu       - real (lx,n-1) upper diagonal matrix elements
c                (may be equivalent to au if no longer needed)
c     r        - real (lx,nx,m) right-hand-side vector elements
c                (may be equivalent to a if no longer needed)
c
c   output argument list:
c     au       - real (lx,n-1) work array
c     a        - real (lx,nx,m) solution vector elements
c
c attributes:
c   language: fortran 77.
c   machine:  cray.
c
      real cl(lx,2:n),cm(lx,n),cu(lx,n-1),r(lx,nx,m),
     &                         au(lx,n-1),a(lx,nx,m)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  march up
      do i=1,l
        fk=1./cm(i,1)
        au(i,1)=fk*cu(i,1)
      enddo
      do j=1,m
        do i=1,l
          fk=1./cm(i,1)
          a(i,1,j)=fk*r(i,1,j)
        enddo
      enddo
      do k=2,n-1
        do i=1,l
          fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)=fk*cu(i,k)
        enddo
        do j=1,m
          do i=1,l
            fk=1./(cm(i,k)-cl(i,k)*au(i,k-1))
            a(i,k,j)=fk*(r(i,k,j)-cl(i,k)*a(i,k-1,j))
          enddo
        enddo
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  march down
      do j=1,m
        do i=1,l
          fk=1./(cm(i,n)-cl(i,n)*au(i,n-1))
          a(i,n,j)=fk*(r(i,n,j)-cl(i,n)*a(i,n-1,j))
        enddo
      enddo
      do k=n-1,1,-1
        do j=1,m
          do i=1,l
            a(i,k,j)=a(i,k,j)-au(i,k)*a(i,k+1,j)
          enddo
        enddo
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      end
