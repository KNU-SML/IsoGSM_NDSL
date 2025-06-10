#define JCAP jcap
#define LNT2 lnt2
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LEVSP levsp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NCLDG ncldg
#define NPES npes
#define NROW nrow
#define NTRAC ntrac
#define NTOTAL ntotal
#define G g
#define PI pi

#include <define.h>
#ifdef NISLQ
c-------------------------------------------------------------------------------
         subroutine nislq_init(nsize,myrank,colradd,rbs22)
c-------------------------------------------------------------------------------
#include <paramodel.h>
#include <constant.h>
#ifdef MP
#include <commpi.h>
#else
c          integer
c           common /comfgrid/ 
c     *     latdef(latg_/2)
c
c           logical iope
c           common /combitl/ iope
#include <comfgrid.h>
#include <comio.h>
#endif
#include <comnislq.h>
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
#ifdef MP
#define LEVSS levsp_
#else
#define LEVSS levs_
#endif
c
c passing variables
c
           integer,intent(in)                    :: nsize, myrank
           real   ,intent(in), dimension(latg2_) :: colradd,rbs22
           integer :: ncld
c
c local variables
c
         integer                               :: myh,my2,i,j,ii,jj,js,k
           real                                  :: hfpi,twopi,dlat,dlon
           integer                               :: n,nm,nr,lonp
           real, dimension(latg2_)               :: colatrad
           integer                               :: j1,j2

          
c
c assign parameter
c
           nx      = lonf_
           my      = latg_
           lev     = levs_
           ncld    = ntotal_
           nlevs   = ncld*lev
           nlevsp  = ncld*LEVSS
c
           lonfull = nx
           lonhalf = nx / 2            
           lonpart = lonhalf/nsize+1
           latfull = my * 2
           lathalf = my 
           latpart = 2*((lathalf/2)/nsize+1)
#ifdef MP
           my_max  = latpart
#else
           my_max  = latg_
#endif
c
c constant
c
           hfpi  = pi_ * 0.5
           twopi = pi_ * 2.0
           myh   = my/2
           my2   = my*2
c
c initialize slq_q
c
           slq_q1=0.
           slq_q2=0.
           slq_q3=0.
           slq_psfc2=0.
           slq_u2=0.
           slq_v2=0.
           slq_w2=0.
c
c informations for latitude
c define latitude
c
           do j = 1,latg2_
#ifdef MP
             jlist1(j)=j
#else
             jlist1(j)=j
#endif
           enddo
c
c define colatitude in radian
c
           do j = 1,latg2_
             colatrad(j)=colradd(j)
           enddo
c
c ------------------- gaussian latitude 
c
c
           do j = 1,myh
             glat(j)=hfpi-colatrad(j)
             glat(my+1-j)=-glat(j)
           enddo
c      c co-latitude
           do j = 1,myh 
             gglat(j) = hfpi - glat(j)
             gglat(my+1-j) = hfpi + glat(j)
           enddo
           do j = my+1,my2
             gglat(j) = twopi - gglat(my2+1-j)
           enddo
c
           gglati(myh+1) = hfpi
           dlat = (gglati(myh+1)-gglat(myh))*2.0
           do j = myh,2,-1
             gglati(j) = gglati(j+1) - dlat
             dlat      = (gglati(j) - gglat(j-1))*2.0
           enddo
           gglati(1) = 0.0
           do j = myh+2,my
             gglati(j) = pi_ - gglati(my+2-j)
           enddo
           gglati(my+1) = pi_
           do j = my+2,my2
             gglati(j) = twopi - gglati(my2+2-j)
           enddo
           gglati(my2+1) = twopi
#ifdef SLDBG
c
           if( iope ) then
             print *,' j interface =1    gglati=',gglati(1)
             do j = 1,my2
               print *,'               j mean =',j,'  gglat= ',gglat(j)
               print *,' j interface =',j+1,' gglati= ',gglati(j+1)
             enddo
           endif
#endif
c
c  determin the longitude with full grid
c
           dlon = twopi / nx
           do i = 1,nx
             gglon(i)=(i-1)*dlon
           enddo
           do i = 2,nx
             ggloni(i)=0.5*(gglon(i-1)+gglon(i))
           enddo
           ggloni(     1)=ggloni(   2)-dlon
           ggloni(nx+1)=ggloni(nx)+dlon
c
#ifdef SLDBG
           if( iope ) then
             print *,' ------ total edge nx=',nx,' -------'
             print *,' i edge number=1    ggloni=',ggloni(1)
             do i = 1,nx
        print *,'               i cell number=',i,'  gglon= ',gglon(i)
               print *,' i edge number=',i+1,' ggloni= ',ggloni(i+1)
             enddo
           endif
#endif
c
c --------------------- for parallel --------------------
c
c in nisl, we do great circle, so nx and my2 are full
c transpose will between (lonfull,lev,latpart)  by (nx  ,lev,my/nsize+1)
c                    and (latfull,lev,lonpart)  by (my*2,lev,nx/2/nsize+1)
c
c
c equally distribute len, no location is considered
c
           sllonlen(1:nsize) = 0
           i=1
           do ii = 1,lonpart
             do n = 1,nsize
               if(i.le.lonhalf) then
          sllonlen(n) = sllonlen(n)+1
                 i=i+1
               endif
             enddo
           enddo
c
c sequential location for longitude
c
           nm=1
           do n = 1,nsize
             sllonstr(n) = nm
             nm = nm + sllonlen(n)
           enddo
c
c check make_list  to have consistent latitude number for each pe here
c
c refer equdis
           sllatlen(1:nsize) = 0
           i=1
           n=1
           do jj = 1,myh
             sllatlen(n)=sllatlen(n)+1
             n=n+i
             if(n.eq.nsize+1) then
               i=-1
               n=n+i
             endif
             if(n.eq.0) then
               i=1
               n=n+i
             endif
           enddo
           sllatlen=sllatlen*2
#ifdef SLDBG
           if(iope) then
             do n = 1,nsize
               print *,'n,sllatlen=',n,sllatlen(n)
             enddo
           endif
#endif
c
c sequential location for latitude
c
           nm=1
           do n = 1,nsize
             sllatstr(n) = nm
             nm = nm + sllatlen(n)
           enddo
c
c   jlistnum=sllatlen(myrank+1)
c   lat_s=sllatstr(myrank+1)
c   lat_e=sllatstr(myrank+1)+sllatlen(myrank+1)-1
#ifdef SLDBG
c
c check by print
c
           if( iope ) then
             do n = 1,nsize
              print *,' pe sllonstr sllonlen ',n,sllonstr(n),sllonlen(n)
         print *,' pe sllatstr sllatlen ',n,sllatstr(n),sllatlen(n)
             enddo
           endif
#endif
c
c make true latitude index follow make_list
c
cc   j=1
cc   do jj=1,my/nsize+1
cc     js = jj
cc     do n=1,nsize
cc       if( j.le.my ) then
cc         truej (js) = j
cc         shflj (j ) = js
cc         j = j + 1
cc       endif
cc       js = js + sllatlen(n)
cc     enddo
cc   enddo
           do jj = 1,latg2_
             j1=jj*2-1
             j2=jj*2
             truej(j1)=latdef(jj)
             truej(j2)=latg_-latdef(jj)+1
           enddo
#ifdef SLDBG
c
c check by print
c
           if(iope) then
c     do j=1,my
c       print *,' true j=',j,' to shaffled j=',shflj(j),   &
c                 ' back to true j=',truej(shflj(j))
c     enddo
             do j = 1,my
               print *,' shfl j=',j,' to true j=',truej(j),
     &                 ' back to shfl j=',shflj(truej(j))
             enddo
             j=0
             do n = 1,nsize
               print *,' --- start pe =',n-1
               do jj = 1,sllatlen(n)
                 j=j+1
              print *,' shaffled j=',j,' to true j=',truej(j)
               enddo
             enddo
           endif
#endif
c
           mysllonlen = sllonlen(myrank+1)
           mysllatlen = sllatlen(myrank+1)
c debug
c  if( n.eq.n ) then
c    call mpe_finalize
c    stop
c  endif
c
           return
           end subroutine nislq_init
c-------------------------------------------------------------------------------
c
c-------------------------------------------------------------------------------
           subroutine cyclic_cell_massadvx(levs,ncld,delt,uc,qq,mass)
c-------------------------------------------------------------------------------
c
c compute local positive advection with mass conservation
c qq is advected by uc which is in radiance/sec from past to next position
c
c author: hann-ming henry juang 2008

#include <paramodel.h>
#include <comnislq.h>
c
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
           real   , parameter                    ::  fa1 = 9./16.     
     &,                                              fa2 = 1./16.
           integer                               ::  levs,mass
           integer                               ::  i,k,im,n,ncld
           real                                  ::  delt, rm2, dist, sc
           real   , dimension(lonfull,levs)      ::  uc
           real   , dimension(lonfull,levs,ncld) ::  qq
           real   , dimension(lonfull)           ::  past,next,da,dxfact
           real   , dimension(lonfull+1)         ::  xpast,xnext,uint
c
c preparations ---------------------------
c
c x is equal grid spacing, so location can be specified by grid point number
c
           im = lonfull
c
      
           do k = 1,levs
c
c 4th order interpolation from mid point to cell interfaces
c
             do i = 3,im-1
               uint(i)=fa1*(uc(i,k)+uc(i-1,k))-fa2*(uc(i+1,k)+uc(i-2,k))
             enddo
             uint(2)=fa1*(uc(2,k)+uc(1 ,k))-fa2*(uc(3,k)+uc(im  ,k))
             uint(1)=fa1*(uc(1,k)+uc(im,k))-fa2*(uc(2,k)+uc(im-1,k))
             uint(im+1)=uint(1)
          uint(im  )=fa1*(uc(im,k)+uc(im-1,k)) -fa2*(uc(1,k)+uc(im-2,k))
c
c compute past and next positions of cell interfaces
c
             do i = 1,im+1
               dist     = uint(i) * delt
               xpast(i) = ggloni(i) - dist
               xnext(i) = ggloni(i) + dist
             enddo
c      
             if( mass.eq.1 ) then
               do i = 1,im
               dxfact(i) = (xpast(i+1)-xpast(i)) / (xnext(i+1)-xnext(i))
               enddo
             endif
c
c  mass positive advection
c
             sc=ggloni(im+1)-ggloni(1)
             do n = 1,ncld
               past(1:im) = qq(1:im,k,n)
             call cyclic_cell_ppm_intp(ggloni,past,xpast,da,im,im,im,sc)
               if(mass.eq.1) da(1:im) = da(1:im) * dxfact(1:im)
             call cyclic_cell_ppm_intp(xnext,da,ggloni,next,im,im,im,sc)
               qq(1:im,k,n) = next(1:im)
             enddo
c
           enddo
c
           return
           end subroutine cyclic_cell_massadvx
c-------------------------------------------------------------------------------
c
c-------------------------------------------------------------------------------
           subroutine cyclic_cell_massadvy(levs,ncld,delt,vc,qq,mass)
c-------------------------------------------------------------------------------
c
c compute local positive advection with mass conserving
c qq will be advect by vc from past to next location with 2*delt
c
c author: hann-ming henry juang 2007
c
c-------------------------------------------------------------------------------

#include <paramodel.h>
#include <comnislq.h>
c-------------------------------------------------------------------------------
           real   , parameter                    ::  fa1 = 9./16.   
     &,                                              fa2 = 1./16.
           integer                               ::  levs,mass  
           integer                                   n,k,j,jm,jm2, ncld
           real                                  ::  delt,sc
           real   , dimension(latfull,levs)      ::  vc
           real   , dimension(latfull,levs,ncld) ::  qq
           real   , dimension(latfull)           ::  var,past,da,next
     &,                                              dyfact
           real   , dimension(latfull+1)         ::  ypast,ynext,dist
c
c preparations ---------------------------
c
           jm   = lathalf
           jm2  = latfull
c
           do k = 1,levs
c
             do j = 1,jm
               var(j)      =-vc(j   ,k) * delt
               var(j+jm)   = vc(j+jm,k) * delt
             enddo
c
             do j = 3,jm2-1
               dist(j)=fa1*(var(j)+var(j-1))-fa2*(var(j+1)+var(j-2))
             enddo
             dist(2)=fa1*(var(2)+var(1  ))-fa2*(var(3)+var(jm2  ))
             dist(1)=fa1*(var(1)+var(jm2))-fa2*(var(2)+var(jm2-1))
             dist(jm2+1)=dist(1)
        dist(jm2  )=fa1*(var(jm2)+var(jm2-1))-fa2*(var(1)+var(jm2-2))
c 
             do j = 1,jm2+1
               ypast(j) = gglati(j) - dist(j)
               ynext(j) = gglati(j) + dist(j)
             enddo
c
             if( mass.eq.1 ) then
               do j = 1,jm2
          dyfact(j) = (ypast(j+1)-ypast(j)) / (ynext(j+1)-ynext(j))
               enddo
             endif
c
c advection all in y
c
             sc=gglati(jm2+1)-gglati(1)
             do n = 1,ncld
               past(1:jm2) = qq(1:jm2,k,n)
          call cyclic_cell_ppm_intp(gglati,past,ypast,da,jm2,jm2,jm2,sc)
               if( mass.eq.1 ) da(1:jm2) = da(1:jm2) * dyfact(1:jm2)
          call cyclic_cell_ppm_intp(ynext,da,gglati,next,jm2,jm2,jm2,sc)
               qq(1:jm2,k,n) = next(1:jm2)
             enddo
c 
           enddo
c
           return
           end subroutine cyclic_cell_massadvy
c-------------------------------------------------------------------------------
c
c-------------------------------------------------------------------------------
           subroutine cyclic_cell_intpx(levs,imp,imf,qq)
c-------------------------------------------------------------------------------
c
c do  mass conserving interpolation from different grid at given latitude
c
c author: hann-ming henry juang 2008
c
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
#include <paramodel.h>
#include <comnislq.h>


           integer                          ::  levs, imp, imf, i,k,im
      real                             ::  two_pi,dxp,dxf,hfdxp,hfdxf,sc
           real   , dimension(lonfull,levs) ::  qq
           real   , dimension(lonfull)      ::  past,next
           real   , dimension(lonfull+1)    ::  xpast,xnext
c
           im = lonfull
c ..................................  
           if( imp.ne.imf ) then
c ..................................
             two_pi = 2.*pi_
             dxp = two_pi / imp
             dxf = two_pi / imf
             hfdxp = 0.5 * dxp
             hfdxf = 0.5 * dxf
c
             do i = 1,imp+1
               xpast(i) = (i-1) * dxp - hfdxp
             enddo
c
             do i = 1,imf+1
               xnext(i) = (i-1) * dxf - hfdxf
             enddo
c
             sc=two_pi
             do k = 1,levs
               do i = 1,imp
                 past(i)=qq(i,k)
               enddo
     
      call cyclic_cell_ppm_intp(xpast,past,xnext,next,im,imp,imf,sc)
               do i = 1,imf
                 qq(i,k)=next(i)
               enddo
             enddo
c .................       
           endif
c .................
           return
           end subroutine cyclic_cell_intpx
c-------------------------------------------------------------------------------
c
c-------------------------------------------------------------------------------
         subroutine cyclic_cell_ppm_intp(pp,qq,pn,qn,lons,lonp,lonn,sc)
c-------------------------------------------------------------------------------
c
c mass conservation in cyclic bc interpolation: interpolate a group
c of grid point  coordiante call pp at interface with quantity qq at
c cell averaged to a group of new grid point coordinate call pn at
c interface with quantity qn at cell average with ppm spline.
c in horizontal with mass conservation is under the condition that
c pp(1)= pp(lons+1)=pn(lons+1)
c
c pp    location at interfac level as input
c qq    quantity at averaged-cell as input
c pn    location at interface of new grid structure as input
c qn    quantity at averaged-cell as output
c lons  numer of cells for dimension
c lonp  numer of cells for input
c lonn  numer of cells for output
c lev   number of vertical layers
c mono  monotonicity o:no, 1:yes
c
c author : henry.juang@noaa.gov
c
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
           integer               ::  lons,lonp,lonn            
           integer               ::  ik,le,kstr,kend           
           integer               ::  i,k, kl, kh, kk, kkl, kkh
           real                  ::  length,sc                 
           real                 ::  dqi,dqimax,dqimin         
           real               ::           tl,tl2,tl3,qql,dql 
           real             ::     th,th2,th3,qqh,dqh         
           real            ::      dpp,dpq,c1,c2              
           real            ::      dpph,dppl,dpqh,dpql
           real   , dimension(lons)    ::  qq,qn
         real   , dimension(lons+1)  ::  pp,pn
           real   , dimension(3*lonp)  ::  locs,mass,hh,dqmono,qmi,qpi
c
c  length = pp(lonp+1) - pp(1)
           length=sc
c
c arrange input array cover output location with cyclic boundary condition
c
           locs(lonp+1:2*lonp) = pp(1:lonp)
           do i = 1,lonp
             locs(i) = locs(i+lonp) - length
             locs(i+2*lonp) = locs(i+lonp) + length
           enddo
c
           find_kstr : do i = 1,3*lonp
             if( pn(1).ge.locs(i) .and. pn(1).lt.locs(i+1) ) then
               kstr = i
               exit find_kstr
             else
               cycle find_kstr
             endif
           enddo find_kstr
           kstr=max(1,kstr)
c
           mass(lonp+1:2*lonp) = qq(1:lonp)
           do i = 1,lonp
             mass(i) = mass(i+lonp)
             mass(i+2*lonp) = mass(i+lonp)
           enddo
c
c prepare grid spacing
c
           do i = lonp+1,2*lonp
             hh(i) = locs(i+1)-locs(i)
           enddo
           do i = 1,lonp
             hh(i) = hh(i+lonp)
             hh (i+2*lonp) = hh(i+lonp)
           enddo
c
c prepare location with monotonic concerns
c
           do i = lonp+1,2*lonp
             dqi = 0.25*(mass(i+1)-mass(i-1))
             dqimax = max(mass(i-1),mass(i),mass(i+1)) - mass(i)
             dqimin = mass(i) - min(mass(i-1),mass(i),mass(i+1))
             dqmono(i) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
           enddo
           do i = 1,lonp
             dqmono(i) = dqmono(i+lonp)
             dqmono(i+2*lonp) = dqmono(i+lonp)
           enddo
c
c compute value at interface with monotone
c
           do i = lonp+1,2*lonp
             qmi(i)=(mass(i-1)*hh(i)+mass(i)*hh(i-1))/(hh(i)+hh(i-1)) 
     &      +(dqmono(i-1)-dqmono(i))/3.0
c    qmi(i)=(mass(i-1)*hh(i)+mass(i)*hh(i-1))/(hh(i)+hh(i-1))
           enddo
           qmi(2*lonp+1)=qmi(lonp+1)
           do i = lonp+1,2*lonp
             qpi(i)=qmi(i+1)
           enddo
c
c do less diffusive
c
           do i = lonp+1,2*lonp
             qmi(i)=mass(i) 
     &    -sign(min(abs(2.*dqmono(i)),abs(qmi(i)-mass(i))),2.*dqmono(i))
           enddo
           do i = lonp+1,2*lonp
             qpi(i)=mass(i)      
     &  +sign(min(abs(2.*dqmono(i)),abs(qpi(i)-mass(i))),2.*dqmono(i))
           enddo
c
c do monotonicity
c
c     if( mono.eq.1 ) then
           do i = lonp+1,2*lonp
             c1=qpi(i)-mass(i)
             c2=mass(i)-qmi(i)
             if( c1*c2.le.0.0 ) then
               qmi(i)=mass(i)
               qpi(i)=mass(i)
             endif
           enddo
           do i = lonp+1,2*lonp
             c1=(qpi(i)-qmi(i))*(mass(i)-0.5*(qpi(i)+qmi(i)))
             c2=(qpi(i)-qmi(i))*(qpi(i)-qmi(i))/6.
             if( c1.gt.c2 ) then
               qmi(i)=3.*mass(i)-2.*qpi(i)
             else if( c1.lt.-c2 ) then
               qpi(i)=3.*mass(i)-2.*qmi(i)
                endif
           enddo
c     endif
c
c extend array with cyclic condition
c
           do i = 1,lonp
             qmi(i)        = qmi(i+lonp)
             qmi(i+2*lonp) = qmi(i+lonp)
             qpi(i)        = qpi(i+lonp)
             qpi(i+2*lonp) = qpi(i+lonp)
           enddo
c
c start interpolation by integral of ppm spline
c
           kkl = kstr
           do i = 1,lonn
             kl = i
             kh = i + 1
c find kkh
             do kk = kkl+1,3*lonp
               if( pn(kh).lt.locs(kk) ) then
                 kkh = kk-1
         go to 100
                endif
             enddo
  100  continue

c mass interpolate
             if( kkh.eq.kkl ) then
               tl=(pn(kl)-locs(kkl))/hh(kkl)
               tl2=tl*tl
               tl3=tl2*tl
               th=(pn(kh)-locs(kkl))/hh(kkl)
               th2=th*th
               th3=th2*th
               qqh=(th3-th2)*qpi(kkl)+
     &          (th3-2.*th2+th)*qmi(kkl)+(-2.*th3+3.*th2)*mass(kkl)
               qql=(tl3-tl2)*qpi(kkl)+                                
     &          (tl3-2.*tl2+tl)*qmi(kkl)+(-2.*tl3+3.*tl2)*mass(kkl)
               qn(i) = (qqh-qql)/(th-tl)
             else if( kkh.gt.kkl ) then
               tl=(pn(kl)-locs(kkl))/hh(kkl)
               tl2=tl*tl
               tl3=tl2*tl
               qql=(tl3-tl2)*qpi(kkl)+                             
     &         (tl3-2.*tl2+tl)*qmi(kkl)+(-2.*tl3+3.*tl2)*mass(kkl)
               dql = mass(kkl)-qql
               th=(pn(kh)-locs(kkh))/hh(kkh)
               th2=th*th
               th3=th2*th
               dqh=(th3-th2)*qpi(kkh)+                              
     &         (th3-2.*th2+th)*qmi(kkh)+(-2.*th3+3.*th2)*mass(kkh)
               dpp  = (1.-tl)*hh(kkl) + th*hh(kkh)
               dpq  = dql*hh(kkl) + dqh*hh(kkh)
               if( kkh-kkl.gt.1 ) then
                 do kk = kkl+1,kkh-1
                   dpp = dpp + hh(kk)
                dpq = dpq + mass(kk)*hh(kk)
                 enddo
               endif
               qn(i) = dpq / dpp
             else
         print *,' Error in cyclic_cell_ppm_intp location not found '
               call abort
             endif
c next one
             kkl = kkh
           enddo
c
           return
           end subroutine cyclic_cell_ppm_intp
c-------------------------------------------------------------------------------
c
c-------------------------------------------------------------------------------
       subroutine vertical_cell_advect(lons,londim,levs,ncld,deltim
     &                                ,ppi,wwi,qql,mass)
c-------------------------------------------------------------------------------
#include <paramodel.h>
#include <comnislq.h>
c-------------------------------------------------------------------------------
           integer         ::  lons,londim,levs,i,k,n,mass,ncld
           real                                  ::  deltim
           real   , dimension(londim,levs+1)     ::  ppi,wwi
           real   , dimension(londim,levs,ncld)  ::  qql
      real   , dimension(lons,levs)         ::  dsfact,rqmm,rqnn,rqda
           real   , dimension(lons,levs+1)       ::  ppii,ppid,ppia
c
           do k = 1,levs+1
             do i = 1,lons
               ppii(i,k)=ppi(i,k)
               ppid(i,k)=ppi(i,k)-wwi(i,k)*deltim
               ppia(i,k)=ppi(i,k)+wwi(i,k)*deltim
             enddo
           enddo
c
           if( mass.eq.1) then
             do k = 1,levs
               do i = 1,lons
         dsfact(i,k)=(ppid(i,k)-ppid(i,k+1))/(ppia(i,k)-ppia(i,k+1))
               enddo
             enddo
           endif
c
           do n = 1,ncld 
             do k = 1,levs
               do i = 1,lons
                 rqmm(i,k) = qql(i,k,n)
               enddo
             enddo
             call vertical_cell_ppm_intp(ppii,rqmm,ppid,rqda,lons,levs)
             if( mass.eq.1 ) then
               do k = 1,levs
                 do i = 1,lons
                   rqda(i,k) = rqda(i,k) * dsfact(i,k)
                 enddo
               enddo
             endif
             call vertical_cell_ppm_intp(ppia,rqda,ppii,rqnn,lons,levs)
             do k = 1,levs
               do i = 1,lons
                 qql(i,k,n)=rqnn(i,k)
                enddo
             enddo
           enddo
c
           return
           end subroutine vertical_cell_advect
c-------------------------------------------------------------------------------
c 
c-------------------------------------------------------------------------------
       subroutine vertical_cell_ppm_intp(pp,qq,pn,qn,lons,levs)
c-------------------------------------------------------------------------------
c
c mass conservation in vertical interpolation: interpolate a group
c of grid point  coordiante call pp at interface with quantity qq at
c cell averaged to a group of new grid point coordinate call pn at
c interface with quantity qn at cell average with ppm spline.
c in vertical with mass conservation is under the condition that
c pp(1)=pn(1), pp(lev+1)=pn(lev+1)
c
c pp    pressure at interfac level as input
c qq    quantity at layer as input
c pn    pressure at interface of new grid structure as input
c qn    quantity at layer as output
c lev  numer of verical layers
c
c author : henry.juang@noaa.gov
c
c-------------------------------------------------------------------------------
#include <paramodel.h>
#include <comnislq.h>
c-------------------------------------------------------------------------------
           integer                      ::  lons,levs            
     &,                                   i,k, kl, kh, kk, kkl, kkh
           real   , dimension(lons,levs)   ::  qq,qn
           real   , dimension(lons,levs+1) ::  pp,pn
c
           real                              ::  dqi,dqimax,dqimin 
     &,                                    tl,tl2,tl3,qql,dql  
     &,                                    th,th2,th3,qqh,dqh
     &,                                    dpp,dpq,c1,c2
           real   , dimension(levs)          ::  hh,dqmono,qmi,qpi
           real   , dimension(0:levs+1)      ::  mass
c
           do i = 1,lons
c
        if( pp(i,1).ne.pn(i,1) .or. pp(i,levs+1).ne.pn(i,levs+1) ) then
        print *,' Error in vertical_cell_ppm_intp for domain values '
        print *,' i,pp1 pn1 ppt pnt ',i,pp(i,1),pn(i,1),pp(i,levs+1)
     &,                                             pn(i,levs+1)
               call abort
             endif
c
c prepare thickness for uniform grid
c
             do k = 1,levs
               hh(k) = pp(i,k+1)-pp(i,k)       
             enddo
c
c prepare location with monotonic concerns
c
             mass(1:levs)=qq(i,1:levs)
c
             mass(0)=(3.*hh(1)+hh(2))*mass(1)-2.*hh(1)*mass(2)
             mass(0)=mass(0)/(hh(1)+hh(2))
             mass(levs+1)=(3.*hh(levs)+hh(levs-1))*mass(levs)-2.*
     &                    hh(levs)*mass(levs-1)
             mass(levs+1)=mass(levs+1)/(hh(levs)+hh(levs-1))
             do k = 1,levs
               dqi = 0.25*(mass(k+1)-mass(k-1))
               dqimax = max(mass(k-1),mass(k),mass(k+1)) - mass(k)
               dqimin = mass(k) - min(mass(k-1),mass(k),mass(k+1))
          dqmono(k) = sign( min( abs(dqimin), dqimin, dqimax ), dqi)
             enddo
c
c compute value at interface with momotone
c
             do k = 2,levs
             qmi(k)=(mass(k-1)*hh(k)+mass(k)*hh(k-1))/(hh(k)+hh(k-1))+
     &                (dqmono(k-1)-dqmono(k))/3.0
             enddo
             do k = 1,levs-1
               qpi(k)=qmi(k+1)
             enddo
             qmi(1)=mass(1)
             qpi(1)=mass(1)
             qmi(levs)=mass(levs)
             qpi(levs)=mass(levs)
c
c do monotonicity
c
c     if( mono.eq.1 ) then
             do k = 1,levs
               c1=qpi(k)-mass(k)
               c2=mass(k)-qmi(k)
               if( c1*c2.le.0.0 ) then
                 qmi(k)=mass(k)
                 qpi(k)=mass(k)
               endif
             enddo
             do k = 1,levs
               c1=(qpi(k)-qmi(k))*(mass(k)-0.5*(qpi(k)+qmi(k)))
               c2=(qpi(k)-qmi(k))*(qpi(k)-qmi(k))/6.
               if( c1.gt.c2 ) then
                 qmi(k)=3.*mass(k)-2.*qpi(k)
               else if( c1.lt.-c2 ) then
                 qpi(k)=3.*mass(k)-2.*qmi(k)
               endif
             enddo
c     endif
c
c start interpolation by integral of ppm spline
c
       kkl = 1
       do k = 1,levs
       kl = k
       kh = k + 1
c find kkh
       do kk = kkl+1,levs+1
c         if( pn(i,kh).ge.pp(i,kk) ) then
         if( pn(i,kh).le.pp(i,kk) ) then    
           kkh = kk-1
           go to 100
         endif
       enddo
c mass interpolate
  100     if( kkh.eq.kkl ) then
           tl=(pn(i,kl)-pp(i,kkl))/hh(kkl)
           tl2=tl*tl
           tl3=tl2*tl
           th=(pn(i,kh)-pp(i,kkl))/hh(kkl)
           th2=th*th
           th3=th2*th
           qqh=(th3-th2)*qpi(kkl)+(th3-2.*th2+th)*qmi(kkl)+  
     &         (-2.*th3+3.*th2)*mass(kkl)
                 qql=(tl3-tl2)*qpi(kkl)+(tl3-2.*tl2+tl)*qmi(kkl)+
     &         (-2.*tl3+3.*tl2)*mass(kkl)
                 qn(i,k) = (qqh-qql)/(th-tl)
               else if( kkh.gt.kkl ) then
                 tl=(pn(i,kl)-pp(i,kkl))/hh(kkl)
                 tl2=tl*tl
                 tl3=tl2*tl
                 qql=(tl3-tl2)*qpi(kkl)+(tl3-2.*tl2+tl)*qmi(kkl)+ 
     &                 (-2.*tl3+3.*tl2)*mass(kkl)
                 dql = qq(i,kkl)-qql
                 th=(pn(i,kh)-pp(i,kkh))/hh(kkh)
                 th2=th*th
                 th3=th2*th
                 dqh=(th3-th2)*qpi(kkh)+(th3-2.*th2+th)*qmi(kkh)+
     &              (-2.*th3+3.*th2)*mass(kkh)
                 dpp= (1.-tl)*hh(kkl) + th*hh(kkh)
                 dpq= dql*hh(kkl) + dqh*hh(kkh)
                 if( kkh-kkl.gt.1 ) then
                   do kk=kkl+1,kkh-1
                     dpp = dpp + hh(kk)
                     dpq = dpq + qq(i,kk)*hh(kk)
                   enddo
                 endif
                 qn(i,k) = dpq / dpp
               else
          print *,' Error in vertical_cell_ppm_intp for no lev found '
                 print *,' i kh kl ',i,kh,kl
                 print *,' pn ',(pn(i,kk),kk=1,levs+1)
                 print *,' pp ',(pp(i,kk),kk=1,levs+1)
                 call abort
               endif
c next one
               kkl = kkh
             enddo
c
1          enddo
c
           return
           end subroutine vertical_cell_ppm_intp
c-------------------------------------------------------------------------------
#endif /* NISLQ end */
