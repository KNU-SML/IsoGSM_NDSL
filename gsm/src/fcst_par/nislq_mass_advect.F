#define JCAP jcap
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LNUV lnuv
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NCLDG ncldg
#define NPES npes
#define NROW nrow
#define NTRAC ntrac
#define NTOTAL ntotal
#define LEVSP levsp
#define G g
#define PI pi
#define RRERTH rrerth
#include <define.h>
#ifdef NISLQ_MASS
        subroutine nislq_mass_advect(deltimm,pt,ut,vt,pdot,q1,q3)
c-------------------------------------------------------------------------------
c
c a routine to do non-iteration semi-Lagrangain finite volume advection
c considering mass advection together with gases
c contact: hann-ming henry juang
c
c <in : top to bottom>
c deltimm  time step from n to n+1
c pt      surface pressure in cb
c ut      horizontal u wind scaled without earth radius cos(phi)^2d(lamda)/dt
c vt      horizontal v wind scaled without earth radius cos(phi)d(phi)/dt
c pdot    vertical wind in dp/dt in cb
c q1      n-1 time step q
c
c <out: bottom to top>
c q3      n+1 time step q (bottom to top)
c
c-------------------------------------------------------------------------------

#include <paramodel.h>
#include <constant.h>
#include <comfgrid.h>
#include <comfspec.vr.h>
#include <comfver.h>
#include <comio.h>
#ifdef MP
#include <commpi.h>
#endif
#include <comreduce.h>
#include <comnislq.h>
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c
#ifdef MP
#define LEVSS levsp_
#else
#define LEVSS levs_
#endif
c passing variables
c
           real                                       ::  deltimm
           real,    dimension(LONF2S,        LATG2S)   ::  pt
           real,    dimension(LONF2S,levs_  ,LATG2S)   ::  ut,vt
           real,    dimension(LONF2S,levs_+1,LATG2S)   ::  pdot
           real,    dimension(LONF2S,levh_  ,LATG2S)   ::  q1
           real,    dimension(LONF2S,levh_  ,LATG2S)   ::  q3
c
c local variables
c
           integer            , parameter                ::  mass=1
           integer :: ncld
#ifndef MP
           integer            , parameter                ::  nsize=1
#else
           integer            , parameter                ::  nsize=nrow_
#endif
           integer ::  jjend,j,j1,j2 
           integer ::  jj,lat        
           integer ::  i,lonsd       
           integer ::  kqp,k,kk
c dp2   
           real  , dimension(LONF2S ,levs_         )  ::  dp2
           real  , dimension(LONF2S ,levs_+1       )  ::  ppi,pdot2
           real  , dimension(LONF2S ,levs_ ,LATG2S )  ::  dp1,dp3
           real  , dimension(LONF2S ,nlevs ,LATG2S )  ::  qp
           real  , dimension(lonf2_ ,nlevsp,LATG2S )  ::  qt
#ifdef MP
           real  , dimension(lonf2_ ,levsp_,latg2p_)  ::  utp,vtp
           real  , dimension(lonf2_ ,nlevsp,latg2p_)  ::  qpp
           real  , dimension(LONF2S ,nlevs ,LATG2S )  ::  qtp
#endif
           real  , dimension(LONF2S ,nlevs         )  ::  qtn
           real  , dimension(lonfull,LEVSS ,latpart)  ::  uulon,vvlon
           real  , dimension(latfull,LEVSS ,lonpart)  ::  vvlat
           real  , dimension(lonfull,nlevsp,latpart)  ::  qqlon,rrlon
           real  , dimension(latfull,nlevsp,lonpart)  ::  qqlat,rrlat

           ncld = ntotal_
c
c initialize
c
        q3=0.
c
c dp2   dp2=0.
c
           dp1=0.;  dp3=0.
           qp=0. ;  qt=0. ;  ppi=0.
#ifdef MP
           utp=0.;  vtp=0.;  qpp=0.;  qtp=0.
#endif
           uulon=0. ;  vvlon=0. ;  vvlat=0.
           qqlon=0. ;  rrlon=0. ;  qqlat=0.  ; rrlat=0.
c
c k-index for q*dp
c
c sldp   kqp=levs_
c
           kqp=0
c
c latitude band
c
#ifdef MP
           jjend=latlen(mype)
#else
           jjend=latg2_
#endif
c
c density (by spectral dynamics)
c
c dp at time step n-1
           call nislq_dp(qm,dp1) 
c dp at time step n+1
           call nislq_dp(z ,dp3)
           do j=1,jjend
#ifdef REDUCE_GRID
#ifdef MP
             lonsd=lonfdp(j,mype)*2
#else
             lonsd=lonfd(latdef(j))*2
#endif
#else
             lonsd=LONF2S
#endif /* REDUCE_GRID end */
#ifdef SLDBG
             if( iope ) print *,'j,lonsd in nislq_advect',j,lonsd
#endif
c
c sldp     c dry air (dp) for SL advection
c sldp     forall(i=1:lonsd,k=1:levs_) qp(i,k,j)=dp1(i,k,j)
c
c moisture (q*dp) for SL advection
c
        do k = 1,levh_
               kk=mod(k-1,levs_)+1
               do i = 1,lonsd
                 qp(i,kqp+k,j)=q1(i,k,j)*dp1(i,kk,j)
               enddo
             enddo
           enddo
#ifdef MP
c
c transpose z-full to x-full
c
           call mpnx2nk(ut,lonf2p_,levs_,utp,lonf2_,levsp_,latg2p_,levs_
     &,                 levsp_,1,1,1)
           call mpnx2nk(vt,lonf2p_,levs_,vtp,lonf2_,levsp_,latg2p_,levs_
     &,                 levsp_,1,1,1)
           call mpnx2nk(qp,lonf2p_,nlevs,qpp,lonf2_,nlevsp,latg2p_,levs_
     &,                 levsp_,1,1,ncld)
#define UU utp
#define VV vtp
#define QQ qpp
#else
#define UU ut
#define VV vt
#define QQ qp
#endif /* MP end */
#ifdef SLDBG
           if( iope ) then
             print *,' enter nislq_advect  with positive definition '
             call print_maxmin_six(qp  ,LONF2S*nlevs
     &,                            LATG2S ,1,LATG2S,'qp   input')
             call print_maxmin_six(ut  ,LONF2S*levs_
     &,                            LATG2S ,1,LATG2S,'ut   input')
             call print_maxmin_six(pdot,LONF2S*(levs_+1),LATG2S 
     &,                            1,LATG2S,'pdot input')
           endif
#endif
c
c first mass conserving interpolation from reduced grid to full grid
c
           do jj = 1,jjend
             j1=2*jj-1    
             j2=2*jj     
#ifdef REDUCE_GRID
#ifdef MP
      lat=latdef(latstr(mype)+jj-1)
#else
      lat=latdef(jj)
#endif
      lonsd=lonfd(lat)
#else
      lat=latdef(jj)
      lonsd=nx
#endif /* REDUCE_GRID end */
c u,v at time step n (convert dynamics to SL grid)
      do k = 1,LEVSS
       do i = 1,lonsd
         uulon(i,k,j1) = UU(i      ,k,jj) * (rrerth_*rbs2(jj))
         uulon(i,k,j2) = UU(lonsd+i,k,jj) * (rrerth_*rbs2(jj))
         vvlon(i,k,j1) = VV(i      ,k,jj) * (rrerth_*sqrt(rbs2(jj)))
         vvlon(i,k,j2) = VV(lonsd+i,k,jj) * (rrerth_*sqrt(rbs2(jj)))
       enddo
      enddo
c dp, q*dp at time step n-1 (convert dynamics to SL grid)
      do k = 1,nlevsp
        do i = 1,lonsd
         qqlon(i,k,j1) = QQ(i      ,k,jj)/sqrt(rbs2(jj))
         qqlon(i,k,j2) = QQ(lonsd+i,k,jj)/sqrt(rbs2(jj))
        enddo
      enddo
#ifdef SLDBG
      if( iope ) then
       print *,'jj,lat,lonsd in nislq_advect',jj,lat,lonsd
       if(lat.eq.7) then
         call print_maxmin_six(uulon(1,1,j1),nx,LEVSS ,1,LEVSS
     &,                        'u in dlamda/dt')
         call print_maxmin_six(vvlon(1,1,j1),nx,LEVSS ,1,LEVSS 
     &,                        'v in dphi/dt  ')
         call print_maxmin_six(qqlon(1,1,j1),nx,nlevsp,1,nlevsp   
     &,                        'q in dq/dt    ')
       endif
      endif
#endif
#undef UU
#undef VV
#undef QQ
#ifdef REDUCE_GRID
c
c reduced grid to full grid
c
        call cyclic_cell_intpx(LEVSS ,lonsd,lonfull,uulon(1,1,j1))
        call cyclic_cell_intpx(LEVSS ,lonsd,lonfull,vvlon(1,1,j1))
        call cyclic_cell_intpx(nlevsp,lonsd,lonfull,qqlon(1,1,j1))
c
        call cyclic_cell_intpx(LEVSS ,lonsd,lonfull,uulon(1,1,j2))
        call cyclic_cell_intpx(LEVSS ,lonsd,lonfull,vvlon(1,1,j2))
        call cyclic_cell_intpx(nlevsp,lonsd,lonfull,qqlon(1,1,j2))
#endif
c
        do k = 1,nlevsp
          do i = 1,lonfull
            rrlon(i,k,j1) = qqlon(i,k,j1)
            rrlon(i,k,j2) = qqlon(i,k,j2)
          enddo
        enddo
c
c first set positive advection in horziontal direction with mass conserving

        call cyclic_cell_massadvx(LEVSS,ncld,deltimm 
     &,                           uulon(1,1,j1),rrlon(1,1,j1),mass)
        call cyclic_cell_massadvx(LEVSS,ncld,deltimm 
     &,                           uulon(1,1,j2),rrlon(1,1,j2),mass)
        enddo
#ifdef SLDBG
           if( iope ) then
             print *,' finish cyclic_cell_massadvx'
             call print_maxmin_six(rrlon,lonfull*nlevsp
     &,                            latpart,1,latpart,'r advx 1st     ')
             call print_maxmin_six(vvlon,lonfull*LEVSS
     &,                            latpart,1,latpart,'vv before we2ns')
           endif
#endif
c ---------------------------------------------------------------------
c mpi para from horizontal full grid to meridional full grid
c ---------------------------------------------------------------------
c
c para vvlon, qqlon, and rrlon to vvlat, qqlat, rrlat
c
           call nislq_transpose_we2ns(vvlon,vvlat,LEVSS ,nsize)
           call nislq_transpose_we2ns(qqlon,qqlat,nlevsp,nsize)
           call nislq_transpose_we2ns(rrlon,rrlat,nlevsp,nsize)
#ifdef SLDBG
           if( iope ) then
             print *,' nislq transport from we to ns '
             call print_maxmin_six(vvlat,latfull*LEVSS
     &,                            mylonlen,1,mylonlen,'v we2ns')
             call print_maxmin_six(qqlat,latfull*nlevsp
     &,                            mylonlen,1,mylonlen,'q we2ns')
             call print_maxmin_six(rrlat,latfull*nlevsp
     &,                            mylonlen,1,mylonlen,'r we2ns')
           endif
#endif
c ---------------------------------------------------------------------
c ------------------- in meridional great circle ----------------------
c ---------------------------------------------------------------------
#ifdef SLDBG
           if ( iope ) then
             print *,' nislq adv loop in y '
             print *,' mylonlen=',mylonlen
           endif
#endif
           do i = 1,mysllonlen
c 
c first set advection in meridional direction in great circle through two poles
c
             call cyclic_cell_massadvy(LEVSS,ncld,deltimm
     &,                                vvlat(1,1,i),rrlat(1,1,i),mass)
c
c second set advection in meridional direction in great circle through two poles
c
             call cyclic_cell_massadvy(LEVSS,ncld,deltimm
     &,                                vvlat(1,1,i),qqlat(1,1,i),mass)

           enddo
c
#ifdef SLDBG
           if( iope ) then
             call print_maxmin_six(rrlat,latfull*nlevsp
     &,                            mylonlen,1,mylonlen,'r advy 1st')
             call print_maxmin_six(qqlat,latfull*nlevsp
     &,                            mylonlen,1,mylonlen,'q advy 2nd')
           endif
#endif
c ----------------------------------------------------------------------
c mpi para from meridional direction to horizontal directory 
c ----------------------------------------------------------------------
c
c para qqlat and rrlat to qqlon and rrlon
c
           call nislq_transpose_ns2we(qqlat,qqlon,nlevsp,nsize)
           call nislq_transpose_ns2we(rrlat,rrlon,nlevsp,nsize)
#ifdef SLDBG
           if( iope ) then
             print *,' nislfv_advq transport from ns to we '
             call print_maxmin_six(qqlon,lonfull*nlevsp
     &,                            latpart,1,latpart,'q ns2we 2nd')
             call print_maxmin_six(rrlon,lonfull*nlevsp
     &,                            latpart,1,latpart,'r ns2we 1st')
           endif
#endif
c ---------------------------------------------------------------
c ---------------- back to east-west direction ------------------
c ---------------------------------------------------------------
c     print *,' nislq adv loop in x for last '
c
           do jj = 1,jjend
             j1=2*jj-1    
             j2=2*jj     
#ifdef REDUCE_GRID
#ifdef MP
             lat=latdef(latstr(mype)+jj-1)
#else
             lat=latdef(jj)
#endif
             lonsd=lonfd(lat)
#else
             lat=latdef(jj)
             lonsd=nx
#endif
c
c second set advection in x for the second of the pair
c
             call cyclic_cell_massadvx(LEVSS,ncld,deltimm 
     &,                                uulon(1,1,j1),qqlon(1,1,j1),mass)
             call cyclic_cell_massadvx(LEVSS,ncld,deltimm
     &,                                uulon(1,1,j2),qqlon(1,1,j2),mass)
             do k = 1,nlevsp
               do i = 1,lonsd
                 rrlon(i,k,j1) = 0.5 * ( qqlon(i,k,j1) + rrlon(i,k,j1) )
                 rrlon(i,k,j2) = 0.5 * ( qqlon(i,k,j2) + rrlon(i,k,j2) )
               enddo
             enddo
#ifdef REDUCE_GRID
c
c full gird to reduced grid
c
             call cyclic_cell_intpx(nlevsp,lonfull,lonsd,rrlon(1,1,j1))
             call cyclic_cell_intpx(nlevsp,lonfull,lonsd,rrlon(1,1,j2))
#endif
c
c convert SL to dynamics grid
c
             do k = 1,nlevsp
               do i = 1,lonsd
                 qt(i      ,k,jj)=rrlon(i,k,j1)*sqrt(rbs2(jj))
                 qt(lonsd+i,k,jj)=rrlon(i,k,j2)*sqrt(rbs2(jj))
               enddo
             enddo
           enddo
#ifdef MP
c
c transpose x-full to z-full
c
           call mpnk2nx(qt,lonf2_,nlevsp,qtp,lonf2p_,nlevs
     &,                 latg2p_,levsp_,levs_,1,1,ncld)
#define QT qtp
#else
#define QT qt
#endif /* MP end */
c --------------------------------------------------------------
c ----------- compute vertical advection and total ------------
c --------------------------------------------------------------
           do j = 1,jjend
#ifdef REDUCE_GRID
#ifdef MP
             lonsd=lonfdp(j,mype)*2
#else
             lonsd=lonfd(latdef(j))*2
#endif
#else
             lonsd=LONF2S
#endif /* REDUCE_GRID end */
c
c    pressure (top to bottom)
c
             do k = 1,levs_+1
               do i = 1,lonsd
#ifdef HYBRID
                 ppi(i,k)=ak5(k)+bk5(k)*pt(i,j)
#else
                 ppi(i,k)=si(levs_+2-k)*pt(i,j)
#endif /* HYBRID end */
                 pdot2(i,k)=pdot(i,k,j)
               enddo
             enddo
c dp2     c dp at time step n
c dp2     forall(i=1:lonsd,k=1:levs_) dp2(i,k)=ppi(i,k+1)-ppi(i,k)
c q*dp weighted by dp(n) before vertical advection
      do k = 1,nlevs
       kk=mod(k-1,levs_)+1
       do i = 1,lonsd
         qtn(i,k)=QT(i,k,j)/dp3(i,kk,j)
cdp2         qtn(i,k)=QT(i,k,j)/dp2(i,kk)
       enddo
      enddo
#undef QT
#ifdef SLDBG
      if( iope ) then
       if(j.eq.1) call print_maxmin_seven(ppi,lonsd,LONF2S
     &,           lev+1,1,lev+1,'ppi')
      endif
#endif
c
c vertical advection with mass conserving positive advection
c
       call vertical_cell_advect(lonsd,LONF2S,lev,ncld,deltimm
     &,                          ppi,pdot2,qtn,mass)
c
c    q update at time step n+1 (bottom to top)
c
      do k = 1,levh_
       kk=mod(k-1,levs_)+1
       do i = 1,lonsd
         q3(i,levh_+1-k,j)=qtn(i,kqp+k)
csldp         q3(i,levh_+1-k,j)=qtn(i,kqp+k)/qtn(i,kk)
       enddo
      enddo
           enddo
c
#ifdef SLDBG
           if( iope ) then
             call print_maxmin_six(qqlon,lonfull*nlevsp
     &,                            latpart,1,latpart,'q advx 2nd  ')
             call print_maxmin_six(rrlon,lonfull*nlevsp
     &,                            latpart,1,latpart,'r after mean')
             call print_maxmin_six(q1,LONF2S*levh_ 
     &,                            LATG2S ,1,LATG2S ,'q1 input    ')
             call print_maxmin_six(q3,LONF2S*levh_
     &,                            LATG2S ,1,LATG2S ,'q3 output   ')
           endif
#endif
c
           return
           end subroutine nislq_mass_advect
#endif /* NISLQ_MASS end */
c
#ifdef NISLQ_MASS
c-------------------------------------------------------------------------------

#include <define.h>
      subroutine nislq_dp(qmmm,dp)
c-------------------------------------------------------------------------------
#include <paramodel.h>
#include <comfver.h>
#ifdef MP
#include <commpi.h>
#endif

#ifdef REDUCE_GRID
#include <comreduce.h>
#ifdef MP
#else
#include <comfgrid.h>
#endif

#endif /* REDUCE_GRID end */
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c
c passing variables
c
       real   , dimension(LNT22S)                 ::  qmmm
       real   , dimension(LONF2S,levs_,LATG2S)    ::  dp
c
c local variables
c
       integer    ::  i,j,k,lonsd
       real       , dimension(LONF2S,LATG2S)          ::  ps
       real       , dimension(LONF2S,levs_+1,LATG2S)  ::  prs
c
       call sph_fft_driver_2d(-1,ps,qmmm,1)
c
c latitude band
c
#ifdef MP
       do j = 1,latlen(mype)
#else
       do j = 1,latg2_
#endif
c grid number for longitude
#ifdef REDUCE_GRID
#ifdef MP
      lonsd=lonfdp(j,mype)*2
#else
      lonsd=lonfd(latdef(j))*2
#endif
#else
      lonsd=LONF2S
#endif /* REDUCE_GRID end */

c pressure
      do k = 1,levs_+1
       do i = 1,lonsd
#ifdef HYBRID
         prs(i,k,j)=ak5(k)+bk5(k)*exp(ps(i,j))
#else /* SIGMA */
         prs(i,k,j)=si(levs_+2-k)*exp(ps(i,j))
#endif /* HYBRID end */
       enddo
      enddo
c dp
      do k = 1,levs_
       do i = 1,lonsd
         dp(i,k,j)=prs(i,k+1,j)-prs(i,k,j)
       enddo
      enddo
      enddo 
c
           return
           end subroutine nislq_dp

#endif 
