#define TWOJ1 twoj1
#define JCAP jcap
#define JCAP1 jcap1
#define JCAP1P jcap1p
#define JCAPP jcapp
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LCAPP lcapp
#define LCAP22P lcap22p
#define LEVH levh
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LEVSP levsp
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LLN22P lln22p
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LNUV lnuv
#define LONF22P lonf22p
#define LONF22 lonf22
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LPNT lpnt
#define LTSTP ltstp
#define NCLDG ncldg
#define NPES npes
#define NTRAC ntrac
#define NTOTAL ntotal
#define LTAG2S latg2s
#include <define.h>
c#ifdef NISLQ
           subroutine sph_fft_trans
c-------------------------------------------------------------------------------
c
c program history log:
c   1988-04-25  joseph sela
c   2001-01-01  masao kanamitsu        seasonal forecast system
c   2009-10-01  jung-eun kim           f90 format with standard physics modules
c   2010-07-01  myung-seo koo          dimension allocatable with namelist input
c
c-------------------------------------------------------------------------------
#include <paramodel.h>
#ifdef REDUCE_GRID
#include <comreduce.h> 
#define LATS latg_/2
#else
#define LATS 1
#endif
c
#define DEFAULT

#ifdef DEFAULT
      common /comlonf/ifax(20,LATS),trigs(lonf_,2,LATS)               
#endif
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
c                                                                               
           do lat = 1,LATS
#ifdef REDUCE_GRID
             lonff=lonfd(lat)
#else
             lonff=lonf_
#endif
#define DEFAULT
#ifdef DCRFT
#undef DEFAULT
             crscale(lat)=1.0
             rcscale(lat)=1./float(lonff)
#endif
#ifdef FFT_SGIMATH
#undef DEFAULT
             scale(lat)=1./float(lonff)
             call dzfftm1dui (lonff,trigs(1,lat))                      
#endif
#ifdef FFTW
#undef DEFAULT
             ifftw_real_to_complex=-1
             ifftw_complex_to_real=1
             ifftw_estimate=0
             scale(lat)=1./float(lonff)
             call rfftw_f77_create_plan(iplan_c_to_r(lat),lonff
     &,           ifftw_complex_to_real,ifftw_estimate)
             call rfftw_f77_create_plan(iplan_r_to_c(lat),lonff
     &,           ifftw_real_to_complex,ifftw_estimate)
c
c ---- plan can be destryed as following, but may not be necessary ----
c ---  call rfftw_f77_destroy_plan(int_plan)
c
#endif
#ifdef RFFTMLT
#undef DEFAULT
             call fftfax (lonff,ifax(1,lat),trigs(1,1,lat))            
#endif
#ifdef ASLES
#undef DEFAULT
             call ldfrmfb(lonff,0,0,0,0,0,ifax(1,lat)
     &,                   trigs(1,1,lat),0,ierr)
#endif
#ifdef DEFAULT
             call    fax (ifax(1,lat), lonff,3)                       
             call fftrig (trigs(1,1,lat),lonff,3)                    
#endif
#ifdef ASLES
             if (ierr .ge. 3000)  then
               write(6,*)' error in asles fft initialization'
               call abort
             endif
#else
#ifdef DEFAULT
             if (ifax(1,lat) .eq. -99)  then
               write(6,'(A,I2,A)')'error in sph_fft_trans. lonff='
     &,                            lonff,' not factorable.'
               call abort
             endif
#endif
#endif
           enddo
c
           return                         
           end subroutine sph_fft_trans
#ifndef RMP
c-------------------------------------------------------------------------------
           subroutine sph_fft_driver(isgn,grid,wave,ntotall)
c-------------------------------------------------------------------------------
c abstract: Spherical harmonics transform on the sphere
c
c     isgn   =  1 : grid to wave
c     isgn   = -1 : wave to grid
c-------------------------------------------------------------------------------
#include <paramodel.h>
#include <comio.h>
#ifdef MP
#include <commpi.h>
#else
#include <comfgrid.h>
#endif
#ifdef REDUCE_GRID
#include <comreduce.h>
#else
#endif
#include <comgpln.h>

#ifdef MP
#define JCAPS jcapp_
#define JCAP1S jcap1p_
#define LLN2S lln2p_
#define LLN22S lln22p_
#define LEVSS levsp_
#define LCAP22S lcap22p_
#define LONF2S lonf2p_
#define LONF22S lonf22p_
#define LATG2S latg2p_
#define LNT2S lnt2p_
#define LNT22S lnt22p_
#define LEVSS levsp_
#define LEVHS levhp_
#define LCAPS lcapp_
#else
#define JCAPS jcap_
#define JCAP1S jcap1_
#define LLN2S lnt2_
#define LLN22S lnt22_
#define LEVSS levs_
#define LCAP22S lonf22_
#define LONF2S lonf2_
#define LONF22S lonf22_
#define LATG2S latg2_
#define LNT2S lnt2_
#define LNT22S lnt22_
#define LEVSS levs_
#define LEVHS levh_
#define LCAPS lonf_
#endif
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
           integer  ::  isgn,ntotall
           integer  ::  nlevs,nlevsp
           integer  ::  lat, lat1, lat2, latdon       
           integer  ::   lcapf, lonff,llstr,llens,llensd
           integer  ::   lan,i,k,jstr,jend
           real, dimension(LONF2S,levs_*ntotall,LATG2S) ::  grid
           real, dimension(LNT22S,levs_*ntotall)        ::  wave
           real, allocatable, dimension(:,:,:)     ::  anl
           real, allocatable, dimension(:,:,:,:)   ::  flp,flm
#ifdef MP
           real, allocatable, dimension(:,:)       ::  pwave
           real, allocatable, dimension(:,:,:)     ::  anf,grs
#endif
c
c setting variables
c
           nlevs =levs_*ntotall
#ifdef MP
           nlevsp=levsp_*ntotall
#else
           nlevsp=nlevs
#endif
c
c allocation
c
           allocate(anl(LCAP22S ,nlevsp,latg2_)
     &,             flp(2,JCAP1S,nlevsp,latg2_)   
     &,             flm(2,JCAP1S,nlevsp,latg2_)     )
#ifdef MP
           allocate(pwave(lln22p_,nlevsp)           )
           allocate(anf(lonf22_ ,nlevsp ,latg2p_)
     &,                 grs(lonf22p_,nlevs  ,latg2p_)   )
#endif
c
c initialize
c
           anl=0.   ;  flp=0. ;  flm=0.
#ifdef MP
           pwave=0. ;  anf=0. ;  grs=0.
#endif
c
c Latitude bands
c
           lat1=1
           lat2=latg2_
           latdon=0
c
c index for wave and grid
c
#ifdef MP
           jstr=latstr(mype)
           jend=latstr(mype)+latlen(mype)-1
           llstr=lwvstr(mype)
           llens=lwvlen(mype)
#define ANLS anf
#define WAVE pwave
#else
           llstr=0
           llens=jcap1_
#define ANLS anl
#define WAVE wave
#endif
c-------------------------------------------------------------------------------
           if(isgn.eq.1) then   
c-------------------------------------------------------------------------------
             wave=0.
c
             do lat = 1,LATG2S
               do k = 1,nlevs
                 do i = 1,LONF2S
#ifdef MP
                   grs(i,k,lat)=grid(i,k,lat)
#else
                   anl(i,k,lat)=grid(i,k,lat)
#endif
                 enddo
               enddo
             enddo
#ifdef MP
c
c mpnk2nx
c
             call mpnx2nk(grs,lonf22p_,nlevs,anf,lonf22_,nlevsp,latg2p_
     &,                   levs_,levsp_,1,1,ntotall)
c
c Latitude bands for MPI
c
             lat1=jstr
             lat2=jend
             latdon=jstr-1
#endif
c
c grid to wave
c
             do lat = lat1,lat2
               lan=lat-latdon
#ifdef REDUCE_GRID
               lcapf=lcapd(latdef(lat))
               lonff=lonfd(latdef(lat))
#else
               lcapf=jcap1_
               lonff=lonf_
#endif
           call ftilonf(ANLS(1,1,lan),ANLS(1,1,lan),2*nlevsp,lcapf
     &,                           lonff,latdef(lat),-1)
             enddo
#ifdef MP
c
c mpnl2ny
c
             call mpny2nl(ANLS,lonf22_,latg2p_,anl,lcap22p_,latg2_
     &,                   nlevsp,1,nlevsp)
c
c Latitude bands
c
             lat1=1
             lat2=latg2_
             latdon=0
#endif
c
c grid to wave (2)
c
             do lat = lat1,lat2
               lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
               llensd=lcapdp(lat,mype)
#else
               llensd=lcapd(lat)
#endif
#else
               llensd=llens
#endif
               call flpflm(flp(1,1,1,lan),flm(1,1,1,lan)
     &,                           anl(1,1,lan),llensd,nlevsp)
             enddo
c
c sph_comp_coeff2
c
             do lat = lat1,lat2
               lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
               llensd=lcapdp(lat,mype)
#else
               llensd=lcapd(lat)
#endif
#else
               llensd=llens
#endif
               call fl2ip(flp(1,1,1,lan),flm(1,1,1,lan),WAVE
     &,                           qww(1,lat),llstr,llensd,lwvdef,nlevsp)
             enddo
#ifdef MP
c
c 2D decomposition
c
         call mpnk2nn(pwave,lln22p_,levsp_,wave,lnt22p_,levs_,ntotall)
#endif
c-------------------------------------------------------------------------------
           else if(isgn.eq.-1) then
c-------------------------------------------------------------------------------
             grid=0.
#ifdef MP
c
c 2D decomposition
c
        call mpnn2nk(wave,lnt22p_,levs_,pwave,lln22p_,levsp_,ntotall)
#endif
c
c sph_sum_coeff
c
             do lat = lat1,lat2
               lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
                llensd=lcapdp(lat,mype)
#else
               llensd=lcapd(lat)
#endif
#else
               llensd=llens
#endif
               call sums2i(WAVE,anl(1,1,lan),qtt(1,lat),llstr
     &,                           llensd,lwvdef,nlevsp)
             enddo
#ifdef MP
c
c mpnl2ny
c
             call mpnl2ny(anl,lcap22p_,latg2_,anf,lonf22_,latg2p_
     &,                         nlevsp,1,nlevsp)
c
c Latitude bands for MPI
c
             lat1=jstr
             lat2=jend
             latdon=jstr-1
#endif
c
c wave to grid
c
             do lat = lat1,lat2
               lan=lat-latdon
#ifdef REDUCE_GRID
               lcapf=lcapd(latdef(lat))
               lonff=lonfd(latdef(lat))
#else
               lcapf=jcap1_
               lonff=lonf_
#endif
               call ftilonf(ANLS(1,1,lan),ANLS(1,1,lan),2*nlevsp
     &,                         lcapf,lonff,latdef(lat),1)
             enddo



#ifdef MP
c
c mpnk2nx
c
             call mpnk2nx(anf,lonf22_,nlevsp,grs,lonf22p_,nlevs,latg2p_
     &,                  levsp_,levs_,1,1,ntotall)
#endif
             do lat = 1,LATG2S
               do k=1,nlevs
                 do i = 1,LONF2S
#ifdef MP
                   grid(i,k,lat)=grs(i,k,lat)
#else
                   grid(i,k,lat)=anl(i,k,lat)
#endif
                 enddo
               enddo
             enddo
c-------------------------------------------------------------------------------
           endif
c
#undef ANLS
#undef WAVE
           deallocate(anl,flp,flm)
#ifdef MP
           deallocate(pwave,anf,grs)
#endif
           end subroutine sph_fft_driver
c
c-------------------------------------------------------------------------------
           subroutine sph_fft_driver_2d(isgn,grid,wave,ntotall)
c-------------------------------------------------------------------------------
c abstract: Spherical harmonics transform on the sphere
c
c     isgn   =  1 : grid to wave
c     isgn   = -1 : wave to grid
c-------------------------------------------------------------------------------
#include <paramodel.h>
#include <comio.h>
#ifdef MP
#include <commpi.h> 
#else
#include <comfgrid.h>
#endif
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
#include <comgpln.h>
c-------------------------------------------------------------------------------
c-------------------------------------------------------------------------------
           integer ::  isgn,ntotall
           integer ::  lat, lat1, lat2, latdon        
           integer ::  lcapf, lonff,llstr,llens,llensd
           integer ::  lan,i,k,jstr,jend
           real, dimension(LONF2S,LATG2S)          ::  grid
           real, dimension(LNT22S)                 ::  wave
           real, allocatable, dimension(:,:)       ::  anl
           real, allocatable, dimension(:,:,:)     ::  flp,flm
#ifdef MP
           real, allocatable, dimension(:)         ::  pwave
           real, allocatable, dimension(:,:)       ::  anf,grs
#endif

#ifdef MP
#define JCAPS jcapp_
#define JCAP1S jcap1p_
#define LLN2S lln2p_
#define LLN22S lln22p_
#define LEVSS levsp_
#define LCAP22S lcap22p_
#define LONF2S lonf2p_
#define LONF22S lonf22p_
#define LATG2S latg2p_
#define LNT2S lnt2p_
#define LNT22S lnt22p_
#define LEVSS levsp_
#define LEVHS levhp_
#define LCAPS lcapp_
#else
#define JCAPS jcap_
#define JCAP1S jcap1_
#define LLN2S lnt2_
#define LLN22S lnt22_
#define LEVSS levs_
#define LCAP22S lonf22_
#define LONF2S lonf2_
#define LONF22S lonf22_
#define LATG2S latg2_
#define LNT2S lnt2_
#define LNT22S lnt22_
#define LEVSS levs_
#define LEVHS levh_
#define LCAPS lonf_
#endif
c
c allocation
c
           allocate(anl(LCAP22S,latg2_)    
     &,             flp(2,JCAP1S,latg2_)   
     &,             flm(2,JCAP1S,latg2_)     )
#ifdef MP
           allocate(pwave(lln22p_)           )
           allocate(anf(lonf22_ ,latg2p_) 
     &,             grs(lonf22p_,latg2p_)   )
#endif
c
c initialize
c
           anl=0.   ;  flp=0. ;  flm=0.
#ifdef MP
           pwave=0. ;  anf=0. ;  grs=0.
#endif
c
c Latitude bands
c
           lat1=1
           lat2=latg2_
           latdon=0
c
c index for wave and grid
c
#ifdef MP
           jstr=latstr(mype)
           jend=latstr(mype)+latlen(mype)-1
           llstr=lwvstr(mype)
           llens=lwvlen(mype)
#define ANLS anf
#define WAVE pwave
#else
           llstr=0
           llens=jcap1_
#define ANLS anl
#define WAVE wave
#endif
c-------------------------------------------------------------------------------
           if(isgn.eq.1) then   
c-------------------------------------------------------------------------------
             wave=0.
c
             do lat = 1,LATG2S
               do i = 1,LONF2S
#ifdef MP
                 grs(i,lat)=grid(i,lat)
#else
                 anl(i,lat)=grid(i,lat)
#endif
                enddo
             enddo
#ifdef MP
c
c mpnk2x
c
             call mpnx2x(grs,lonf22p_,1,anf,lonf22_,1,latg2p_,1,1,1)
c
c Latitude bands for MPI
c
             lat1=jstr
             lat2=jend
             latdon=jstr-1
#endif
c
c grid to wave
c
             do lat = lat1,lat2
               lan=lat-latdon
#ifdef REDUCE_GRID
               lcapf=lcapd(latdef(lat))
               lonff=lonfd(latdef(lat))
#else
               lcapf=jcap1_
               lonff=lonf_
#endif
               call ftilonf(ANLS(1,lan),ANLS(1,lan),2,lcapf
     &,                           lonff,latdef(lat),-1)
             enddo
#ifdef MP
c
c mpnl2ny
c
             call mpny2nl(ANLS,lonf22_,latg2p_,anl,lcap22p_
     &,                   latg2_,1,1,1)
c
c Latitude bands
c
             lat1=1
             lat2=latg2_
             latdon=0
#endif
c
c grid to wave (2)
c
             do lat = lat1,lat2
               lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
               llensd=lcapdp(lat,mype)
#else
               llensd=lcapd(lat)
#endif
#else
                llensd=llens
#endif
               call flpflm(flp(1,1,lan),flm(1,1,lan)
     &,                 anl(1,lan),llensd,1)
             enddo
c
c sph_comp_coeff2
c
             do lat = lat1,lat2
               lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
               llensd=lcapdp(lat,mype)
#else
               llensd=lcapd(lat)
#endif
#else
               llensd=llens
#endif
               call fl2ip(flp(1,1,lan),flm(1,1,lan),WAVE
     &,                 qww(1,lat),llstr,llensd,lwvdef,1)
             enddo
#ifdef MP
c
c 2D decomposition
c
             call mpn2nn(pwave,lln22p_,wave,lnt22p_,1)
#endif
c-------------------------------------------------------------------------------
           else if(isgn.eq.-1) then   
c-------------------------------------------------------------------------------
             grid=0.
#ifdef MP
c
c 2D decomposition
c
             call mpnn2n(wave,lnt22p_,pwave,lln22p_,1)
#endif
c
c sph_sum_coeff
c
             do lat = lat1,lat2
               lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
               llensd=lcapdp(lat,mype)
#else
                llensd=lcapd(lat)
#endif
#else
               llensd=llens
#endif
               call sums2i(WAVE,anl(1,lan),qtt(1,lat),llstr
     &,                         llensd,lwvdef,1)
             enddo
#ifdef MP
c
c mpnl2ny
c
             call mpnl2ny(anl,lcap22p_,latg2_,anf,lonf22_,latg2p_,1,1,1)
c
c Latitude bands for MPI
c
             lat1=jstr
             lat2=jend
             latdon=jstr-1
#endif
c
c wave to grid
c
             do lat = lat1,lat2
               lan=lat-latdon
#ifdef REDUCE_GRID
               lcapf=lcapd(latdef(lat))
               lonff=lonfd(latdef(lat))
#else
               lcapf=jcap1_
               lonff=lonf_
#endif
               call ftilonf(ANLS(1,lan),ANLS(1,lan),2,lcapf
     &,                         lonff,latdef(lat),1)
             enddo
#ifdef MP
c
c mpnk2nx
c
             call mpx2nx(anf,lonf22_,1,grs,lonf22p_,1,latg2p_,1,1,1)
#endif
             do lat = 1,LATG2S
               do i = 1,LONF2S
#ifdef MP
                 grid(i,lat)=grs(i,lat)
#else
                 grid(i,lat)=anl(i,lat)
#endif
               enddo
             enddo
c-------------------------------------------------------------------------------
           endif
c
#undef ANLS
#undef WAVE
           deallocate(anl,flp,flm)
#ifdef MP
           deallocate(pwave,anf,grs)
#endif
           end subroutine sph_fft_driver_2d
#endif /* ~RMP end */

c#endif /* End of NISLQ */
