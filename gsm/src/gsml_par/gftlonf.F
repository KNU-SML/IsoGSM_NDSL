#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF lonf
#define NPES npes
      subroutine gftlonf
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c
c $Id: gftlonf.F,v 1.2 2009/10/02 00:35:42 kana Exp $
c                                                                               
#include <paramodel.h>
#ifdef REDUCE_GRID
#include <comreduce.h>
#define LATS latg_/2
#else
#define LATS 1
#endif
#define DEFAULT
#ifdef DCRFT
#undef DEFAULT
      common /comlonf/crscale(LATS),rcscale(LATS)
#endif
#ifdef FFT_SGIMATH
#undef DEFAULT
      common /comlonf/scale(LATS),trig(lonf_+15,LATS)                        
#endif
#ifdef FFTW
#undef DEFAULT
      common /comlonf/ scale(LATS)
      common /comlonfi/iplan_c_to_r(LATS),iplan_r_to_c(LATS)
#endif
#ifdef DEFAULT
      common /comlonf/ifax(20,LATS),trigs(lonf_,2,LATS)                        
#endif
c                                                                               
      do lat=1,LATS
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
       call rfftw_f77_create_plan(iplan_c_to_r(lat),lonff,
     +                            ifftw_complex_to_real,ifftw_estimate)
       call rfftw_f77_create_plan(iplan_r_to_c(lat),lonff,
     +                            ifftw_real_to_complex,ifftw_estimate)
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
c            print *,'new fft initializaiton starts'
            call ldfrmfb(lonff,0,0,0,0,0,ifax(1,lat),
     1                                      trigs(1,1,lat),0,ierr)
#endif
#ifdef DEFAULT
       call    fax (ifax(1,lat), lonff,3)                                         
       call fftrig (trigs(1,1,lat),lonff,3)                                         
#endif
#ifdef ASLES
            if (ierr .ge. 3000)  print 150
            if (ierr .ge. 3000)  stop
150         format (' error in asles fft initialization')
#else
#ifdef DEFAULT
            if (ifax(1,lat) .eq. -99)  then
              print 120 ,lonff                                   
#ifdef MP
#ifdef RSM
                call rmpabort
#else
                call mpabort
#endif
#else
                call abort
#endif
            endif
120    format (' error in gftlonf. lonff=',i12,' not factorable. ')
#endif
#endif
      enddo
c
      return                                                                    
      end                                                                       
