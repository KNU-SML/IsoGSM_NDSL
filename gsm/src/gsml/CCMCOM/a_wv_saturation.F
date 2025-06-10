      module wv_saturation      

        private
        save
 
        public aqsat ! Returns saturation vapor pressure
 
      contains

      subroutine aqsat(t       ,p       ,es      ,qs        ,ii      ,
     1                 ilen    ,kk      ,kstart  ,kend      )
!------------------------------------------------------------
#define RD rd
#define RV rv
#ifdef CRAY_THREAD
cfpp$ noconcur r
cfpp$ expand(fpvs,fpvs0)
#endif
                                                                                
#include <paramodel.h>
                                                                                
        real,parameter :: eps=rd_/rv_
        real,parameter :: epsm1=eps-1.
        real,parameter :: qwmx=1.e4
        real,parameter :: epsq=2.e-12
!------------------------------Arguments--------------------------------
!
! Input arguments
!
        integer, intent(in) :: ii       ! I dimension of arrays t, p, es, qs
        integer, intent(in) :: kk       ! K dimension of arrays t, p, es, qs
        integer, intent(in) :: ilen     ! Length of vectors in I direction which
        integer, intent(in) :: kstart   ! Starting location in K direction
        integer, intent(in) :: kend     ! Ending location in K direction
        real,    intent(in) :: t(ii,kk) ! Temperature
        real,    intent(in) :: p(ii,kk) ! Pressure
!
! Output arguments
!
        real, intent(out) :: es(ii,kk)  ! Saturation vapor pressure
        real, intent(out) :: qs(ii,kk)  ! Saturation specific humidity
!
c
c local
c
        real pres ! [cb]
        real qw
        integer i,k


        do k = kstart, kend
        do i = 1, ilen

          pres = p(i,k) * 1.e-3   ! Pa -> cb

#ifdef ICE
          qw      = fpvs(t(i,k))
#else
          qw      = fpvs0(t(i,k))
#endif
          es(i,k) = qw * 1.e3     ! cb -> Pa

          qw      = eps * qw / (pres + epsm1 * qw)
          if (qw.lt.0.0)  qw = qwmx
          qs(i,k) = max(qw,epsq)

        end do
        end do

        return
      end subroutine aqsat

      end module wv_saturation      
