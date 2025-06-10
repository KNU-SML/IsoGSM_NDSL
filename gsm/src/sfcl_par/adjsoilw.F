#define LATG latg
#define LONF lonf
#define LSOIL lsoil
      subroutine adjsoilw(idim,jdim,lsoil,
     1                    intobsp,year,month,day,hour,fhour,
     2                    obspcp,mdlpcp,land,sfctmp,snow,
     3                    soilm,mdlrnf)
c
      implicit none
c
c	take observed precip - analyzed precip
c	and use it to correct surface moisture
c
c	assume that bges and observed precip are correct
c v1.0  based on Wesley's adjbges program
c
#include <paramodel.h>
c
c   input:
c     idim       : x-dimension
c     jdim       : y-dimension
c     lsoil      : number of soil layers
c     intobsp    : interval of obs precip in hour unit
c     year, month, day, hour :
c     fhour      : forecast hour
c     obspcp     : observed accumulated precipitation in m
c     mdlpcp     : model accumulated precipitation in m
c     land       : land sea mask
c     sfctmp     : surface temperature
c     soilm      : volmetric soil moisture
c     mdlrnf     : accumulated runoff m*sec  (?   RUNOF=Runoff [kg/m**2] )
c     snow       : water equivalent snow depth in  mm
c
c   output:
c     soilm      : volmetric soil moisture
c     mdlrnf     : runoff m*sec
c     snow       : water equivalent snow depth  mm
c
      integer idim, jdim, lsoil
      integer intobsp
      real fhour                                   ! forecast hour
      real sfctmp(idim,jdim)                       ! sfc temp
      real soilm(idim,jdim,lsoil)                  ! soil moisture
      real snow(idim,jdim)                         ! sfc temp
      real temp(idim,jdim), tempv(idim,jdim,lsoil) ! temp array
c
      real jakhak(idim,jdim)
c
      real obspcp(idim,jdim) ! observed precip
      real mdlpcp(idim,jdim) ! model precip
      real land(idim,jdim)   ! land mask, 1=land 0=water
      real mdlrnf(idim,jdim) ! runoff
      real offend(idim,jdim) ! offset residual 
c
      real soilcap(lsoil_)   ! soil moisture capacity for each layer (cm)
      integer year, month, day, hour
c
      integer i,j,k
      real offset, smoist
      real rmax,rmin
c
      real smcmax, smcmin
c
#ifdef DBG 
      real soilmprev(idim,jdim,lsoil_)
      real snowprev (idim,jdim)
#endif
      real ssnow
c
#ifdef OSULSM1
      data soilcap/0.1,1.9/
#endif
#ifdef OSULSM2
      data soilcap/0.1,1.9/
#endif
#ifdef NOALSM1
      data soilcap/0.1,0.3,0.6,1.0/   ! thicknesses
#endif
c
c
      data smcmin,smcmax/0.029,0.4698/
c
c
#ifdef DBG
      write(*,*) 'obspcp:',obspcp(1,1), obspcp(100,1)
      write(*,*) 'mdlpcp:',mdlpcp(1,1), mdlpcp(100,1)
      write(*,*) 'mdlrnf:',mdlrnf(1,1), mdlrnf(100,1)
      write(*,*) 'snow  :',snow  (1,1), snow  (100,1)
      call stats(obspcp,idim,jdim,'opcp')
      call stats(mdlpcp,idim,jdim,'mpcp')
      call stats(mdlrnf,idim,jdim,'mrnf')
#endif
c
c  accumulated precipitation = meter
c  note that runoff is also in meter unit multiplied by time step (mm I think, perhaps notation applies to osu/Rean-2 ?)
c
#ifdef DBG
      do k = 1, lsoil
        rmin = soilm(1,1,k)
        rmax = soilm(1,1,k)
        do j = 1, jdim
          do i = 1, idim
             rmin = min(rmin, soilm(i,j,k))
             rmax = max(rmax, soilm(i,j,k))
          enddo
        enddo
        write(*,*) 'old level=',k,' max/min=',rmax,rmin
      enddo
      do k=1,lsoil
        do j=1,jdim
          do i=1,idim
            soilmprev(i,j,k)=soilm(i,j,k)
          enddo
        enddo
      enddo
      do j=1,jdim
        do i=1,idim
          snowprev(i,j)=snow(i,j)
        enddo
      enddo
#endif
c
c	now to do the correction
c
c	v1.0
c
c	   obspcp .ge. mdlpcp
c		=> assume extra precip would have gone into first soil layer
c		=> if first soil layer is saturated, go to runoff
c  	=> we leave the model to runoff excess soil moisture beyond the
c      soil capacity (which is a function of geographical location).
c	   surface temp < 0c .. add snow and no change to soil moisture
c
c    obspcp .le. mdlpcp
c   => take away equivalent water from first soil layer
c   => if soil moisture at first layer is zero, take away from second layer
c      and so on.
c   => if soil is frozen, take away from snow.
c      if no snow, do nothing.
c
      do j = 1, jdim
        do i = 1, idim
c
c  check for water point c if (land(i,j).eq.0) goto 200
c
          if(land(i,j).eq.0.) go to 200
*
c
c  figure out the precip difference
c
          offset = obspcp(i,j) - mdlpcp(i,j)
c
c  check for surface temperature
c
          if (sfctmp(i,j).ge.273.16) then
c
c  figure out the precip surplus
c
c    this is not precipitation rate, but accumulated precipitation 
c    in meter unit
c
c
             if (offset.gt.0.0.and.soilm(i,j,1).ge.smcmax) then
c
c intobsp=24 runoff ~.03 * 86400 = 2592 * .15 = 389mm?(/1000.0)
c
               mdlrnf(i,j)=mdlrnf(i,j)+offset*float(intobsp)*3600.
               offset = 0.0
             endif
c
             if (offset.eq.0.0) goto 200
c
             do k = 1, lsoil
               smoist = soilm(i,j,k)*soilcap(k) + offset 
               if (smoist.lt.smcmin*soilcap(k)) then
                  offset = smoist - (smcmin*soilcap(k))
                  soilm(i,j,k) = smcmin
               else if (smoist.gt.(soilcap(k)*smcmax)) then
                  offset = smoist - (soilcap(k)*smcmax)
                  soilm(i,j,k) = smcmax
               else
                  offset = 0.0
                  soilm(i,j,k) = smoist / soilcap(k)
                 goto 200
               endif
             enddo
          else
             ssnow = snow(i,j) + offset*1000.
             if(ssnow.lt.0.) then
                offset= ssnow/1000.0
                snow(i,j)=0.
             endif
          endif
  200     continue
          offend(i,j)=offset
        enddo
      enddo
c
#ifdef DBG
      do k = 1, lsoil
        rmin = soilm(1,1,k)
        rmax = soilm(1,1,k)
        do j = 1, jdim
          do i = 1, idim
            rmin = min(rmin, soilm(i,j,k))
            rmax = max(rmax, soilm(i,j,k))
          enddo
        enddo
        write(*,*) 'new level=',k,' max/min=',rmax,rmin
      enddo
      do k=1,lsoil
        do j = 1, jdim
          do i = 1, idim
            soilmprev(i,j,k)=soilm(i,j,k)-soilmprev(i,j,k)
          enddo
        enddo
      enddo
      do j = 1, jdim
        do i = 1, idim
          snowprev(i,j)=snow(i,j)-snowprev(i,j)
        enddo
      enddo
      do k = 1, lsoil
        rmin = soilmprev(1,1,k)
        rmax = soilmprev(1,1,k)
        do j = 1, jdim
          do i = 1, idim
            rmin = min(rmin, soilmprev(i,j,k))
            rmax = max(rmax, soilmprev(i,j,k))
          enddo
        enddo
        write(*,*) 'soilm correction level=',k,' max/min=',rmax,rmin
      enddo
      rmin = snowprev(1,1)
      rmax = snowprev(1,1)
      do j = 1, jdim
        do i = 1, idim
          rmin = min(rmin, snowprev(i,j))
          rmax = max(rmax, snowprev(i,j))
        enddo
      enddo
      write(*,*) 'snow correction max/min=',rmax,rmin
#endif
      return
      end
