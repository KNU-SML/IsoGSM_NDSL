#define IGRD igrd
#define IGRD12 igrd12
#define JGRD jgrd
#define JGRD12 jgrd12
#define PI pi
      subroutine gulfofca(rdata,slmask)
#include <paramodel.h>
#ifdef RSM
#define ILOTS igrd12_
#define JLOTS jgrd12_
      parameter(idim=ILOTS/2,jdim=JLOTS*2)
c
c add gulf of California sst by tsfcinc
c
      dimension rdata(idim,jdim),slmask(idim,jdim)
      dimension rlat(idim,jdim),rlon(idim,jdim)
c
      common /comrloc/
     1 projo,trutho,oriento,delxo,delyo,cenlato,cenlono,rftgrdo,btmgrdo
     2,cproj,ctruth,corient,cdelx,cdely,ccenlat,ccenlon,clftgrd,cbtmgrd
c
c     parameter (xlon=-108.,xlat=23.)
      parameter (tsfcinc=2.)
c
c  compute latitude and longitude
c
      call setgrd(rlat,rlon,delxo,delyo,dlamda0)
      call shalfo(rlat,1)
      call shalfo(rlon,1)
c
      io=idim
      jo=jdim
c
      do j=1,jo
        do i=1,io
          rlat(i,j) = rlat(i,j) * 180. / pi_
          rlon(i,j) = rlon(i,j) * 180. / pi_
          if(rlon(i,j).gt.180.) then
            rlon(i,j)=rlon(i,j)-360.
          endif
        enddo
      enddo
c
c     print *,'rlat(1,1)=',rlat(1,1)
c     print *,'rlon(1,1)=',rlon(1,1)
c
      do j=1,jdim
         do i=1,jdim
           if(slmask(i,j).eq.0.) then
            if(rlat(i,j).ge.29.4.and.rlat(i,j).lt.32.) then
               if(rlon(i,j).gt.-115.0.and.rlon(i,j).lt.-105) then
                  rdata(i,j)=rdata(i,j)+tsfcinc
                  print *,'lat=',rlat(i,j),' lon=',rlon(i,j)
               endif
            endif
            if(rlat(i,j).ge.27.0.and.rlat(i,j).lt.29.4) then
               if(rlon(i,j).gt.-114.0.and.rlon(i,j).lt.-105) then
                  rdata(i,j)=rdata(i,j)+tsfcinc
                  print *,'lat=',rlat(i,j),' lon=',rlon(i,j)
               endif
            endif
            if(rlat(i,j).ge.24.7.and.rlat(i,j).lt.27.0) then
               if(rlon(i,j).gt.-112.0.and.rlon(i,j).lt.-105) then
                  rdata(i,j)=rdata(i,j)+tsfcinc
                  print *,'lat=',rlat(i,j),' lon=',rlon(i,j)
               endif
            endif
            if(rlat(i,j).ge.24.0.and.rlat(i,j).lt.24.7) then
               if(rlon(i,j).gt.-111.0.and.rlon(i,j).lt.-105) then
                  rdata(i,j)=rdata(i,j)+tsfcinc
                  print *,'lat=',rlat(i,j),' lon=',rlon(i,j)
               endif
            endif
            if(rlat(i,j).ge.23.0.and.rlat(i,j).lt.24.0) then
               if(rlon(i,j).gt.-110.0.and.rlon(i,j).lt.-105) then
                  rdata(i,j)=rdata(i,j)+tsfcinc
                  print *,'lat=',rlat(i,j),' lon=',rlon(i,j)
               endif
            endif
          endif
         enddo
      enddo
#endif
      return
      end
