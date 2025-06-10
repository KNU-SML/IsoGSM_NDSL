      subroutine i2oini(proj ,orient ,truth ,cotru ,
     1                  delx ,dely ,rlat1 ,rlon1 ,ii, jj,
     2                  projo,oriento,trutho,cotruo,
     3                  x00,y00,x00o,y00o,
     4                  rlat1o,rlon1o,jo,glao,glon,glat,gaul)
#include <paramodel.h>
c
      logical*1 lflag
c
      dimension glat(jj),glon(ii),gaul(jj)
      dimension glao(jo)
      dimension gauo(jo)
c
      if(abs(proj).eq.4.) then  ! gaussian grid lat/lon data
       do j = 1,jj
         call gaulat(gaul,jj)   ! n-->s ( 0 --> 180)
         if(proj.eq.-4.) then      ! n--> s
            glat(j) = 90.-gaul(j)
         elseif(proj.eq.4) then    ! s--> n
            glat(jj-j+1) = 90.-gaul(j)
         endif
       enddo
      endif
      if(abs(projo).eq.4.) then  ! gaussian grid lat/lon data
         call gaulat(gauo,jo)   ! n-->s ( 0 --> 180)
         do j = 1,jo
            if(projo.eq.-4.) then      ! n--> s
               glao(j) = 90.-gauo(j)
            elseif(projo.eq.4) then    ! s--> n
               glao(jo-j+1) = 90.-gauo(j)
            endif
         enddo
      endif
      if(abs(proj).eq.3.) then  ! lat/lon data
         do j = 1,jj
            if(proj.eq.3.) then
               glat(j) = rlat1 + dely/1000. * float(j-1)
            elseif(proj.eq.-3.) then
               glat(j) = rlat1 - dely/1000. * float(j-1)
            endif
         enddo
      endif
c
      if(abs(proj).eq.3.) then
         do i = 1,ii
            glon(i) = rlon1 + delx/1000. * float(i-1)
            if(glon(i).lt.0.) glon(i) = glon(i) + 360.
         enddo
      elseif(abs(proj).eq.4.) then
         delx=360./float(ii)
         do i = 1,ii
            glon(i) = delx * float(i-1)
            if(glon(i).lt.0.) glon(i) = glon(i) + 360.
         enddo
      else
         rlon1r=rlon1
         call ll2xyr(1,proj,orient,truth,cotru,rlat1,rlon1,x00,y00)
#ifdef DBG
         print *, ' input grid rlat1 rlon1 x00 y00 ',
     1                     rlat1,rlon1,x00,y00
#endif
      endif
#ifdef DBG
      print*,'glat ',(glat(j),j=1,jj)
      print*,'glon ',(glon(i),i=1,ii)
#endif
c
      if(abs(projo).ne.3..and.projo.ne.4.) then
         call ll2xyr(1,projo,oriento,trutho,cotruo,rlat1o,rlon1o,
     1               x00o,y00o)
#ifdef DBG
         print *, ' output grid rlat1  rlon1  x00  y00 ',
     1                      rlat1o,rlon1o,x00o,y00o
#endif
      endif
c
      return
      end
