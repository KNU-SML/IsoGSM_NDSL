#define PI pi
#define RERTH reath
      subroutine rgetarea(kgds,proj,orient,truth,cotru,
     1                   delx,dely,rlat1,rlat2,rlon1,rlon2)
#include <paramodel.h>
c
c  get area of the grib record
c
cfpp$ noconcur r
      dimension kgds(25)
      logical ijordr
c
      rad=pi_/180.
#ifdef DBG
      write(6,*) ' kgds( 1-10)=',(kgds(j),j= 1,10)
      write(6,*) ' kgds(11-20)=',(kgds(j),j=11,20)
      write(6,*) ' kgds(21-  )=',(kgds(j),j=21,25)
#endif
c
      idrt=kgds(1)
      rlat1=kgds(4)*1.e-3
      rlon1=kgds(5)*1.e-3

          if(idrt.eq.0) then
            proj=3
            rlat2=kgds(07)*1.e-3  ! latitude of last point
            rlon2=kgds(08)*1.e-3  ! longitude of last point
            delx=kgds(09)  ! dx (milledegree) on truth latitude
            dely=kgds(10)  ! dy (milledegree) on truth latitude
c           orient=kgds(11)  ! scanning mode ( 0 : n->s, 64 : s->n)
            if(kgds(11).eq.0.) proj = -3.
c  increase accuracy of dlon if possible (3-19-10)
            delxx=(rlon2-rlon1)/float(kgds(02)-1)
            if(abs(delxx-delx/1000.).lt.0.001) then
               if(abs(delx).eq.abs(dely)) then
                  dely=dely/abs(dely)*delxx*1000.
               endif
               delx=delx/abs(delx)*delxx*1000.
            endif
          endif
          if( idrt.eq.4 ) then      ! gaussian projection
            proj=4
            rlat2=kgds(07)*1.e-3  ! latitude of last point
            rlon2=kgds(08)*1.e-3  ! longitude of last point
            delx=kgds(09)  ! dx (milledegree) on truth latitude
            dely=kgds(10)  ! dy (milledegree) on truth latitude
c           orient=kgds(11)  ! scanning mode ( 0 : n->s, 64 : s->n)
            if(kgds(11).eq.0.) proj = -4.
          endif
          if( idrt.eq.1 ) then      ! mercater projection
            proj=0
            rlat2=kgds(07)*1.e-3  ! latitude of last point
            rlon2=kgds(08)*1.e-3  ! longitude of last point

c regional mercator projection that crosses 0E.
c hyodae/kana 6/3/2008
	      if(rlon1.gt.0.and.rlon1>rlon2) then
                  rlon2=rlon2+360.
               endif

            delx=(rlon2-rlon1)/float(kgds(02)-1)
c           delx=kgds(12)  ! dx (meter) on truth latitude
            dely=kgds(13)  ! dy (meter) on truth latitude
            orient=rlon1
            truth=kgds(09)*1.e-3 ! truth latitude
            cotru=truth          ! co-truth latitude
            delx=delx*cos(truth*rad)*rerth_*rad
          endif
          if( idrt.eq.5 ) then      ! polar projection
            truth=60.0              ! truth latitude
            cotru=60.0              ! co-truth latitude
            orient=kgds(07)*1.e-3   ! orientation
            delx=kgds(08)  ! dx (meter) on 60 deg
            dely=kgds(09)  ! dy (meter) on 60 deg
            iproj=kgds(10) ! polar projection (first bit 0:north;1=south)
            if( iproj.eq.0 ) then
              proj=1.0
            else
              proj=-1.0
            endif
          endif
          if( idrt.eq.3 ) then      ! lambert projection
            orient=kgds(07)*1.e-3   ! orientation
            delx=kgds(08)  ! dx (meter) on 60 deg
            dely=kgds(09)  ! dy (meter) on 60 deg
            iproj=kgds(10) !  projection (first bit 0:north;1=south)
            if( iproj.eq.0 ) then
              proj=2.0
            else
              proj=-2.0
            endif
            truth=kgds(12)*1.e-3 ! the 1st lat from pole to cut
            cotru=kgds(13)*1.e-3 ! the 2nd lat from pole to cut
          endif
c
      return
      end
