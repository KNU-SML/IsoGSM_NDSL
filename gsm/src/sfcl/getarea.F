      subroutine getarea(kgds,dlat,dlon,rslat,rnlat,wlon,elon,ijordr)
#include <paramodel.h>
c
c  get area of the grib record
c
      dimension kgds(22)
      logical ijordr
c
#ifdef DBG
      write(6,*) ' kgds( 1-10)=',(kgds(j),j= 1,10)
      write(6,*) ' kgds(11-20)=',(kgds(j),j=11,20)
      write(6,*) ' kgds(21-  )=',(kgds(j),j=21,22)
#endif
c
      if(kgds(1).eq.0) then
c
c  lat/lon grid
c
        write(6,*) 'lat/lon grid'
        dlat=float(kgds(10))/1000.0
        dlon=float(kgds( 9))/1000.0
        f0lon=float(kgds(5))/1000.0
        f0lat=float(kgds(4))/1000.0
        kgds11=kgds(11)
c
c  increase accuracy of dlon if possible (11-24-00)
c
        dlonx=360./float(kgds(2))
        if(abs(dlonx-dlon).lt.0.001) then
          if(abs(dlon).eq.abs(dlat)) then
            dlat=dlat/abs(dlat)*dlonx
          endif
          dlon=dlon/abs(dlon)*dlonx
        endif
        if(kgds11.ge.128) then
          wlon=f0lon-dlon*(kgds(2)-1)
          elon=f0lon
          if(dlon*kgds(2).gt.359.99) then
            wlon=f0lon-dlon*kgds(2)
          endif
          dlon=-dlon
          kgds11=kgds11-128
        else
          wlon=f0lon
          elon=f0lon+dlon*(kgds(2)-1)
          if(dlon*kgds(2).gt.359.99) then
            elon=f0lon+dlon*kgds(2)
          endif
        endif
        if(kgds11.ge.64) then
          rnlat=f0lat+dlat*(kgds(3)-1)
          rslat=f0lat
          kgds11=kgds11-64
        else
          rnlat=f0lat
          rslat=f0lat-dlat*(kgds(3)-1)
          dlat=-dlat
        endif
        if(kgds11.ge.32) then
          ijordr=.false.
        else
          ijordr=.true.
        endif

        if(wlon.gt.180.) wlon=wlon-360.
        if(elon.gt.180.) elon=elon-360.
        wlon=nint(wlon*1000.)/1000.
        elon=nint(elon*1000.)/1000.
        rslat=nint(rslat*1000.)/1000.
        rnlat=nint(rnlat*1000.)/1000.
        return
c
c  mercator projection
c
      elseif(kgds(1).eq.1) then
        write(6,*) 'mercator grid'
        write(6,*) 'cannot process'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
c
c  gnomonic projection
c
      elseif(kgds(1).eq.2) then
        write(6,*) 'gnomonic grid'
        write(6,*) 'error!! gnomonic projection not coded'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
c
c  lambert conformal
c
      elseif(kgds(1).eq.3) then
        write(6,*) 'lambert conformal'
        write(6,*) 'cannot process'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
      elseif(kgds(1).eq.4) then
c
c  gaussian grid
c
        write(6,*) 'gaussian grid'
        dlat=99.
        dlon=float(kgds( 9))/1000.0
        f0lon=float(kgds(5))/1000.0
        f0lat=99.
        kgds11=kgds(11)
        if(kgds11.ge.128) then
          wlon=f0lon
          elon=f0lon
          if(dlon*kgds(2).gt.359.99) then
            wlon=f0lon-dlon*kgds(2)
          endif
          dlon=-dlon
          kgds11=kgds11-128
        else
          wlon=f0lon
          elon=f0lon+dlon*(kgds(2)-1)
          if(dlon*kgds(2).gt.359.99) then
            elon=f0lon+dlon*kgds(2)
          endif
        endif
        if(kgds11.ge.64) then
          rnlat=99.
          rslat=99.
          kgds11=kgds11-64
        else
          rnlat=99.
          rslat=99.
          dlat=-99.
        endif
        if(kgds11.ge.32) then
          ijordr=.false.
        else
          ijordr=.true.
        endif
        return
c
c  polar strereographic
c
      elseif(kgds(1).eq.5) then
        write(6,*) 'polar stereographic grid'
        write(6,*) 'cannot process'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
        return
c
c  oblique lambert conformal
c
      elseif(kgds(1).eq.13) then
        write(6,*) 'oblique lambert conformal grid'
        write(6,*) 'cannot process'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
c
c  spherical coefficient
c
      elseif(kgds(1).eq.50) then
        write(6,*) 'spherical coefficient'
        write(6,*) 'cannot process'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
        return
c
c  space view perspective (orthographic grid)
c
      elseif(kgds(1).eq.90) then
        write(6,*) 'space view perspective grid'
        write(6,*) 'cannot process'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
        return
c
c  unknown projection.  abort.
c
      else
        write(6,*) 'error!! unknown map projection'
        write(6,*) 'kgds(1)=',kgds(1)
        print *,'error!! unknown map projection'
        print *,'kgds(1)=',kgds(1)
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
c
      return
      end
