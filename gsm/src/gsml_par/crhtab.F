      subroutine crhtab(rhcl,ier)
c---------------------------------------------------------------------
c..  cld-rh relations obtained from mitchell-hahn procedure, here read
c     cld/rh tuning tables for day 0,1,...,5 and merge into 1 file..
c                         .............k.a.c.   mar 93
c     use only one table (day 1) for all fcst hrs....k.a.c. feb 94
c...    4 cld types .... kac  feb96
c...    smooth out last bunch of bins of the tables...kac aug97
c    output:
c        rhcl - tuning tables for all forecast days
c        ier  - =1 if tables available.. =-1 if no tables
c
c     this routine is called only by master in mp mode
c--------------------------------------------------------------------
#include <paramodel.h>
crh1t parameter (mcld=3,nseal=2,ida=6,
cmcl4 parameter (mcld=4,nseal=2,ida=1,
      parameter (mcld=3,nseal=2,ida=1,
     1           nbin=100,nlon=2,nlat=4)
      dimension rhfd(nbin,nlon,nlat,mcld,nseal)
      dimension rrhfd(nbin,nlon,nlat,mcld,nseal)
      dimension rtnffd(nbin,nlon,nlat,mcld,nseal)
      dimension rrnffd(nbin,nlon,nlat,mcld,nseal)
      dimension rhcf(nbin,nlon,nlat,mcld,nseal)
      dimension rtnfcf(nbin,nlon,nlat,mcld,nseal)
      dimension kpts(nlon,nlat,mcld,nseal)
      dimension kkpts(nlon,nlat,mcld,nseal)
      dimension rhc(nlon,nlat,mcld,nseal)
      dimension rhcl (nbin,nlon,nlat,mcld,nseal,ida)
      dimension rhcla(nbin,nlon,nlat,mcld,nseal)
      dimension icdays(15),idate(4)
c...........................  begin here  ..............
      ier = 1
      do 8000 itim=1,ida
      icfq = 43 + itim-1
#ifdef RSM
      open(unit=icfq,file='r_tune',form='unformatted',status='old')
#else
      open(unit=icfq,file='tune',form='unformatted',status='old')
#endif
      rewind icfq
cmcl3       nclds=1,2,3 (l,m,h)..jsl=1,2 (land,sea)
cmcl4       mcld=1,2,3,4 (bl,l,m,h)
      binscl = 1./nbin
      do 1000 m=1,nseal
       do 1000 l=1,mcld
        do 1000 k=1,nlat
         do 1000 j=1,nlon
          do 1000 i=1,nbin
           rrhfd(i,j,k,l,m) = 0.
           rrnffd(i,j,k,l,m) = 0.
 1000 continue
      do 1001 m=1,nseal
       do 1001 l=1,mcld
        do 1001 k=1,nlat
         do 1001 j=1,nlon
          kkpts(j,k,l,m) = 0
 1001 continue
c....  read the data off the rotating file
      read (icfq,err=898,end=899) nbdayi,icdays
#ifdef DBG 
      print 11,nbdayi
#endif
      go to 123
  898 print 988,itim
#ifdef MP
#ifdef RSM
      call rmpabort
#else
      call mpabort
#endif
#else
      call abort
#endif
  899 print 989,itim
#ifdef MP
#ifdef RSM
      call rmpabort
#else
      call mpabort
#endif
#else
      call abort
#endif
123   continue
      do 53 ld=1,nbdayi
       id = icdays(ld) / 10000
       im = (icdays(ld)-id*10000) / 100
       iy = icdays(ld)-id*10000-im*100
#ifndef NOPRINT
       print 51,id,im,iy
#endif
   53 continue
      read (icfq,err=998,end=999) fhour,idate
#ifndef NOPRINT
      print 3003,idate,fhour,itim
#endif
      go to 223
  998 print 988,itim
#ifdef MP
#ifdef RSM
      call rmpabort
#else
      call mpabort
#endif
#else
      call abort
#endif
  999 print 989,itim
#ifdef MP
#ifdef RSM
      call rmpabort
#else
      call mpabort
#endif
#else
      call abort
#endif
223   continue
      do 1300 kd=1,nbdayi
       read (icfq) rhfd
       read (icfq) rtnffd
       read (icfq) kpts
       do 1002 m=1,nseal
        do 1002 l=1,mcld
         do 1002 k=1,nlat
          do 1002 j=1,nlon
           do 1002 i=1,nbin
            rrhfd(i,j,k,l,m) = rrhfd(i,j,k,l,m) + rhfd(i,j,k,l,m)
            rrnffd(i,j,k,l,m) = rrnffd(i,j,k,l,m)+rtnffd(i,j,k,l,m)
 1002  continue
       do 1003 m=1,nseal
        do 1003 l=1,mcld
         do 1003 k=1,nlat
          do 1003 j=1,nlon
           kkpts(j,k,l,m) = kkpts(j,k,l,m) + kpts(j,k,l,m)
 1003  continue
 1300 continue
c
      do 1004 m=1,nseal
       do 1004 l=1,mcld
        do 1004 k=1,nlat
         do 1004 j=1,nlon
          do 1004 i=1,nbin
           rhcf(i,j,k,l,m) = rrhfd(i,j,k,l,m)
           rtnfcf(i,j,k,l,m) = rrnffd(i,j,k,l,m)
 1004 continue
      do 1005 m=1,nseal
       do 1005 l=1,mcld
        do 1005 k=1,nlat
         do 1005 j=1,nlon
          kpts(j,k,l,m) = kkpts(j,k,l,m)
 1005 continue
c.....  compute the cumulative frequency distribution..
      do 200 n=1,nseal
       do 200 k=1,mcld
        do 200 l=1,nlat
         do 200 j=1,nlon
          do 190 i=2,nbin
           rhcf(i,j,l,k,n) = rhcf(i-1,j,l,k,n) + rhcf(i,j,l,k,n)
           rtnfcf(i,j,l,k,n)=rtnfcf(i-1,j,l,k,n) + rtnfcf(i,j,l,k,n)
  190     continue
  200 continue
      do 300 n=1,nseal
       do 300 l=1,nlat
        do 300 j=1,nlon
         do 300 k=1,mcld
          do 300 i=1,nbin
           if (kpts(j,l,k,n).gt.0) then
            rhcf(i,j,l,k,n) = rhcf(i,j,l,k,n) / kpts(j,l,k,n)
            rtnfcf(i,j,l,k,n) = rtnfcf(i,j,l,k,n) / kpts(j,l,k,n)
c
c...  cause we mix calculations of rh retune with cray and ibm words
c      the last value of rhcf is close to but ne 1.0,
c      so we reset it in order that the 360 loop gives compleat tabl
c...  rtnfcf caused couple of problems, seems to be ever so slightly
c      gt 1.0
c
            if (i.eq.nbin) then
             rhcf(i,j,l,k,n) = 1.0
            end if
            if (rtnfcf(i,j,l,k,n).ge.1.0) then
             rtnfcf(i,j,l,k,n) = 1.0
            end if
           else
            rhcf(i,j,l,k,n) = -0.1
            rtnfcf(i,j,l,k,n) = -0.1
           end if
  300 continue
#ifdef DBG
      do 255 nsl=1,nseal
       do 255 kcl=1,mcld
        print 264,kcl,nsl
        print 265,((kpts(i,l,kcl,nsl),i=1,nlon),l=1,nlat)
  255 continue
#endif
      do 360 nsl=1,nseal
       do 360 k=1,mcld
        do 360 l=1,nlat
         do 360 j=1,nlon
          if (kpts(j,l,k,nsl).le.0) go to 317
          do 320 i=1,nbin
           icrit = i
           if (rhcf(i,j,l,k,nsl).ge.rtnfcf(1,j,l,k,nsl)) go to 350
  320     continue
c... no critical rh
  317     icrit=-1
#ifdef DBG
          print 210,l,j,nsl
          print 202
          do 135 i=1,nbin
           print 203,rhcf(i,j,l,k,nsl),rtnfcf(i,j,l,k,nsl)
  135     continue
#endif
  350     rhc(j,l,k,nsl) = icrit * binscl
  360 continue
#ifdef DBG
      do 1210 nsl=1,nseal
       do 1210 k=1,mcld
        print 1221,k,nsl
        do 1210 l=1,nlat
         print 211,(rhc(j,l,k,nsl),j=1,nlon)
 1210 continue
#endif
      do 450 nsl=1,nseal
       do 450 ken=1,mcld
        do 450 l=1,nlat
         do 450 jl=1,nlon
          do 400 i=1,nbin
           rhcl(i,jl,l,ken,nsl,itim) = -0.1
  400     continue
  450 continue
      do 751 nsl=1,nseal
       do 751 ken=1,mcld
        do 751 l=1,nlat
         do 751 jl=1,nlon
          if (kpts(jl,l,ken,nsl).le.0) go to 751
          do 753 i=1,nbin
           do 755 j=1,nbin
            if (rhcf(j,jl,l,ken,nsl).ge.rtnfcf(i,jl,l,ken,nsl)) then
             rhcl(i,jl,l,ken,nsl,itim) = j*binscl
             go to 753
            end if
  755      continue
  753     continue
  751 continue
      do 3000 lon=1,nlon
       do 3000 lat=1,nlat
        do 3000 nc=1,mcld
         do 3000 nsl=1,nseal
         isat = 0
         do 67 it=1,nbin
          cfrac = binscl * (it-1)
          if (rhcl(it,lon,lat,nc,nsl,itim).lt.0.) then
           print 1941,it,nsl,nc,lat,lon
          end if
          if (it.lt.nbin.and.rtnfcf(it,lon,lat,nc,nsl).ge.1.) then
           if (isat.le.0) then
            isat = it
            rhsat = rhcl(it,lon,lat,nc,nsl,itim)
            clsat = cfrac
           end if
           rhcl(it,lon,lat,nc,nsl,itim) =
     1               rhsat + (1.-rhsat)*(cfrac-clsat)/(1.-clsat)
          end if
          if (it.eq.nbin) rhcl(it,lon,lat,nc,nsl,itim) = 1.
   67    continue
 3000 continue
c
c smooth out the table as it reaches rh=1.0, via linear interpolation
c      between location of rh ge .98 and the nbin bin (where rh=1.0)
c previously rh=1.0 occurred for many of the latter bins in the
c      table, thereby giving a cloud value of less then 1.0 for rh=1.0
c
      nb=nbin-2
      do 4000 lon=1,nlon
       do 4000 lat=1,nlat
        do 4000 nc=1,mcld
         do 4000 nsl=1,nseal
         do 4167 it=1,nbin
          rhcla(it,lon,lat,nc,nsl)=rhcl(it,lon,lat,nc,nsl,itim)
 4167    continue
         do 4067 it=1,nb
          ibs=it
          cfrac=binscl*ibs
          if (rhcl(it,lon,lat,nc,nsl,itim).ge..98) then
           do 4068 kt=ibs,nbin
            cstem=binscl*kt
            rhcla(kt,lon,lat,nc,nsl) =
     1       rhcl(ibs,lon,lat,nc,nsl,itim)+
     2       (rhcl(nbin,lon,lat,nc,nsl,itim)
     3                      -rhcl(ibs,lon,lat,nc,nsl,itim))*
     3       (cstem-cfrac)/(1.-cfrac)
 4068     continue
          go to 4000
          end if
 4067    continue
 4000 continue
c... restore table data to preferred location..
      do 4200 lon=1,nlon
       do 4200 lat=1,nlat
        do 4200 nc=1,mcld
         do 4200 nsl=1,nseal
          do 4200 it=1,nbin
           rhcl(it,lon,lat,nc,nsl,itim)= rhcla(it,lon,lat,nc,nsl)
 4200 continue
 8000 continue
      do 8001 ken=1,ida
       icfq = 42 + ken
       rewind icfq
 8001 continue
      return
   11 format(1h ,' days on file =',i5)
   51 format(1h ,' archv data from da,mo,yr=',3i4)
  202 format(1h0,' model rh ',' obs rtcld')
  203 format(2f10.2)
  210 format(1h ,' no crit rh for lat=',i3,' and lon band=',i3,
     1           ' land(=1) sea=',i3)
  211 format(1h ,15f6.2)
  264 format(1h ,' number of gg points used in each area..by latitude',
     1           '..for cloud type=',i4,'sealand=',i2)
  265 format(1h ,15i8)
  988 format(1h ,'....error reading tables for time=',i4)
  989 format(1h ,'....e.o.f reading tables for time=',i4)
 1221 format(1h0,' critical rh for lon,lat arrays for cld type=',i3,
     1           ' land(=1) sea=',i3)
 1941 format(1h ,' neg rhcl for it,nsl,nc,lat,lon=',5i4,'...stoppp..')
 3003 format(5x,'...last date/time and current itim',/,10x,
     1       4i15,f7.1,i6)
      end
