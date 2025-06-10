      subroutine fixrd(lugb,fngrib,kkpds,
     &                 iy,im,id,ih,fh,imdata,jmdata,
     &                 imax,jmax,kkpds5,kkgds1,
     &                 proji,orienti,truthi,cotrui,
     &                 delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i,
     &                 dlon,dlat,wlon,rnlat,
     &                 data,iret)
#include <paramodel.h>
#ifdef CRAY_THREAD
cfpp$ noconcur r
#endif
c
c read in grib climatology files.
c
c interpolate climatology to the dates
c
c grib file should allow all the necessary parameters to be extracted fr
c the description records.
c
c
c  nvalid:  analysis later than (current date - nvalid) is regarded as
c           valid for current analysis
c
#ifdef ADJSOILW
      parameter(nvalid=1)
#else
      parameter(nvalid=5)
#endif
c
c  iret=0  ... successfully read
c      =-1 ... record with unmatched dates exist
c
      character*128 fngrib
c
      real data(imdata*jmdata)
      real, allocatable :: data_next(:)
c
      logical*1 lbms(imdata*jmdata)
c
      parameter(mbuf=1024*128*64)
      character*1 cbuf(mbuf)
c
      integer kpds(25),kgds(22),kens(5)
      integer jpds(25),jgds(22),jens(5)
      integer kkpds(25),kpds0(25)
c
      character*80 asgnstr
c
      data msk1/3840000/,msk2/2400000/
c
#ifdef REAL4_W3LIB
      integer*4 lugb4,msk14,msk24,mnum4,mbuf4
      integer*4 nlen4,nnum4,iret4
      integer*4 ndata4
      real*4, allocatable :: data4(:)
      integer*4 lskip4,lgrib4,lret4
      integer*4 n4,jpds4(25),jgds4(22),jens4(5)
      integer*4 k4,kpds4(25),kgds4(22),kens4(5)
#endif
c
#ifdef FASTBAREAD
      character, allocatable :: bbuf(:)
#endif
c
      logical lclim
c
c julian day of the middle of each month
c
      real dayhf(13)
      data dayhf/ 15.5, 45.0, 74.5,105.0,135.5,166.0,
     1           196.5,227.5,258.0,288.5,319.0,349.5,380.5/
      save dayhf
c
#ifdef REAL4_W3LIB
      allocate (data4(imdata*jmdata))
#endif
c
      call vday(iy,im,id,ih,fh,jy,jm,jd,jh,rjday)
c
      close(lugb)
      call numchar(fngrib,nch)
#ifdef ASSIGN
      if(lugb.lt.10) then
        write(asgnstr,'(22hassign -s unblocked u:,i1)') lugb
      else
        write(asgnstr,'(22hassign -s unblocked u:,i2)') lugb
      endif
      call assign('assign -R')
      call assign(asgnstr)
      open(unit=lugb,file=fngrib(1:nch),status='old',form='unformatted',
     1     err=910)
#else
#ifdef FASTBAREAD
      call nainit(lugb,fngrib(1:nch),isize,iret)
      allocate (bbuf(isize))
      call naopen(lugb,fngrib(1:nch),bbuf,isize,iret)
#else
      call baopen(lugb,fngrib(1:nch),iret)
#endif
      if(iret.ne.0) go to 910
#endif
      go to 911
  910 continue
      print *,'fixrd: error in opening file ',fngrib(1:nch)
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
  911 continue
      write(6,*) ' file ',fngrib(1:nch),' opened. unit=',lugb
c
c  get grib index buffer
c
      mnum=0
#ifdef REAL4_W3LIB
      lugb4=lugb
      msk14=msk1
      msk24=msk2
      mnum4=mnum
      mbuf4=mbuf
#ifdef FASTBAREAD
      call ngetgir(bbuf,isize,msk14,msk24,mnum4,mbuf4,
     1            cbuf,nlen4,nnum4,iret4)
#else
      call getgir(lugb4,msk14,msk24,mnum4,mbuf4,
     1            cbuf,nlen4,nnum4,iret4)
#endif
      nlen=nlen4
      nnum=nnum4
      iret=iret4
#else
#ifdef FASTBAREAD
      call ngetgir(bbuf,isize,msk1,msk2,mnum,mbuf,cbuf,nlen,nnum,iret)
#else
      call getgir(lugb,msk1,msk2,mnum,mbuf,cbuf,nlen,nnum,iret)
#endif
#endif
#ifdef DBG
      print *,'nlen=',nlen,' nnum=',nnum
#endif
      if(iret.ne.0) then
        print *,'error.  cbuf length too short in ngetgir'
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
      if(nnum.eq.0) then
        print *,'error.  not a grib file. detected in ngetgir'
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
      if(nlen.eq.0) then
        print *,'error.  nlen=0.  detected in ngetgir'
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
      if(nlen.gt.imdata*jmdata) then
        print *,'error.  nlen .gt. imdata*jmdata  detected in ngetgir'
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
c  find file type climatology or analysis
c
      do i=1,25
        jpds(i)=-1
      enddo
      do i=1,22
        jgds(i)=-1
      enddo
      do i=1,5
        jens(i)=-1
      enddo
      do i=1,25
        jpds(i)=kkpds(i)
      enddo
c
c fix for ecmwf grib
c
      if ( (jpds(5).eq.11.or.jpds(5).eq.91) .and. 
     &     (jpds(7).eq.-1.or.jpds(7).eq.0) ) then
        jpds(5)=-1
        jpds(7)=0
        write(6,*) 'Changing jpds5 and jpds7 to allow ',
     1             'for ECMWF and NCEP ice and sst'
      endif
c
#ifdef DBG
      write(6,*) ' Searching the following fields with kpds=.'
      write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
      write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
      write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
      n=0
#ifdef REAL4_W3LIB
      nlen4=nlen
      nnum4=nnum
      n4=n
      do i=1,25
        jpds4(i)=jpds(i)
      enddo
      do i=1,22
        jgds4(i)=jgds(i)
      enddo
      do i=1,5
        jens4(i)=jens(i)
      enddo
      call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,
     &             k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
      k=k4
      do i=1,25
        kpds(i)=kpds4(i)
      enddo
      do i=1,22
        kgds(i)=kgds4(i)
      enddo
      do i=1,5
        kens(i)=kens4(i)
      enddo
      lskip=lskip4
      lgrib=lgrib4
      iret=iret4
#else
      call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,
     &             k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
c
#ifdef DBG
      write(6,*) 'first grib record in the  grib file',fngrib(1:nch)
      write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
      write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
      write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
c
      if(lgrib.eq.0) then
        write(6,*) ' error in getgbss.  No matching records.'
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
      do i=1,25
        kpds0(i)=kpds(i)
      enddo
      kpds0(4)=-1
      kpds0(18)=-1
c
clu_rev: manipulate kpds(13) and kpds(15) for fixed fields
c
c neet to clean up  mk      
c
      if(kpds(5).eq.236 .or. kpds(5).eq.255) then
         kpds(13) = 4
         kpds(15) = 1
         kpds(16) = 51
         print *,' quick fix for noah fixed fields'
      endif
c
      if(kpds(16).eq.51) then
        write(6,*) ' climatology file.'
        lclim=.true.
      else
        write(6,*) ' analysis file.'
        lclim=.false.
      endif
c
c   handling climatology file
c
      if(lclim) then
c
c       find average type
c         weekly,biweekly,monthly,seasonal,annual
c
c  kpds(13)=4 & kpds(15)=1 .. annual mean
c  kpds(13)=3 & kpds(15)=1 .. monthly mean
c  kpds(13)=2 & kpds(15)=7 .. weekly mean
c  kpds(13)=2 & kpds(15)=14.. bi-weekly mean
c
        if(kpds(13).eq.2.and.kpds(15).eq.7) then
          write(6,*) ' this is weekly mean climatology'
          write(6,*) ' cannot process.'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
        elseif(kpds(13).eq.2.and.kpds(15).eq.14) then
          write(6,*) ' this is bi-weekly mean climatology'
          write(6,*) ' cannot process.'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
        elseif(kpds(13).eq.3.and.kpds(15).le.1) then
          write(6,*) ' this is monthly mean climatology'
          monend=12
          do mm=1,monend
            mmm=mm
            mmp=mm+1
            if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
              mon1=mmm
              mon2=mmp
              go to 20
            endif
          enddo
          print *,'wrong rjday',rjday
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
   20     continue
          ijmax=0
          do nn=1,2
            if(nn.eq.2) then
              allocate (data_next(ijmax))
            endif
            lskip = -1
            n=0
            do i=1,25
              jpds(i)=kpds0(i)
            enddo
            jpds(24)=-1
            jpds(25)=-1
            if(nn.eq.1) jpds( 9)=mon1
            if(nn.eq.2) jpds( 9)=mon2
            if(jpds(9).eq.13) jpds(9)=1
#ifdef DBG
            write(6,*) ' Searching the following fields with kpds=.'
            write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
            write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
            write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
#ifdef REAL4_W3LIB
            nlen4=nlen
            nnum4=nnum
            n4=n
            do i=1,25
              jpds4(i)=jpds(i)
            enddo
            do i=1,22
              jgds4(i)=jgds(i)
            enddo
            do i=1,5
              jens4(i)=jens(i)
            enddo
            call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,
     &                   k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
            k=k4
            do i=1,25
              kpds(i)=kpds4(i)
            enddo
            do i=1,22
              kgds(i)=kgds4(i)
            enddo
            do i=1,5
              kens(i)=kens4(i)
            enddo
            lskip=lskip4
            lgrib=lgrib4
            iret=iret4
#else
            call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,
     &                   k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
c
            if(lgrib.eq.0) then
#ifdef DBG
              write(6,*) ' No matching record found.'
              write(6,*) ' The last header records read:'
              write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
              write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
              write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
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
            write(6,*) ' Matching record found for mon=',jpds(9)
#ifdef DBG
            write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
            write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
            write(6,*) ' kpds(21-  )=',(kpds(j),j=21,22)
#endif
            do i=1,imdata*jmdata
              lbms(i)=.true.
            enddo
#ifdef REAL4_W3LIB
            lugb4=lugb
#ifdef FASTBAREAD
            call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2                data4,6)
#else
            call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2                data4,6)
#endif
            do i=1,25
              kpds(i)=kpds4(i)
            enddo
            do i=1,22
              kgds(i)=kgds4(i)
            enddo
            ndata=ndata4
            if(nn.eq.1) then
              do i=1,ndata
                data(i)=data4(i)
              enddo
            else
              do i=1,ndata
                data_next(i)=data4(i)
              enddo
            endif
#else
            if(nn.eq.1) then
#ifdef FASTBAREAD
              call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data     ,6)
#else
              call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data     ,6)
#endif
            else
#ifdef FASTBAREAD
              call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data_next,6)
#else
              call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data_next,6)
#endif
            endif
#endif
            if(ndata.eq.0) then
              write(6,*) ' error in nrdgb'
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
            imax=kgds(2)
            jmax=kgds(3)
            ijmax=imax*jmax
#ifdef DBG
            write(6,*) 'imax,jmax,ijmax=',imax,jmax,ijmax
#endif
          enddo
c
          wei1=(dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
          wei2=(rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))
          if(mon2.eq.13) mon2=1
#ifdef DBG
          print *,'rjday,mon1,mon2,wei1,wei2=',
     1             rjday,mon1,mon2,wei1,wei2
#endif
          do i=1,ijmax
            data(i)=wei1*data(i)+wei2*data_next(i)
          enddo
          deallocate (data_next)
        elseif(kpds(13).eq.4.and.kpds(15).eq.3) then
          write(6,*) ' this is seasonal mean climatology'
          monend=4
          is=im/3+1
          if(is.eq.5) is=1
          is1=mon1/3+1
          if(is1.eq.5) is1=1
          is2=mon2/3+1
          if(is2.eq.5) is2=1
          do mm=1,monend
            mmm=mm*3-2
            mmp=(mm+1)*3-2
            if(rjday.ge.dayhf(mmm).and.rjday.lt.dayhf(mmp)) then
              mon1=mmm
              mon2=mmp
              go to 30
            endif
          enddo
          print *,'wrong rjday',rjday
#ifdef MP
#ifdef RSM
          call rmpabort
#else
          call mpabort
#endif
#else
          call abort
#endif
   30     continue
          ijmax=0
          do nn=1,2
            if(nn.eq.2) then
              allocate (data_next(ijmax))
            endif
            do i=1,25
              jpds(i)=kpds0(i)
            enddo
            jpds(24)=-1
            jpds(25)=-1
            n=0
            if(nn.eq.1) then
              isx=is1
            else
              isx=is2
            endif
            if(isx.eq.1) jpds(9)=12
            if(isx.eq.2) jpds(9)=3
            if(isx.eq.3) jpds(9)=6
            if(isx.eq.4) jpds(9)=9
            if(jpds(9).eq.13) jpds(9)=1
#ifdef DBG
            write(6,*) ' Searching the following fields with kpds=.'
            write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
            write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
            write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
#ifdef REAL4_W3LIB
            nlen4=nlen
            nnum4=nnum
            n4=n
            do i=1,25
              jpds4(i)=jpds(i)
            enddo
            do i=1,22
              jgds4(i)=jgds(i)
            enddo
            do i=1,5
              jens4(i)=jens(i)
            enddo
            call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,
     &                   k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
            k=k4
            do i=1,25
              kpds(i)=kpds4(i)
            enddo
            do i=1,22
              kgds(i)=kgds4(i)
            enddo
            do i=1,5
              kens(i)=kens4(i)
            enddo
            lskip=lskip4
            lgrib=lgrib4
            iret=iret4
#else
            call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,
     &                   k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
c
            if(lgrib.eq.0) then
#ifdef DBG
              write(6,*) ' No matching record found.'
              write(6,*) ' The last header records read:'
              write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
              write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
              write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
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
            do i=1,imdata*jmdata
              lbms(i)=.true.
            enddo
#ifdef REAL4_W3LIB
            lugb4=lugb
#ifdef FASTBAREAD
            call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2                data4,6)
#else
            call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2                data4,6)
#endif
            do i=1,25
              kpds(i)=kpds4(i)
            enddo
            do i=1,22
              kgds(i)=kgds4(i)
            enddo
            ndata=ndata4
            if(nn.eq.1) then
              do i=1,ndata
                data(i)=data4(i)
              enddo
            else
              do i=1,ndata
                data_next(i)=data4(i)
              enddo
            endif
#else
            if(nn.eq.1) then
#ifdef FASTBAREAD
              call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data     ,6)
#else
              call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data     ,6)
#endif
            else
#ifdef FASTBAREAD
              call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data_next,6)
#else
              call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data_next,6)
#endif
            endif
#endif
            if(ndata.eq.0) then
              write(6,*) ' error in nrdgb'
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
            imax=kgds(2)
            jmax=kgds(3)
            ijmax=imax*jmax
          enddo
          wei1=(dayhf(mon2)-rjday)/(dayhf(mon2)-dayhf(mon1))
          wei2=(rjday-dayhf(mon1))/(dayhf(mon2)-dayhf(mon1))
          if(mon2.eq.13) mon2=1
#ifdef DBG
          print *,'rjday=',rjday
          print *,'mon1 =',mon1 ,' mon2=',mon2
          print *,'wei1 =',wei1 ,' wei2=',wei2
          print *,'ses1 =', is1 ,' ses2=', is2
#endif
          do i=1,ijmax
            data(i)=wei1*data(i)+wei2*data_next(i)
          enddo
          deallocate (data_next)
        elseif(kpds(13).eq.4.and.kpds(15).eq.1) then
          write(6,*) ' this is annual mean climatology'
          monend=-1
#ifdef REAL4_W3LIB
          lugb4=lugb
#ifdef FASTBAREAD
          call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2              data4,6)
#else
          call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     2              data4,6)
#endif
          do i=1,25
            kpds(i)=kpds4(i)
          enddo
          do i=1,22
            kgds(i)=kgds4(i)
          enddo
          ndata=ndata4
          do i=1,ndata
            data(i)=data4(i)
          enddo
#else
#ifdef FASTBAREAD
          call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#else
          call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#endif
#endif
          if(ndata.eq.0) then
            write(6,*) ' error in nrdgb'
            write(6,*) ' kpds=',kpds
            write(6,*) ' kgds=',kgds
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
          imax=kgds(2)
          jmax=kgds(3)
          ijmax=imax*jmax
        else
          write(6,*) ' climatology file average period unknown.'
          write(6,*) ' kpds(13)=',kpds(13),' kpds(15)=',kpds(15)
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
      else
c-------------------------------------------------------------------
c
c  handling analysis file
c
c  find record for the given hour/day/month/year
c
        monend=0
        nrept=0
        do i=1,25
          kpds(i)=kpds0(i)
        enddo
        iyr=jy
        imo=jm
        idy=jd
        ihr=jh
   50   continue
        jpds( 8)=iyr
        jpds( 9)=imo
        jpds(10)=idy
        if(ihr.eq.-1) then
c
c t.r.    jpds(11)=0
          jpds(11)=-1
        else
c
c t.r.    jpds(11)=ihr
          jpds(11)=-1
        endif
#ifdef DBG
        write(6,*) ' Searching the following fields with kpds=.'
        write(6,*) ' jpds( 1-10)=',(jpds(j),j= 1,10)
        write(6,*) ' jpds(11-20)=',(jpds(j),j=11,20)
        write(6,*) ' jpds(21-  )=',(jpds(j),j=21,25)
#endif
        n=0
#ifdef REAL4_W3LIB
        nlen4=nlen
        nnum4=nnum
        n4=n
        do i=1,25
          jpds4(i)=jpds(i)
        enddo
        do i=1,22
          jgds4(i)=jgds(i)
        enddo
        do i=1,5
          jens4(i)=jens(i)
        enddo
        call getgbss(cbuf,nlen4,nnum4,n4,jpds4,jgds4,jens4,
     &               k4,kpds4,kgds4,kens4,lskip4,lgrib4,iret4)
        k=k4
        do i=1,25
          kpds(i)=kpds4(i)
        enddo
        do i=1,22
          kgds(i)=kgds4(i)
        enddo
        do i=1,5
          kens(i)=kens4(i)
        enddo
        lskip=lskip4
        lgrib=lgrib4
        iret=iret4
#else
        call getgbss(cbuf,nlen,nnum,n,jpds,jgds,jens,
     &               k,kpds,kgds,kens,lskip,lgrib,iret)
#endif
c
        if(lgrib.ne.0) then
          if(nrept.le.iabs(nvalid).and.nrept.ne.0) then
            write(6,*) '<warning:cycl> grib record date does not match',
     1                 ' but within ',iabs(nvalid),' days.'
            write(6,*) '<warning:cycl> assume it is valid.'
          endif
          write(6,*) ' Matching record found'
#ifdef DBG
          write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
          write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
          write(6,*) ' kpds(21-  )=',(kpds(j),j=21,22)
          write(6,*) ' Now reading grib record'
#endif
#ifdef REAL4_W3LIB
          lugb4=lugb
#ifdef FASTBAREAD
          call nrdgb(bbuf,isize,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,
     1              data4,6)
#else
          call rdgb(lugb4,lgrib4,lskip4,kpds4,kgds4,ndata4,lbms,data4,6)
#endif
          do i=1,25
            kpds(i)=kpds4(i)
          enddo
          do i=1,22
            kgds(i)=kgds4(i)
          enddo
          ndata=ndata4
          do i=1,ndata
            data(i)=data4(i)
          enddo
#else
#ifdef FASTBAREAD
          call nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#else
          call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,data,6)
#endif
#endif
          if(ndata.eq.0) then
            write(6,*) ' error in nrdgb'
            write(6,*) ' kpds=',kpds
            write(6,*) ' kgds=',kgds
            write(6,*) ' lgrib,lskip=',lgrib,lskip
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
#ifdef DBG
          write(6,*) ' matching grib record successfully read.'
#endif
          imax=kgds(2)
          jmax=kgds(3)
          ijmax=imax*jmax
        else
          if(nrept.eq.0) then
            write(6,*) ' no matching dates found.  start searching',
     1                 ' the nearest matching past date.'
          endif
c
c  no matching ih found. search nearest hour
c
          if(ihr.eq.6) then
            ihr=0
            write(6,*) ' hour set to 0'
            go to 50
          elseif(ihr.eq.12) then
            ihr=0
            write(6,*) ' hour set to 0'
            go to 50
          elseif(ihr.eq.18) then
            ihr=12
            write(6,*) ' hour set to 12'
            go to 50
          elseif(ihr.eq.0.or.ihr.eq.-1) then
            idy=idy-1
            if(idy.eq.0) then
              imo=imo-1
              if(imo.eq.0) then
                iyr=iyr-1
                if(iyr.lt.0) iyr=99
                imo=12
              endif
              idy=31
              if(imo.eq.4.or.imo.eq.6.or.imo.eq.9.or.imo.eq.11) idy=30
              if(imo.eq.2) then
                if(mod(iyr,4).eq.0) then
                  idy=29
c
c                 Check 100/400 year leap
                  if(mod(jy,100).eq.0.and.mod(jy,400).ne.0) then
                     mondy=28
                  endif
c
                else
                  idy=28
                endif
              endif
            endif
            ihr=-1
            write(6,*) ' dates decremented to:',iyr,imo,idy,ihr
            nrept=nrept+1
            if(nrept.gt.nvalid) then
              iret=-1
              write(6,*) ' <warning:surface> search range exceeded.'
              write(6,*) ' <warning:surface> terminating search.'
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
            go to 50
          else
#ifdef DBG
            write(6,*) ' Matching date could not be fond.'
            write(6,*) ' The last header records read:'
            write(6,*) ' kpds( 1-10)=',(kpds(j),j= 1,10)
            write(6,*) ' kpds(11-20)=',(kpds(j),j=11,20)
            write(6,*) ' kpds(21-  )=',(kpds(j),j=21,25)
#endif
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
        endif
      endif
c
c  finished reading grib file. 
c
c  substitution to the output array
c
   80 continue
      call rmaxmin(data,imax*jmax,1,1,1,'input data')
c
#ifdef DBG
c     if(kkpds(5).eq.225) then
c       print *,'input uninterpolated vegtype'
c       call nntprt(data,imax,jmax)
c     endif
c     if(kkpds(5).eq.87) then
c       print *,'input uninterpolated vegfrac'
c       call nntprt(data,imax,jmax)
c     endif
c     if(kkpds(5).eq.230) then
c       print *,'input uninterpolated soiltype'
c       call nntprt(data,imax,jmax)
c     endif
#endif
c
#ifdef DBG
      write(6,*) 'imax,jmax,ijmax=',imax,jmax,ijmax
#endif
      proji=0.
      orienti=0.
      truthi=0.
      cotrui=0.
      call rgetarea(kgds,proji,orienti,truthi,cotrui,
     1              delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i)
#ifdef DBG
      write(6,*) 'proji,orienti,truthi,cotrui='
      write(6,*)  proji,orienti,truthi,cotrui
      write(6,*) 'delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i='
      write(6,*)  delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i
#endif
#ifndef RSM
      call getarea(kgds,dlat,dlon,rslat,rnlat,wlon,elon,ijordr)
#ifdef DBG
      write(6,*) 'dlat,dlon,rslat,rnlat,wlon,elon,ijordr ='
      write(6,*)  dlat,dlon,rslat,rnlat,wlon,elon,ijordr
#endif
      call subst(data,imax,jmax,dlon,dlat,ijordr)
#endif
c
      write(6,*) 'imax,jmax of grib=',imax,jmax
c
      kkpds5=kkpds(5)
      kkgds1=kgds(1)
c
#ifdef REAL4_W3LIB
      deallocate (data4)
#endif
#ifdef FASTBAREAD
      deallocate (bbuf)
#endif
      close(lugb)
c
      return
      end
