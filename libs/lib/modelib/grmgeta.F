      subroutine grmgeta(lugb,nkeys,fngrib,ifile,irec,maxrec,
     1                  kpd5,kpd6,ipres,iy,im,id,ih,ifh,
     2                  io,jo,nrec,datao,fgribn,lgot,lptr,jret)
c
c$$$  subprogram program documentation block
c
c sub program:  grmget      searches multiple grib files for matching 
c                           multiple fields and retrievs the fields
c   prgmmr: kanamitsu        org: w/np51     date: 96-03-15
c
c abstract: program searches multiple grib files that matches multiple
c           variable id, level and date keys and retrieve multiple 
c           unpacked fields interpolated to a specified grid.
c
c program history log:
c   92-03-15  kanamitsu
c
c usage:  call grgrep(lugb,nkeys,fngrib,maxrec,kpd5,kpd6,ipres,
c                     iy,im,id,ih,ifh,
c                     io,jo,nrec,datao,fgribn,lgot,lptr,jret)
c   input argument list
c      lugb-  integer.  unit number used to read multiple grib file
c     nkeys-  integer.  number of keys
c    fngrib-  character*80 array of any length.  end of array should be
c             identified as ten or more blanks.  grib file name
c    ifile -  index of the file to be read.  modified on output, index of
c             file last read.
c    irec  -  records to be skipped.  zero to open
c             the first grib file and create index buffer.  note that
c             if field is not found, it opens new file even irec>0,
c             so user needs to be aware which file you are currently reading.  
c             modified on output, last read record number.
c    maxrec-  maximum number of record expected to be retrieved
c             this is used to avoid going out of dimension, in case
c             you do not know how many records matches the key.
c             if this is set to one, program returns the field as soon as
c             it finds one.  convenient for processing all the fields
c             one by one without knowing what they are.
c      kpd5-  integer array of length nkeys.  variable id
c      kpd6-  integer array of length nkeys.  vertical level id
c     ipres-  integer array of length nkeys.  pressure or level
c        iy-  integer array of length nkeys.  year to be searched
c        im-  integer array of length nkeys.  month to be searched
c        id-  integer array of length nkeys.  day to be searched
c        ih-  integer array of length nkeys.  hour to be searched
c       ifh-  integer array of length nkeys.  forecast hour to be searched
c
c    a value of -1 can be specified as a wild card for the varables
c    kpd5,kpd6,ipres,iy,im,id,ih,ifh.
c
c   output argument list
c        io-  integer. output array x-dimension
c        jo-  integer. output array y-dimension
c      nrec-  actual number of records retrieved
c     datao-  real array of (io*jo*nrec).  output fields
c    fgribn-  character*80 array of length nkeys.  grib file names where
c             the field is found
c      lgot-  logical array of length nkeys.  true for field found, false for not
c             found
c      lptr-  integer.  print output unit number
c      jret-  integer.  return code. 0=normal return. -1=field not found
c
c  subroutine called:
c    assign - file assign for cray.  dummy routine for other
c    abort  - terminates program with abort.  system routine
c    getgir - creates grib index buffer
c    getgbss- searchs grib index buffer for field
c    maxmin - find max/min and prints
c    rdgb   - read grib file and unpacks
c    getarea- obtain area definition from grib description record
c    ll2ll  - interploation from one lat/lon grid to the other lat/lon grid
c    subst  - substitute grib unpacked field in to standard format array
c             (from zero degree east to west, north to south)
c
c  attributes:
c    language:  fortran 77
c
c    machine: sgi cray
c
      character*80 fngrib(*)
c
      character*80 fgribn(*)
c
      dimension iy(*),im(*),id(*),ih(*),ifh(*)
      dimension kpd5(*),kpd6(*),ipres(*)
c
      parameter(mdata=360*181)
c
      logical lopened
c
      logical lbms(mdata)
c
      real datao(*)
c
      real data(mdata),work(mdata)
c
      data msk1/32000/,msk2/4000/
      parameter(mbuf=1024*128*64)
      character*1 cbuf(mbuf)
c
      integer kpds(25),kgds(22),kens(5)
      integer jpds(25),jgds(22),jens(5)
c
      character*80 asgnstr
c
      logical lgot(*)
c
      logical lwork1,lwork2
c
#include <machine.h>
c
#ifndef CRAY
      integer*4 lugb4,msk14,msk24,mnum4,mbuf4
      integer*4 nlen4,nnum4,iret4
      integer*4 ndata4
      real*4 data4(mdata*2)
      integer*4 lskip4,lgrib4,lret4
      integer*4 n4,jpds4(25),jgds4(22),jens4(5)
      integer*4 k4,kpds4(25),kgds4(22),kens4(5)
#endif
c
      write(lptr,*) '============================================='
c
      write(lptr,*) 'kpd5,kpd6,ipres,iy,im,id,ih,ifh'
      do n=1,nkeys
        write(lptr,*) kpd5(n),kpd6(n),ipres(n),
     1                iy(n),im(n),id(n),ih(n),ifh(n)
      enddo
c
      do mm=1,nkeys
        lgot(mm)=.false.
      enddo
c
c  kount=count for how may key elements are found
c  nrec=number of fields retrieved.  this is different from nkeys when
c       wild card -1 is used in the key element
c  nfile=index of the file to be read
c  kfile=number of grib files read
c
      nrec=0
      kount=0
      nfile=ifile
      kfile=1
 1000 continue
      if(kfile.gt.1) then
        irec=0
      endif
c
c  exit for exhosting grib files
c
      if(fngrib(nfile)(1:10).eq.'          ') then
        write(6,*) 'field not found'
        mm=1
        dowhile(lgot(mm))
          mm=mm+1
        enddo
        if(mm.le.nkeys) then
          write(6,*) 'search keys:'
          write(6,*) 'idvar=',kpd5(mm),' levtyp=',kpd6(mm),
     1               ' lev=',ipres(mm),' year=',iy(mm),' month=',im(mm),
     2               ' day=',id(mm),' hour=',ih(mm),
     3               ' fhour=',ifh(mm)
          write(6,310) fngrib(nfile-1)
  310     format('last searched file=',a50)
          write(lptr,*) 'other fields not found are...'
          do n=1,nkeys
            if(.not.lgot(n)) then
              write(lptr,*) 'idvar=',kpd5(mm),' levtyp=',kpd6(mm),
     1                      ' lev=',ipres(mm),' year=',iy(mm),
     2                      ' month=',im(mm),' day=',id(mm),
     3                      ' hour=',ih(mm),' fhour=',ifh(mm)
            endif
          enddo
        endif
        stop
      endif
c
c  open grib input file
c
      inquire(file=fngrib(nfile),opened=lopened)
      if(.not.lopened) then
        if(lugb.lt.10) then
          write(asgnstr,'(22hassign -s unblocked u:,i1,)') lugb
        else
          write(asgnstr,'(22hassign -s unblocked u:,i2,)') lugb
        endif
        call assign(asgnstr)
        open(lugb,file=fngrib(nfile),status='old',form='unformatted',
     1       err=910)
        go to 911
  910   continue
        write(6,*) ' error in opening file '
        print *,'error in opening file '
        write(lptr,*) ' fn=',fngrib(nfile)
        write(lptr,*) ' unit=',lugb 
c       call abort
        jret=99
        return
  911   continue
        write(lptr,*) ' grib file opened'
        write(lptr,*) ' fn=',fngrib(nfile)
        write(lptr,*) ' unit=',lugb 
c
c  create grib index buffer
c
        mnum=0
c
#ifndef CRAY
        lugb4=lugb
        msk14=msk1
        msk24=msk2
        mnum4=mnum
        mbuf4=mbuf
        call getgir(lugb4,msk14,msk24,mnum4,mbuf4,
     1              cbuf,nlen4,nnum4,iret4)
        nlen=nlen4
        nnum=nnum4
        iret=iret4
#else
        call getgir(lugb,msk1,msk2,mnum,mbuf,cbuf,nlen,nnum,iret)
#endif
c
        write(lptr,*) 'nlen=',nlen,' nnum=',nnum
        if(iret.ne.0) then
          write(6,*) 'error.  cbuf length too short in getgir'
          print *,'error.  cbuf length too short in getgir'
          call abort
        endif
        if(nnum.eq.0) then
          write(6,*) 'error. not a grib file. detected in getgir'
          print *,'error.  not a grib file. detected in getgir'
          call abort
        endif
        if(nlen.eq.0) then
          write(6,*) 'error. nlen=0. detected in getgir'
          print *,'error.  nlen=0.  detected in getgir'
          call abort
        endif
      endif
c
c  find matching field in this file
c
      do mm=1,nkeys
c
      if(lgot(mm)) go to 2000
c
      nwild=0
 3000 continue
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
c
      jpds( 5)=kpd5(mm)
      jpds( 6)=kpd6(mm)
      jpds( 7)=ipres(mm)
      jpds( 8)=iy(mm)
      jpds( 9)=im(mm)
      jpds(10)=id(mm)
      jpds(11)=ih(mm)
      jpds(14)=ifh(mm)
c
      if(nwild.eq.0) then
        if(mm.eq.1) then
          n=irec
        else
          mn=mm-1
          lwork1=kpd5(mn).eq.-1.or.kpd6(mn).eq.-1.or.ipres(mn).eq.-1.or.
     1          iy(mn).eq.-1.or.im(mn).eq.-1.or.id(mn).eq.-1.or.
     2          ih(mn).eq.-1.or.ifh(mn).eq.-1
          lwork2=kpd5(mn).eq.kpd5(mm).and.kpd6(mn).eq.kpd6(mm).and.
     4          ipres(mn).eq.ipres(mm).and.iy(mn).eq.iy(mm).and.
     5          im(mn).eq.im(mm).and.id(mn).eq.id(mm).and.
     3          ih(mn).eq.ih(mm).and.ifh(mn).eq.ifh(mm)
          if(lwork1.and.lwork2) then
            n=nextn
          else
            n=irec
          endif
        endif
      else
        n=nextn
      endif
c
#ifndef CRAY
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
      nextn=k
      irec=nextn
c
      if(lgrib.eq.0) go to 2000
c
c  found field, start reading
c
#ifndef CRAY
              lugb4=lugb
              call rdgb(lugb4,lgrib4,lskip4,
     1                  kpds4,kgds4,ndata4,lbms,
     2                  data4,6)
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
              call rdgb(lugb,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                  data,6)
#endif
c
      if(ndata.eq.0) then
        write(6,*) ' error in rdgb'
        write(6,*) ' kpds=',kpds
        write(6,*) ' kgds=',kgds
        write(6,*) ' lgrib,lskip=',lgrib,lskip
        call abort
      endif
      nrec=nrec+1
      if(nrec.gt.maxrec) then
        write(6,*) 'error!!  number of retrieved records exceeded',
     1             ' maxrec ',maxrec
        call abort
      endif
      write(lptr,*) 'field found ... '
      write(lptr,*) 'search keys:'
      write(lptr,*) 'idvar=',kpd5(mm),' levtyp=',kpd6(mm),
     1           ' lev=',ipres(mm),' year=',iy(mm),' month=',im(mm),
     2           ' day=',id(mm),' hour=',ih(mm),
     3           ' fhour=',ifh(mm)
      write(lptr,*) 'matched keys:'
      write(lptr,*) 'idvar=',kpds(5),' levtyp=',kpds(6),
     1           ' lev=',kpds(7),' year=',kpds(8),' month=',kpds(9),
     2           ' day=',kpds(10),' hour=',kpds(11),
     3           ' fhour=',kpds(14)
      write(lptr,*) 'in file:',fngrib(nfile)(1:50)
c     call maxmin(data,ndata,1,ndata,1,1,lptr)
c
      kpd5(mm)=kpds(5)
      kpd6(mm)=kpds(6)
      ipres(mm)=kpds(7)
      io=kgds(2)
      jo=kgds(3)
      iy(mm)=kpds(8)
      im(mm)=kpds(9)
      id(mm)=kpds(10)
      ih(mm)=kpds(11)
      ifh(mm)=kpds(14)
c
      imax=io
      jmax=jo
      ijmax=imax*jmax
c
      if(imax*jmax.gt.mdata) then
        write(6,*) 'error imax*jmax too large'
        call abort
      endif
c
      call getarea(kgds,dlat,dlon,rslat,rnlat,wlon,elon,ijordr,lptr)
      write(lptr,*) 'imax,jmax,ijmax,dlon,dlat,ijordr,wlon,rnlat='
      write(lptr,*)  imax,jmax,ijmax,dlon,dlat,ijordr,wlon,rnlat
      if(abs(dlat).ne.abs(dlon)) then
        write(lptr,*) 'dlat.ne.dlon'
      endif
      call subst(data,imax,jmax,ijmax,dlon,dlat,ijordr,work)
      write(lptr,*) ' maxmin of the field'
      call maxmin(data,ijmax,1,ijmax,1,1,lptr)
c
      dlon=360./float(io)
      dlat=180./float(jo-1)
      if(kgds(1).eq.4) then
        call gau2ll(data,imax,jmax,
     1              0.,90.,dlon,dlat,datao((mm-1)*io*jo+1),io,jo,
     2              undef,lptr)
c
      else
        if(imax.ne.io.or.jmax.ne.jo) then
          dloni=360./float(imax)
          dlati=180./float(jmax-1)
          call ll2ll(data,imax,jmax,dloni,dlati,0.,90.,.true.,
     1               0.,90.,dlon,dlat,datao((mm-1)*io*jo+1),io,jo,
     2               1.e30,lptr)
        else
          do ij=1,ijmax
            datao((mm-1)*io*jo+ij)=data(ij)
          enddo
        endif
      endif
c
      fgribn(mm)=fngrib(nfile)
      lgot(mm)=.true.
c
c  if maxrec.eq.1 return (specal case)
c
      if(maxrec.eq.1) return
c
      lwork1=kpd5(mm).eq.-1.or.kpd6(mm).eq.-1.or.ipres(mm).eq.-1.or.
     1     iy(mm).eq.-1.or.im(mm).eq.-1.or.id(mm).eq.-1.or.
     2     ih(mm).eq.-1.or.ifh(mm).eq.-1
      nwild=nwild+1
      if(lwork1) go to 3000
c
      write(lptr,*) '--------------------------------------------'
      kount=kount+1
      if(kount.eq.nkeys) then
c       close(lugb)
        return
      endif
c
 2000 continue
c
      enddo
c
      mm=mm-1
      write(lptr,*) ' file:',fngrib(nfile)(1:50),' scanned for:'
      write(lptr,*) 'idvar=',kpd5(mm),' levtyp=',kpd6(mm),
     1              ' lev=',ipres(mm),' year=',iy(mm),' month=',im(mm),
     2           ' day=',id(mm),' hour=',ih(mm),
     3           ' fhour=',ifh(mm)
      write(lptr,*) 'but field was not found in this portion.'
      nfile=nfile+1
      kfile=kfile+1
      ifile=nfile
      go to 1000
c
      end
