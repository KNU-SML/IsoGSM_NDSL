      subroutine grmget(lugb,nkeys,fngrib,ifile,irec,maxrec,
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
c    maxminf - find max/min and prints
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
c
c     call maxminf(data,ndata,1,ndata,1,1,lptr)
c
      kpd5(mm)=kpds(5)
      kpd6(mm)=kpds(6)
      ipres(mm)=kpds(7)
c
      io=kgds(2)
      jo=kgds(3)
c
      imax=kgds(2)
      jmax=kgds(3)
c
      iy(mm)=kpds(8)
      im(mm)=kpds(9)
      id(mm)=kpds(10)
      ih(mm)=kpds(11)
      ifh(mm)=kpds(14)
c
c     imax=io
c     jmax=jo
c
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
      call maxminf(data,ijmax,1,ijmax,1,1,lptr)
c
      if(kgds(1).eq.4) then
				dlon=360./float(io)
				dlat=180./float(jo-1)
        call gau2ll(data,imax,jmax,
     1              0.,90.,dlon,dlat,datao((mm-1)*io*jo+1),io,jo,
     2              undef,lptr)
c
      else
        if(imax.ne.io.or.jmax.ne.jo) then
		      dloni=360./float(imax)
          dlati=180./float(jmax-1)
c
c  special consideration for gaussian input
c
          if(io.eq.128.and.jo.eq.64.or.
     1       io.eq.192.and.jo.eq.94.or.
     2       io.eq.384.and.jo.eq.190) then
             call la2ga(data,imax,jmax,dloni,dlati,0.,90.,
     1                 datao((mm-1)*io*jo+1),io,jo,.false.,dummy,dummy)
					else
           call ll2ll(data,imax,jmax,dloni,dlati,0.,90.,.true.,
     1                0.,90.,dlon,dlat,datao((mm-1)*io*jo+1),io,jo,
     2                1.e30,lptr)
          endif
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
      subroutine gau2ll(gauin,imxin,jmxin,
     1                  xlonw,xlatn,dxlon,dxlat,regout,imxout,jmxout,
     2                  undef,luptr)
c
      save
c                                                                       
c  interpolation from lat/lon grid to other lat/lon grid                
c                                                                       
      dimension gauin (imxin,jmxin)                                   
c                                                                       
      dimension regout(imxout,jmxout)                                     
      dimension gaul(500),regl(500)                                     
      dimension iindx1(1000)                                            
      dimension iindx2(1000)                                            
      dimension jindx1(500)                                             
      dimension jindx2(500)                                             
      dimension ddx(1000)                                               
      dimension ddy(500)                                                
c                                                                       
      data ifp/0/                                                       
c                                                                       
      if(ifp.ne.0) go to 111                                            
      ifp=1                                                             
c
      write(luptr,*) 'imxin=',imxin,' jmxin=',jmxin
      write(luptr,*) 'xlatn=',xlatn,' xlonw=',xlonw,
     1               ' dxlat=',dxlat,' dxlon=',dxlon
      write(luptr,*) 'imxout=',imxout,' jmxout=',jmxout
c                                                                       
      call gaulat(gaul,jmxin)                                           
c                                                                       
      do 20 j=1,jmxout                                                  
      regl(j)=xlatn-float(j-1)*dxlat                                       
   20 continue                                                          
c                                                                       
      dxin =360./float(imxin)                                          
c                                                                       
      do 30 i=1,imxout                                                  
      alamd=xlonw+float(i-1)*dxlon
      if(alamd.lt.0.) alamd=360.+alamd
      i1=alamd/dxin+1.001                                               
      if(i1.gt.imxin) i1=1                                              
      iindx1(i)=i1                                                      
      i2=i1+1                                                           
      if(i2.gt.imxin) i2=1                                              
      iindx2(i)=i2                                                      
      ddx(i)=(alamd-float(i1-1)*dxin)/dxin                              
   30 continue                                                          
c                                                                       
      do 40 j=1,jmxout                                                  
      aphi=regl(j)                                                      
      do 50 jj=1,jmxin                                                  
      if(aphi.lt.gaul(jj)) go to 50                                     
      j2=jj                                                             
      go to 42                                                          
   50 continue                                                          
      j2=jmxin                                                              
   42 continue                                                          
      if(j2.gt.2) go to 43                                              
      j1=1                                                              
      j2=2                                                              
      go to 44                                                          
   43 continue                                                          
      if(j2.le.jmxin) go to 45                                          
      j1=jmxin-1                                                        
      j2=jmxin                                                          
      go to 44                                                          
   45 continue                                                          
      j1=j2-1                                                           
   44 continue                                                          
      jindx1(j)=j1                                                      
      jindx2(j)=j2                                                      
      ddy(j)=(aphi-gaul(j1))/(gaul(j2)-gaul(j1))                        
   40 continue                                                          
c                                                                       
  111 continue                                                          
c                                                                       
c     write(luptr,*) 'iindx1'                                                  
c     write(luptr,*) (iindx1(n),n=1,imxout)                                    
c     write(luptr,*) 'iindx2'                                                  
c     write(luptr,*) (iindx2(n),n=1,imxout)                                    
c     write(luptr,*) 'jindx1'                                                  
c     write(luptr,*) (jindx1(n),n=1,jmxout)                                    
c     write(luptr,*) 'jindx2'                                                  
c     write(luptr,*) (jindx2(n),n=1,jmxout)                                    
c     write(luptr,*) 'ddy'                                                     
c     write(luptr,*) (ddy(n),n=1,jmxout)                                       
c     write(luptr,*) 'ddx'                                                     
c     write(luptr,*) (ddx(n),n=1,jmxout)                                       
c                                                                       
      do 60 j=1,jmxout                                                  
      y=ddy(j)                                                          
      j1=jindx1(j)                                                      
      j2=jindx2(j)                                                      
      do 60 i=1,imxout                                                  
      x=ddx(i)                                                          
      i1=iindx1(i)                                                      
      i2=iindx2(i)                                                      
      if(gauin(i1,j1).eq.undef.or.gauin(i2,j1).eq.undef.or.
     1   gauin(i1,j2).eq.undef.or.gauin(i2,j2).eq.undef) then
        regout(i,j)=undef
      else
        regout(i,j)=(1.-x)*(1.-y)*gauin(i1,j1)+(1.-y)*x*gauin(i2,j1)+     
     1           (1.-x)*y*gauin(i1,j2)+x*y*gauin(i2,j2)                 
      endif
   60 continue                                                          
c                                                                       
      sum1=0.                                                           
      sum2=0.                                                           
      do 70 i=1,imxin                                                   
      sum1=sum1+gauin(i,1)                                              
      sum2=sum2+gauin(i,jmxin)                                          
   70 continue                                                          
      sum1=sum1/float(imxin)                                            
      sum2=sum2/float(imxin)                                            
      do i=1,imxin                                                   
        if(gauin(i,1).eq.undef) sum1=undef
        if(gauin(i,jmxin).eq.undef) sum2=undef
      enddo
c                                                                       
      do 80 i=1,imxout                                                  
      if(abs(regl(1)).eq.90.) then
        regout(i,     1)=sum1                                             
      endif
      if(abs(regl(jmxout)).eq.90.) then
        regout(i,jmxout)=sum2                                             
      endif
   80 continue                                                          
c                                                                       
      return                                                            
      end                                                               
      subroutine getarea(kgds,dlat,dlon,rslat,rnlat,wlon,elon,ijordr,
     1                   lptr)
c
c  get area of the grib record
c
      dimension kgds(22)
      logical ijordr
c
      if(kgds(1).eq.0) then
c
c  lat/lon grid
c
        write(lptr,*) 'lat/lon grid'
        dlat=float(kgds(10))/1000.0
        dlon=float(kgds( 9))/1000.0
        f0lon=float(kgds(5))/1000.0
        f0lat=float(kgds(4))/1000.0
        kgds11=kgds(11)
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
        call abort
c
c  gnomonic projection
c
      elseif(kgds(1).eq.2) then
        write(6,*) 'gnomonic grid'
        write(6,*) 'error!! gnomonic projection not coded'
        call abort
c
c  lambert conformal
c
      elseif(kgds(1).eq.3) then
        write(6,*) 'lambert conformal'
        write(6,*) 'cannot process'
        call abort
      elseif(kgds(1).eq.4) then
c
c  gaussian grid
c
        write(lptr,*) 'gaussian grid'
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
        call abort
        return
c
c  oblique lambert conformal
c
      elseif(kgds(1).eq.13) then
        write(6,*) 'oblique lambert conformal grid'
        write(6,*) 'cannot process'
        call abort
c
c  spherical coefficient
c
      elseif(kgds(1).eq.50) then
        write(6,*) 'spherical coefficient'
        write(6,*) 'cannot process'
        call abort
        return
c
c  space view perspective (orthographic grid)
c
      elseif(kgds(1).eq.90) then
        write(6,*) 'space view perspective grid'
        write(6,*) 'cannot process'
        call abort
        return
c
c  unknown projection.  abort.
c
      else
        write(6,*) 'error!! unknown map projection'
        write(6,*) 'kgds(1)=',kgds(1)
        print *,'error!! unknown map projection'
        print *,'kgds(1)=',kgds(1)
        call abort
      endif
c
      return
      end
      subroutine subst(data,imax,jmax,ijmax,dlon,dlat,ijordr,work)
c
      logical ijordr
c
      dimension data(ijmax)
c
      dimension work(ijmax)
c
      if(.not.ijordr.or.
     1  (ijordr.and.(dlat.gt.0..or.dlon.lt.0.))) then
        if(.not.ijordr) then
          ij=0
          do j=1,jmax
            do i=1,imax
              ij=(j-1)*imax+i
              ji=(i-1)*jmax+j
              work(ij)=data(ji)
            enddo
          enddo
        else
          do j=1,jmax
            do i=1,imax
              ij=(j-1)*imax+i
              work(ij)=data(ij)
            enddo
          enddo
        endif
        do j=1,jmax
          do i=1,imax
            if(dlat.gt.0..and.dlon.gt.0.) then
              ij=imax*jmax-imax*j+i
            elseif(dlat.gt.0..and.dlon.lt.0.) then
              ij=imax*jmax-(j-1)*imax-imax+i-1
            elseif(dlat.lt.0..and.dlon.lt.0.) then
              ij=imax*(j-1)+imax-i+1
            endif
            ijo=(j-1)*imax+i
            data(ij)=work(ijo)
          enddo
        enddo
      endif
      return
      end
      subroutine maxminf(f,idim,jdim,imax,jmax,kmax,lptr)
c
      dimension f(idim,jdim,kmax)
c
      do 10 k=1,kmax
c
      fmax=f(1,1,k)
      fmin=f(1,1,k)
c
      do 20 j=1,jmax
      do 20 i=1,imax
      if(fmax.le.f(i,j,k)) then
      fmax=f(i,j,k)
      iimax=i
      jjmax=j
      endif
      if(fmin.ge.f(i,j,k)) then
      fmin=f(i,j,k)
      iimin=i
      jjmin=j
      endif
   20 continue
c
      write(lptr,100) k,fmax,iimax,jjmax,fmin,iimin,jjmin
  100 format(2x,'level=',i2,' max=',e10.4,' at i=',i5,' j=',i5,
     1                      ' min=',e10.4,' at i=',i5,' j=',i5)
c
   10 continue
c
      return
      end
      subroutine ll2ll(regin,imxin,jmxin,dloin,dlain,wlon,rnlat,lglobl,
     1                 xlonw,xlatn,dxlon,dxlat,regout,imxout,jmxout,
     2                 undef,lptr)
c
      save
c
c  interpolation from lat/lon or gaussian grid to other lat/lon grid
c
      dimension regin (imxin ,jmxin )
      dimension regout(imxout,jmxout)
c
      logical lglobl
c
      dimension rinlat(500),outlat(500)
      dimension rinlon(1000)
c
      dimension iindx1(1000)
      dimension iindx2(1000)
      dimension jindx1(500)
      dimension jindx2(500)
c
      dimension ddx(1000)
      dimension ddy(500)
c
      data ifpi,jfpi,ifpo,jfpo,rfp1,rfp2/4*0,2*0./
c
      write(lptr,*) 'imxin=',imxin,' jmxin=',jmxin
      write(lptr,*) 'dloin=',dloin,' dlain=',dlain,
     1               'wlon=',wlon,' rnlat=',rnlat,' lglobl=',lglobl
      write(lptr,*) 'xlatn=',xlatn,' xlonw=',xlonw,
     1               ' dxlat=',dxlat,' dxlon=',dxlon
      write(lptr,*) 'imxout=',imxout,' jmxout=',jmxout
c
      if(imxin.eq.1.or.jmxin.eq.1) then
      do 1 j=1,jmxout
      do 1 i=1,imxout
        regout(i,j)=0.0
    1 continue
      return
      endif
c
      if(dloin.eq.0.0.or.dlain.eq.0.0) then
        print *,'dloin or dlain is zero'
        call abort
      endif
c
      if(ifpi.eq.imxin .and.jfpi.eq.jmxin .and.
     1   ifpo.eq.imxout.and.jfpo.eq.jmxout.and.
     2   rfp1.eq.wlon.and.rfp2.eq.rnlat) go to 111
c
      ifpi=imxin
      jfpi=jmxin
      ifpo=imxout
      jfpo=jmxout
      rfp1=wlon
      rfp2=rnlat
c
c     write(lptr,*) 'dloin=',dloin
c     write(lptr,*) 'dlain=',dlain
c     write(lptr,*) 'wlon=',wlon
c     write(lptr,*) 'rnlat=',rnlat
c
      do 5 j=1,jmxin
        rinlat(j)=rnlat-float(j-1)*dlain
    5 continue
c
c     write(lptr,*) 'rinlat='
c     write(lptr,*) (rinlat(j),j=1,jmxin)
c
      do 7 j=1,jmxout
        outlat(j)=xlatn-float(j-1)*dxlat
    7 continue
c
c     write(lptr,*) 'outlat='
c     write(lptr,*) (outlat(j),j=1,jmxout)
c
      do 15 i=1,imxin
        rinlon(i)=wlon+float(i-1)*dloin
   15 continue
c
c     write(lptr,*) 'rinlon='
c     write(lptr,*) (rinlon(i),i=1,imxin)
c
c  find i-index for interplation
c
      do 30 i=1,imxout
        alamd=xlonw+float(i-1)*dxlon
        if(wlon.lt.0.0) then
          if(alamd.gt.180.0) alamd=alamd-360.0
        else
          if(alamd.lt.0.) alamd=360.+alamd
        endif
        if(.not.lglobl.and.alamd.gt.rinlon(imxin)) then
          i1=-999
          i2=-999
          ddx(i)=0.
        elseif(.not.lglobl.and.alamd.lt.rinlon(1)) then
          i1=-999
          i2=-999
          ddx(i)=0.
        else
          do 35 ii=1,imxin
            if(alamd.gt.rinlon(ii)) go to 35
            ix=ii
            go to 32
   35     continue
          if(.not.lglobl) then
            i1=imxin
            i2=imxin-1
          else
            i1=360.0/dloin+0.50
            if(i1.gt.imxin) i1=imxin
            i2=1
          endif
          go to 34
   32     continue
          if(ix.ge.2) go to 33
          if(.not.lglobl) then
            i1=2
            i2=1
          else
            i1=360.0/dloin+0.50
            if(i1.gt.imxin) i1=imxin
            i2=1
          endif
          go to 34
   33     continue
          i2=ix
          i1=i2-1
   34     continue
          denom=rinlon(i2)-rinlon(i1)
          if(denom.lt.0.0) denom=denom+360.0
          rnume=alamd-rinlon(i1)
          if(rnume.lt.0.0) rnume=rnume+360.0
          ddx(i)=rnume/denom
        endif
        iindx1(i)=i1
        iindx2(i)=i2
   30 continue
c
c  find j-index for interplation
c
      jq=1
      do 40 j=1,jmxout
        aphi=outlat(j)
        if(.not.lglobl.and.aphi.gt.rinlat(1)) then
          j1=-999
          j2=-999
          ddy(j)=0.
        elseif(.not.lglobl.and.aphi.lt.rinlat(jmxin)) then
          j1=-999
          j2=-999
          ddy(j)=0.
        else
          do 50 jj=1,jmxin
            if(aphi.lt.rinlat(jj)) go to 50
            jq=jj
            go to 42
   50     continue
          j1=jmxin
          j2=jmxin
          go to 44
   42     continue
          if(jq.ge.2) go to 43
          j1=1
          j2=1
          go to 44
   43     continue
          j2=jq
          j1=jq-1
   44     continue
          if(j2.ne.j1) then
             ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
          else
            if(j1.eq.1) then
              if(abs(90.0-rinlat(j1)).gt.0.0010) then
                ddy(j)=(aphi-rinlat(j1))/(90.0-rinlat(j1))
              else
                ddy(j)=0.0
              endif
            endif
            if(j1.eq.jmxin) then
              if(abs(-90.0-rinlat(j1)).gt.0.001) then
                ddy(j)=(aphi-rinlat(j1))/(-90.0-rinlat(j1))
              else
                ddy(j)=0.0
              endif
            endif
          endif
        endif
        jindx1(j)=j1
        jindx2(j)=j2
   40 continue
c
c     write(lptr,*) 'iindx1'
c     write(lptr,*) (iindx1(n),n=1,imxout)
c     write(lptr,*) 'iindx2'
c     write(lptr,*) (iindx2(n),n=1,imxout)
c     write(lptr,*) 'jindx1'
c     write(lptr,*) (jindx1(n),n=1,jmxout)
c     write(lptr,*) 'jindx2'
c     write(lptr,*) (jindx2(n),n=1,jmxout)
c     write(lptr,*) 'ddy'
c     write(lptr,*) (ddy(n),n=1,jmxout)
c     write(lptr,*) 'ddx'
c     write(lptr,*) (ddx(n),n=1,jmxout)
c
  111 continue
c
      if(lglobl) then
        sum1=0.0
        sum2=0.0
        do 80 i=1,imxin
          sum1=sum1+regin(i,1)
          sum2=sum2+regin(i,jmxin)
   80   continue
        sumn=sum1/float(imxin)
        sums=sum2/float(imxin)
      endif
c
c  quasi-bilinear interpolation
c
      do 70 j=1,jmxout
        y=ddy(j)
        j1=jindx1(j)
        j2=jindx2(j)
        do 70 i=1,imxout
          x=ddx(i)
          i1=iindx1(i)
          i2=iindx2(i)
          if(i1.eq.-999.or.i2.eq.-999.or.j1.eq.-999.or.j2.eq.-999) then
           	 regout(i,j)=undef
          elseif(regin(i1,j1).eq.undef.or.regin(i2,j1).eq.undef.or.
     1           regin(i1,j2).eq.undef.or.regin(i2,j2).eq.undef) then
           	 regout(i,j)=undef
          else
            wi1j1=(1.0-x)*(1.0-y)
            wi2j1=     x *(1.0-y)
            wi1j2=(1.0-x)*      y
            wi2j2=     x *      y
            wsum  =wi1j1+wi2j1+wi1j2+wi2j2
            wsumiv=1./wsum
            if(j1.ne.j2) then
              regout(i,j)=(wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+
     1                     wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))*wsumiv
            else
              if(j1.eq.1) then
                regout(i,j)=(wi1j1*sumn        +wi2j1*sumn        +
     1                     wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))*wsumiv
              endif
              if(j1.eq.jmxin) then
                regout(i,j)=(wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+
     1                       wi1j2*sums        +wi2j2*sums        )*wsumiv
              endif
            endif
          endif
   70 continue
c
      return
      end
      subroutine gaulat(gaul,k)                                         
c                                                                       
      dimension a(500)                                                  
      dimension gaul(1)                                                 
c                                                                       
      esp=1.e-6                                                       
      c=(1.e0-(2.e0/3.14159265358979e0)**2)*0.25e0                      
      fk=k                                                              
      kk=k/2                                                            
      call bsslz1(a,kk)                                                 
      do 30 is=1,kk                                                     
      xz=cos(a(is)/sqrt((fk+0.5e0)**2+c))                               
      iter=0                                                            
   10 pkm2=1.e0                                                         
      pkm1=xz                                                           
      iter=iter+1                                                       
      if(iter.gt.10) go to 70                                           
      do 20 n=2,k                                                       
      fn=n                                                              
      pk=((2.e0*fn-1.e0)*xz*pkm1-(fn-1.e0)*pkm2)/fn                     
      pkm2=pkm1                                                         
   20 pkm1=pk                                                           
      pkm1=pkm2                                                         
      pkmrk=(fk*(pkm1-xz*pk))/(1.e0-xz**2)                              
      sp=pk/pkmrk                                                       
      xz=xz-sp                                                          
      avsp=abs(sp)                                                      
      if(avsp.gt.esp) go to 10                                          
      a(is)=xz                                                          
   30 continue                                                          
      if(k.eq.kk*2) go to 50                                            
      a(kk+1)=0.e0                                                      
      pk=2.e0/fk**2                                                     
      do 40 n=2,k,2                                                     
      fn=n                                                              
   40 pk=pk*fn**2/(fn-1.e0)**2                                          
   50 continue                                                          
      do 60 n=1,kk                                                      
      l=k+1-n                                                           
      a(l)=-a(n)                                                        
   60 continue                                                          
c                                                                       
      radi=180./(4.*atan(1.))                                           
      do 211 n=1,k                                                      
c     gaul(n)=90.-acos(a(n))*radi                                           
      gaul(n)=acos(a(n))*radi                                           
  211 continue                                                          
c
c     print *,'gaussian lat (deg) for jmax=',k 
c     print *,(gaul(n),n=1,k)                                           
c                                                                       
      return                                                            
   70 write(6,6000)                                                     
 6000 format(//5x,14herror in gauaw//)                                  
      stop                                                              
      end                                                               
      subroutine bsslz1(bes,n)                                          
c                                                                       
      dimension bes(n)                                                  
      dimension bz(50)                                                  
c                                                                       
      data pi/3.14159265358979e0/                                       
      data bz         / 2.4048255577e0, 5.5200781103e0,                 
     $  8.6537279129e0,11.7915344391e0,14.9309177086e0,18.0710639679e0, 
     $ 21.2116366299e0,24.3524715308e0,27.4934791320e0,30.6346064684e0, 
     $ 33.7758202136e0,36.9170983537e0,40.0584257646e0,43.1997917132e0, 
     $ 46.3411883717e0,49.4826098974e0,52.6240518411e0,55.7655107550e0, 
     $ 58.9069839261e0,62.0484691902e0,65.1899648002e0,68.3314693299e0, 
     $ 71.4729816036e0,74.6145006437e0,77.7560256304e0,80.8975558711e0, 
     $ 84.0390907769e0,87.1806298436e0,90.3221726372e0,93.4637187819e0, 
     $ 96.6052679510e0,99.7468198587e0,102.888374254e0,106.029930916e0, 
     $ 109.171489649e0,112.313050280e0,115.454612653e0,118.596176630e0, 
     $ 121.737742088e0,124.879308913e0,128.020877005e0,131.162446275e0, 
     $ 134.304016638e0,137.445588020e0,140.587160352e0,143.728733573e0, 
     $ 146.870307625e0,150.011882457e0,153.153458019e0,156.295034268e0/ 
      nn=n                                                              
      if(n.le.50) go to 12                                              
      bes(50)=bz(50)                                                    
      do 5 j=51,n                                                       
    5 bes(j)=bes(j-1)+pi                                                
      nn=49                                                             
   12 do 15 j=1,nn                                                      
   15 bes(j)=bz(j)                                                      
      return                                                            
      end                                                               
      subroutine la2ga(regin ,imxin ,jmxin , dloin,dlain, rlon, rlat,
     1                 gauout,imxout,jmxout,lmask,rslmsk,slmask)
c
c  interpolation from lat/lon or gaussian grid to other lat/lon grid
c
      save
c
      dimension regin (imxin ,jmxin )
      dimension gauout(imxout,jmxout)
c
      dimension rslmsk(imxin,jmxin)
      dimension slmask(imxout,jmxout)
c
      dimension gaulo (500)
      dimension gauli (500)
      dimension rinlat(500),outlat(500)
      dimension rinlon(1000)
c
      dimension iindx1(1000)
      dimension iindx2(1000)
      dimension jindx1(500)
      dimension jindx2(500)
c
      dimension ddx(1000)
      dimension ddy(500)
c
      logical lmask
c
      data ifpi,jfpi,ifpo,jfpo,rfp1,rfp2/4*0,2*0./
c
      if(imxin.eq.1.or.jmxin.eq.1) then
        do j=1,jmxout
           do i=1,imxout
            gauout(i,j)=0.
           enddo
        enddo
        return
      endif
c
      if(dloin.eq.0..or.dlain.eq.0.) then
        print *,'dloin or dlain is zero .... check data cards'
      endif
c
      if(ifpi.eq.imxin .and.jfpi.eq.jmxin .and.
     1   ifpo.eq.imxout.and.jfpo.eq.jmxout.and.
     2   rfp1.eq.rlon.and.rfp2.eq.rlat) go to 111
c
      ifpi=imxin
      jfpi=jmxin
      ifpo=imxout
      jfpo=jmxout
      rfp1=rlon
      rfp2=rlat
c
c     print *,'dloin=',dloin
c     print *,'dlain=',dlain
c     print *,'rlon=',rlon
c     print *,'rlat=',rlat
 
      do j=1,jmxin
        if(rlat.gt.0.) then
          rinlat(j)=rlat-float(j-1)*dlain
        else
          rinlat(j)=rlat+float(j-1)*dlain
        endif
      enddo
c
c     print *,'rinlat='
c     print *,(rinlat(j),j=1,jmxin)
c
c    compute gaussian latitude for output grid
c
      call gaulat(gaulo,jmxout)
c
      do j=1,jmxout
        outlat(j)=90.-gaulo(j)
      enddo
c
c     print *,'outlat='
c     print *,(outlat(j),j=1,jmxout)
c
      do i=1,imxin
        rinlon(i)=rlon+float(i-1)*dloin
      enddo
c
c     print *,'rinlon='
c     print *,(rinlon(i),i=1,imxin)
c
c  find i-index for interplation
c
      do 30 i=1,imxout
      alamd=float(i-1)*360./float(imxout)
      if(rlon.lt.0.) then
      if(alamd.gt.180.) alamd=alamd-360.
      endif
      do 35 ii=1,imxin
      if(alamd.gt.rinlon(ii)) go to 35
      ix=ii
      go to 32
   35 continue
      i1=360./dloin+0.5
      i2=1
      go to 34
   32 continue
      if(ix.ge.2) go to 33
      i1=360./dloin+0.5
      i2=1
      go to 34
   33 continue
      i2=ix
      i1=i2-1
   34 continue
      iindx1(i)=i1
      iindx2(i)=i2
      denom=rinlon(i2)-rinlon(i1)
      if(denom.lt.0.) denom=denom+360.
      rnume=alamd-rinlon(i1)
      if(rnume.lt.0.) rnume=rnume+360.
      ddx(i)=rnume/denom
   30 continue
c
c  find j-index for interplation
c
      jq=1
      do 40 j=1,jmxout
      aphi=outlat(j)
      do 50 jj=1,jmxin
      jx=jj
      if(rlat.lt.0.) jx=jmxin-jj+1
      if(aphi.lt.rinlat(jx)) go to 50
      jq=jx
      go to 42
   50 continue
      if(rlat.gt.0.) then
        j1=jmxin
        j2=jmxin
      else
        j1=1
        j2=1
      endif
      go to 44
   42 continue
      if(rlat.gt.0.) then
         if(jq.ge.2) go to 43
         j1=1
         j2=1
      else
         if(jq.lt.jmxin) go to 43
         j1=jmxin
         j2=jmxin
      endif
      go to 44
   43 continue
      if(rlat.gt.0.) then
      j2=jq
      j1=jq-1
      else
      j1=jq
      j2=jq+1
      endif
   44 continue
      jindx1(j)=j1
      jindx2(j)=j2
      if(j2.ne.j1) then
         ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
      else
      if(j1.eq.1.and.rlat.gt.0..or.j1.eq.jmxin.and.rlat.lt.0.) then
            if(abs(90.-rinlat(j1)).gt.0.001) then
               ddy(j)=(aphi-rinlat(j1))/(90.-rinlat(j1))
            else
               ddy(j)=0.0
            endif
      endif
      if(j1.eq.1.and.rlat.lt.0..or.j1.eq.jmxin.and.rlat.gt.0.) then
            if(abs(-90.-rinlat(j1)).gt.0.001) then
               ddy(j)=(aphi-rinlat(j1))/(-90.-rinlat(j1))
            else
               ddy(j)=0.0
            endif
         endif
      endif
   40 continue
c
c     print *,'la2ga'
c     print *,'iindx1'
c     print *,(iindx1(n),n=1,imxout)
c     print *,'iindx2'
c     print *,(iindx2(n),n=1,imxout)
c     print *,'jindx1'
c     print *,(jindx1(n),n=1,jmxout)
c     print *,'jindx2'
c     print *,(jindx2(n),n=1,jmxout)
c     print *,'ddy'
c     print *,(ddy(n),n=1,jmxout)
c     print *,'ddx'
c     print *,(ddx(n),n=1,jmxout)
c
  111 continue
c
      sum1=0.
      sum2=0.
      wei1=0.
      wei2=0.
      do i=1,imxin
        sum1=sum1+regin(i,1) * rslmsk(i,1)
        sum2=sum2+regin(i,jmxin) * rslmsk(i,jmxin)
        wei1=wei1+rslmsk(i,1)
        wei2=wei2+rslmsk(i,jmxin)
      enddo
      if(rlat.gt.0.) then
        if(wei1.gt.0.) then
          sumn=sum1/wei1
        else
          sumn=0.
        endif
        if(wei2.gt.0.) then
          sums=sum2/wei2
        else
          sums=0.
        endif
      else
        if(wei1.gt.0.) then
          sums=sum1/wei1
        else
          sums=0.
        endif
        if(wei2.gt.0.) then
          sumn=sum2/wei2
        else
          sumn=0.
        endif
      endif
c
c  quasi-bilinear interpolation
c
      ifill=0
      do 70 j=1,jmxout
      y=ddy(j)
      j1=jindx1(j)
      j2=jindx2(j)
      do 70 i=1,imxout
      x=ddx(i)
      i1=iindx1(i)
      i2=iindx2(i)
c
      if(lmask) then
        if(slmask(i,j).eq.rslmsk(i1,j1).and.
     1     slmask(i,j).eq.rslmsk(i2,j1).and.
     2     slmask(i,j).eq.rslmsk(i1,j2).and.
     3     slmask(i,j).eq.rslmsk(i2,j2)) then
          wi1j1=(1.-x)*(1.-y)  
          wi2j1=    x *(1.-y) 
          wi1j2=(1.-x)*      y 
          wi2j2=    x *      y
        elseif(slmask(i,j).eq.1.) then
          wi1j1=(1.-x)*(1.-y)  *rslmsk(i1,j1)
          wi2j1=    x *(1.-y)  *rslmsk(i2,j1)
          wi1j2=(1.-x)*      y *rslmsk(i1,j2)
          wi2j2=    x *      y *rslmsk(i2,j2)
        elseif(slmask(i,j).eq.0.) then
          wi1j1=(1.-x)*(1.-y)  *(1.-rslmsk(i1,j1))
          wi2j1=    x *(1.-y)  *(1.-rslmsk(i2,j1))
          wi1j2=(1.-x)*      y *(1.-rslmsk(i1,j2))
          wi2j2=    x *      y *(1.-rslmsk(i2,j2))
        endif
      else
        wi1j1=(1.-x)*(1.-y)  
        wi2j1=    x *(1.-y) 
        wi1j2=(1.-x)*      y 
        wi2j2=    x *      y
      endif
c
      wsum  =wi1j1+wi2j1+wi1j2+wi2j2
      if(wsum.ne.0.) then
        wsumiv = 1./wsum
c
        if(j1.ne.j2) then
            gauout(i,j)=(wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+
     1                   wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))*wsumiv
        else
           if(j1.eq.1.and.rlat.gt.0..or.j1.eq.jmxin.and.rlat.lt.0.) then
            gauout(i,j)=(wi1j1*sumn        +wi2j1*sumn        +
     1                   wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))*wsumiv
           endif
           if(j1.eq.1.and.rlat.lt.0..or.j1.eq.jmxin.and.rlat.gt.0.) then
            gauout(i,j)=(wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+
     1                   wi1j2*sums        +wi2j2*sums        )*wsumiv
           endif
        endif
      else
        if(.not.lmask) then
          write(6,*) ' wrong logic in interpolation'
          call abort
        endif
        ifill=ifill+1
        if(ifill.le.2) then
          write(6,*) 'i1,i2,j1,j2=',i1,i2,j1,j2
          write(6,*) 'rslmsk=',rslmsk(i1,j1),rslmsk(i1,j2),
     1                         rslmsk(i2,j1),rslmsk(i2,j2)
          write(6,*) 'i,j=',i,j,' slmask(i,j)=',slmask(i,j)
        endif
        do jx=j1,jmxin
          do ix=i1,imxin
            if((slmask(i,j).eq.1..and.
     1          slmask(i,j).eq.rslmsk(ix,jx)).or.
     2         (slmask(i,j).eq.0..and.
     3          slmask(i,j).eq.rslmsk(ix,jx))) then
              gauout(i,j)=regin(ix,jx)
              go to 71
            endif
          enddo
          do ix=i1,1,-1
            if((slmask(i,j).eq.1..and.
     1          slmask(i,j).eq.rslmsk(ix,jx)).or.
     2         (slmask(i,j).eq.0..and.
     3          slmask(i,j).eq.rslmsk(ix,jx))) then
              gauout(i,j)=regin(ix,jx)
              go to 71
            endif
          enddo
        enddo
        do jx=j1,1,-1
          do ix=i1,imxin
            if((slmask(i,j).eq.1..and.
     1          slmask(i,j).eq.rslmsk(ix,jx)).or.
     2         (slmask(i,j).eq.0..and.
     3          slmask(i,j).eq.rslmsk(ix,jx))) then
              gauout(i,j)=regin(ix,jx)
              go to 71
            endif
          enddo
          do ix=i1,1,-1
            if((slmask(i,j).eq.1..and.
     1          slmask(i,j).eq.rslmsk(ix,jx)).or.
     2         (slmask(i,j).eq.0..and.
     3          slmask(i,j).eq.rslmsk(ix,jx))) then
              gauout(i,j)=regin(ix,jx)
              go to 71
            endif
          enddo
        enddo
        write(6,*) ' error!!! no filling value found in la2ga'
        call abort
      endif
c
   71 continue
   70 continue
c
      if(ifill.gt.1) then
        write(6,*) ' unable to interpolate.  filled with nearest',
     1             ' point value at ',ifill,' points'
      endif
c
      print *,'max min from la2ga'
      call maxminf(gauout,imxout,jmxout,imxout,jmxout,1,6)
c
      return
      end
