      subroutine nrdgb(bbuf,isize,lgrib,lskip,kpds,kgds,ndata,lbms,
     1                 data,luptr)
c
c  read grib file
c  input
c    bbuf  - total grib record in array
c    isize - size of grib record
c    lgrib - length of grib record
c    lskip - bytes to skip for grib record
c  output
c    kpds(22) - unpacked product definition section
c    kgds(22) - unpacked grid definition section
c    ndata    - number of data points
c    lbms(ndata) - logical bit map
c    data(ndata) - data unpacked
c
c     parameter(llgrib=2500*1250*4)
c     character grib(lgrib)*1
      character bbuf(isize)
c     character grib(llgrib)*1
      character grib(isize)*1
      integer kpds(25),kgds(22),kptr(16)
      logical*1 lbms(*)
      real data(*)
      ndata=0
      call nbaread(bbuf,isize,lskip,lgrib,lread,grib)
      if(lread.lt.lgrib) then
        write(luptr,*) ' error in rdgb.  lread.lt.lgrib'
        call abort
      endif
      call w3fi63(grib,kpds,kgds,lbms,data,kptr,iret)
      if(iret.ne.0) then
        write(luptr,*) ' error in rdgb.  iret.ne.0 from w3fi63'
        write(luptr,*) ' iret=',iret
        call abort
      endif
      ndata=kptr(10)
      return
      end
