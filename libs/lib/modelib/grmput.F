      subroutine grmput(lugb,fngrib,f,ibms,lbms,io,jo,typ,
     1                  kpd5,kpd6,ipres,iy,im,id,ih,if1,if2,numver,
     2                  iret)
c
      save
c subprogram:    grmput      create grib file
c   prgmmr: kanamitsu          org: w/nmc23    date: 92-10-31
c
c abstract: create a grib record from a full field.
c   at present, only global latlon grids and gaussian grids
c   and regional polar projections are allowed.
c
c program history log:
c   92-10-31  iredell
c   94-05-04  juang (for gsm and rsm use)
c   96-03-18  kanamitsu
c
c usage:    call grmput(lugb,fngrib,f,ibms,lbms,io,jo,typ,
c    &                  kpd5,kpd6,ipres,iy,im,id,ih,if1,if2,numvar)
c
c   input argument list:
c     lugb     - integer unit number
c     fngrib   - character*80 output grib file name
c     f        - real (io*jo) field data to pack into grib message
c     lbms     - logical (im*jm) bitmap to use if ibms=1
c     io       - integer longitudinal dimension
c     jo       - integer latitudinal dimension
c     typ      - character*3 type of data 'avr'/'acc'/'clm' or other
c                avr=average, acc=accumulation, clm=climatology
c                if not any of those, then inatantaneous forecast/analysis
c     kpd5     - integer parameter and unit indicator
c     kpd6     - integer type of level indicator
c     il1      - integer first level value (0 for single level)
c     ipres    - integer second level value
c     iy       - integer year
c     im       - integer month
c     id       - integer day
c     ih       - integer hour
c     if1      - integer first time period
c     if2      - integer second time period (0 for single period)
c     numver   - number of samples went into the average/accumulation
c                used only when typ='clm' 
c
c   output argument list:
c     iret     - integer error code (0 for success)
c
c subprograms called:
c   assign     - file assign used in cray
c   gribit     - create grib packed field
c   wryte      - write grib record sequentially
c
c attributes:
c   language: fortran 77
c
c$$$
c
      save
c
      parameter(kpd13=1)
c
      parameter(lenpds=28,lengds=32,mxbit=16)
      parameter(maxio=360,maxjo=181)
c
      character*80 fngrib
c
      logical lbms(io,jo)
c
      dimension jdumy(5)
c
      character*3 typ
c
      character grib(30+lenpds+lengds+maxio*maxjo*(mxbit+1)/8)
c
      dimension ids(255)
c
c  id for pressure level
c
      data icen/7/,icen2/1/,igen/0/
c
      character*80 asgnstr
c
      data ifp/0/
c
      if(io*jo.gt.maxio*maxjo) then
        write(6,*) 'error io*jo too large'
        call abort
      endif
c
      if(typ(1:3).eq.'avr') then
        kpd16=3
      elseif(typ(1:3).eq.'acc') then
        kpd16=4
      elseif(typ(1:3).eq.'clm') then
        kpd16=51
      else
        kpd16=10
      endif
c
c  open output grib file
c
      if(ifp.eq.0) then
        ifp=1
        if(lugb.lt.10) then
          write(asgnstr,'(22hassign -s unblocked u:,i1,)') lugb
        else
          write(asgnstr,'(22hassign -s unblocked u:,i2,)') lugb
        endif
        call assign(asgnstr)
        open(lugb,file=fngrib,form='unformatted',err=910)
        go to 911
  910   continue
        write(6,*) ' error in opening file '
        print *,'error in opening file '
        call abort
  911   continue
        write(6,*) ' grib file:',fngrib(1:50),' opened. unit=',lugb
        call idsdef(2,ids)
      endif
c
c  write out grib messages sequentially
c
      call gribit(f,lbms,0,io,jo,mxbit,1.,
     &            lenpds,2,icen,igen,ibms,kpd5,kpd6,0,ipres,
     &            iy,im,id,ih,kpd13,if1,if2,kpd16,
     &            numavr,0,icen2,ids(kpd5),jdumy,
     &            0.,0.,0.,0.,0.,0.,0.,0.,
     &            grib,lgrib,ierr)
      call wryte(lugb,lgrib,grib)
      print *,' grib written to ',lugb,' of length ',lgrib
c
      return
      end
      subroutine gribit(f,lbm,idrt,im,jm,mxbit,colat1,
     &                  ilpds,iptv,icen,igen,ibms,ipu,itl,il1,il2,
     &                  iyr,imo,idy,ihr,iftu,ip1,ip2,itr,
     &                  ina,inm,icen2,ids,iens,
     &                  xlat1,xlon1,xlat2,xlon2,delx,dely,oritru,proj,
     &                  grib,lgrib,ierr)
c$$$  subprogram documentation block
c
c subprogram:    gribit      create grib message
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31
c
c abstract: create a grib message from a full field.
c   at present, only global latlon grids and gaussian grids
c   and regional polar projections are allowed.
c
c program history log:
c   92-10-31  iredell
c   94-05-04  juang (for gsm and rsm use)
c
c usage:    call gribit(f,lbm,idrt,im,jm,mxbit,colat1,
c    &                  ilpds,iptv,icen,igen,ibms,ipu,itl,il1,il2,
c    &                  iyr,imo,idy,ihr,iftu,ip1,ip2,itr,
c    &                  ina,inm,icen2,ids,iens,
c    &                  xlat1,xlon1,delx,dely,oritru,proj,
c    &                  grib,lgrib,ierr)
c   input argument list:
c     f        - real (im*jm) field data to pack into grib message
c     lbm      - logical (im*jm) bitmap to use if ibms=1
c     idrt     - integer data representation type
c                (0 for latlon or 4 for gaussian or 5 for polar)
c     im       - integer longitudinal dimension
c     jm       - integer latitudinal dimension
c     mxbit    - integer maximum number of bits to use (0 for no limit)
c     colat1   - real first colatitude of grid if idrt=4 (radians)
c     ilpds    - integer length of the pds (usually 28)
c     iptv     - integer parameter table version (usually 2)
c     icen     - integer forecast center (usually 7)
c     igen     - integer model generating code
c     ibms     - integer bitmap flag (0 for no bitmap)
c     ipu      - integer parameter and unit indicator
c     itl      - integer type of level indicator
c     il1      - integer first level value (0 for single level)
c     il2      - integer second level value
c     iyr      - integer year
c     imo      - integer month
c     idy      - integer day
c     ihr      - integer hour
c     iftu     - integer forecast time unit (1 for hour)
c     ip1      - integer first time period
c     ip2      - integer second time period (0 for single period)
c     itr      - integer time range indicator (10 for single period)
c     ina      - integer number included in average
c     inm      - integer number missing from average
c     icen2    - integer forecast subcenter
c                (usually 0 but 1 for reanal or 2 for ensemble)
c     ids      - integer decimal scaling
c     iens     - integer (5) ensemble extended pds values
c                (application,type,identification,product,smoothing)
c                (used only if icen2=2 and ilpds>=45)
c     xlat1    - real first point of regional latitude (radians)
c     xlon1    - real first point of regional longitude (radians)
c     xlat2    - real last  point of regional latitude (radians)
c     xlon2    - real last  point of regional longitude (radians)
c     delx     - real dx on 60n for regional (m)
c     dely     - real dy on 60n for regional (m)
c     proj     - real polar projection flag 1 for north -1 for south
c                     mercater projection 0
c     oritru   - real orientation of regional polar projection or
c                     truth for regional mercater projection
c
c   output argument list:
c     grib     - character (lgrib) grib message
c     lgrib    - integer length of grib message
c                (no more than 100+ilpds+im*jm*(mxbit+1)/8)
c     ierr     - integer error code (0 for success)
c
c subprograms called:
c   gtbits     - compute number of bits and round data appropriately
c   w3fi72     - engrib data into a grib1 message
c
c attributes:
c   language: cray fortran
c
c$$$
      integer iens(5)
      real f(im*jm)
      logical lbm(im*jm)
      character grib(*)
      parameter(nibm=360*181)
      integer ibm(nibm),ipds(100),igds(100),ibds(100)
      real fr(nibm)
      character pds(255)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  determine grid parameters
      pi=acos(-1.)
      nf=im*jm
      if(idrt.eq.0) then
        if(im.eq.144.and.jm.eq.73) then
          igrid=2
        elseif(im.eq.360.and.jm.eq.181) then
          igrid=3
        else
          igrid=255
        endif
        iresfl=128
        iscan=0
        lat1=nint(90.e3)
        lon1=0
        lati=nint(180.e3/(jm-1))
        loni=nint(360.e3/im)
        igds09=-lat1
        igds10=-loni
        igds11=lati
        igds12=loni
        igds13=0
        igds14=0
      elseif(idrt.eq.4) then
        if(im.eq.192.and.jm.eq.94) then
          igrid=98
        elseif(im.eq.384.and.jm.eq.190) then
          igrid=126
        else
          igrid=255
        endif
        iresfl=128
        iscan=0
        lat1=nint(90.e3-180.e3/pi*colat1)
        lon1=0
        lati=jm/2
        loni=nint(360.e3/im)
        igds09=-lat1
        igds10=-loni
        igds11=lati
        igds12=loni
        igds13=iscan
        igds14=0
      elseif(idrt.eq.5) then    ! polar projection
        igrid=255
        lat1=nint(180.e3/acos(-1.) * xlat1)
        lon1=nint(180.e3/acos(-1.) * xlon1)
        iresfl=0
        igds09=nint(oritru*1.e3)
        igds10=delx
        igds11=dely
        if( nint(proj).eq.1  ) igds12=0         ! north polar proj
        if( nint(proj).eq.-1 ) igds12=128       ! south polat proj
        iscan=64
        igds13=iscan
        igds14=0
      elseif(idrt.eq.1) then    ! mercater projection
        igrid=255
        lat1=nint(180.e3/acos(-1.) * xlat1)
        lon1=nint(180.e3/acos(-1.) * xlon1)
        iresfl=0
        igds09=nint(180.e3/acos(-1.) * xlat2)
        igds10=nint(180.e3/acos(-1.) * xlon2)
        igds11=delx
        igds12=dely
        igds13=nint(oritru*1.e3)
        iscan=64
        igds14=iscan
      else
        ierr=40
        return
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  reset time range parameter in case of overflow
      if(itr.ge.2.and.itr.le.5.and.ip2.ge.256) then
        jp1=ip2
        jp2=0
        jtr=10
      else
        jp1=ip1
        jp2=ip2
        jtr=itr
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  fill pds parameters
      ipds(01)=ilpds    ! length of pds
      ipds(02)=iptv     ! parameter table version id
      ipds(03)=icen     ! center id
      ipds(04)=igen     ! generating model id
      ipds(05)=igrid    ! grid id
      ipds(06)=1        ! gds flag
      ipds(07)=ibms     ! bms flag
      ipds(08)=ipu      ! parameter unit id           kpds(5)
      ipds(09)=itl      ! type of level id            kpds(6)
      ipds(10)=il1      ! level 1 or 0                kpds(7)
      ipds(11)=il2      ! level 2                     kpds(7)
      ipds(12)=iyr      ! year                        kpds(8)
      ipds(13)=imo      ! month                       kpds(9)
      ipds(14)=idy      ! day                         kpds(10)
      ipds(15)=ihr      ! hour                        kpds(11)
      ipds(16)=0        ! minute                      kpds(12)
      ipds(17)=iftu     ! forecast time unit id       kpds(13)
      ipds(18)=ip1      ! time period 1               kpds(14)
      ipds(19)=ip2      ! time period 2 or 0          kpds(15)
      ipds(20)=itr      ! time range indicator        kpds(16)
      ipds(21)=ina      ! number in average           kpds(17)
      ipds(22)=inm      ! number missing
      ipds(23)=20       ! century
      ipds(24)=icen2    ! forecast subcenter
      ipds(25)=ids      ! decimal scaling
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  fill gds and bds parameters
      igds(01)=0        ! number of vertical coords
      igds(02)=255      ! vertical coord flag
      igds(03)=idrt     ! data representation type
      igds(04)=im       ! east-west points
      igds(05)=jm       ! north-south points
      igds(06)=lat1     ! latitude of origin
      igds(07)=lon1     ! longitude of origin
      igds(08)=iresfl   ! resolution flag
      igds(09)=igds09   ! latitude of end or orientation
      igds(10)=igds10   ! longitude of end or dx in meter on 60n
      igds(11)=igds11   ! lat increment or gaussian lats or dy in meter
      igds(12)=igds12   ! longitude increment or projection
      igds(13)=igds13   ! scanning mode or lat of intercut on earth for mercater
      igds(14)=igds14   ! not used or scanning mode for mercater
      igds(15)=0     ! not used
      igds(16)=0     ! not used
      igds(17)=0     ! not used
      igds(18)=0     ! not used
      ibds(1)=0       ! bds flags
      ibds(2)=0       ! bds flags
      ibds(3)=0       ! bds flags
      ibds(4)=0       ! bds flags
      ibds(5)=0       ! bds flags
      ibds(6)=0       ! bds flags
      ibds(7)=0       ! bds flags
      ibds(8)=0       ! bds flags
      ibds(9)=0       ! bds flags
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  fill bitmap and count valid data.  reset bitmap flag if all valid.
      nbm=nf
      if(ibms.ne.0) then
        nbm=0
        do i=1,nf
          if(lbm(i)) then
            ibm(i)=1
            nbm=nbm+1
          else
            ibm(i)=0
          endif
        enddo
        if(nbm.eq.nf) ipds(7)=0
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  round data and determine number of bits
      if(nbm.eq.0) then
        do i=1,nf
          fr(i)=0.
        enddo
        nbit=0
      else
        call gtbits(ipds(7),ids,nf,ibm,f,fr,fmin,fmax,nbit)
        if(mxbit.gt.0) nbit=min(nbit,mxbit)
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  create product definition section
      call w3fi68(ipds,pds)
      if(icen2.eq.2.and.ilpds.ge.45) then
        ilast=45
        call pdsens(iens,kprob,xprob,kclust,kmembr,ilast,pds)
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  create grib message
      call w3fi72(0,fr,0,nbit,1,ipds,pds,
     &            1,255,igds,0,0,ibm,nf,ibds,
     &            nfo,grib,lgrib,ierr)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
      subroutine gtbits(ibm,ids,len,mg,g,ground,gmin,gmax,nbit)
c$$$  subprogram documentation block
c
c subprogram:    gtbits      compute number of bits and round field.
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31
c
c abstract: the number of bits required to pack a given field
c   at a particular decimal scaling is computed using the field range.
c   the field is rounded off to the decimal scaling for packing.
c   the minimum and maximum rounded field values are also returned.
c   grib bitmap masking for valid data is optionally used.
c
c program history log:
c   92-10-31  iredell
c
c usage:    call gtbits(ibm,ids,len,mg,g,gmin,gmax,nbit)
c   input argument list:
c     ibm      - integer bitmap flag (=0 for no bitmap)
c     ids      - integer decimal scaling
c                (e.g. ids=3 to round field to nearest milli-value)
c     len      - integer length of the field and bitmap
c     mg       - integer (len) bitmap if ibm=1 (0 to skip, 1 to keep)
c     g        - real (len) field
c
c   output argument list:
c     ground   - real (len) field rounded to decimal scaling
c     gmin     - real minimum valid rounded field value
c     gmax     - real maximum valid rounded field value
c     nbit     - integer number of bits to pack
c
c subprograms called:
c   isrchne  - find first value in an array not equal to target value
c
c attributes:
c   language: cray fortran
c
c$$$
      dimension mg(len),g(len),ground(len)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  round field and determine extremes where bitmap is on
      ds=10.**ids
      if(ibm.eq.0) then
        ground(1)=nint(g(1)*ds)/ds
        gmax=ground(1)
        gmin=ground(1)
        do i=2,len
          ground(i)=nint(g(i)*ds)/ds
          gmax=max(gmax,ground(i))
          gmin=min(gmin,ground(i))
        enddo
      else
        i1=isrchne(len,mg,1,0)
        if(i1.gt.0.and.i1.le.len) then
          ground(i1)=nint(g(i1)*ds)/ds
          gmax=ground(i1)
          gmin=ground(i1)
          do i=i1+1,len
            if(mg(i).ne.0) then
              ground(i)=nint(g(i)*ds)/ds
              gmax=max(gmax,ground(i))
              gmin=min(gmin,ground(i))
            endif
          enddo
        else
          gmax=0.
          gmin=0.
        endif
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  compute number of bits
      nbit=log((gmax-gmin)*ds+0.9)/log(2.)+1.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
      subroutine idsdef(iptv,ids)
c$$$  subprogram documentation block
c
c subprogram: idsdef         sets default decimal scalings
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31
c
c abstract: sets decimal scalings defaults for various parameters.
c   a decimal scaling of -3 means data is packed in kilo-si units.
c
c program history log:
c   92-10-31  iredell
c
c usage:    call idsdef(iptv,ids)
c   input arguments:
c     iptv         paramter table version (only 1 or 2 is recognized)
c   output arguments:
c     ids          integer (255) decimal scalings
c                  (unknown decimal scalings will not be set)
c
c attributes:
c   language: cray fortran
c
c$$$
      dimension ids(255)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do i=1,255
        ids(i)=0
      enddo
c
      if(iptv.eq.1.or.iptv.eq.2) then
        ids(001)=-1     ! pressure (pa)
        ids(002)=-1     ! sea-level pressure (pa)
        ids(003)=5      ! pressure tendency (pa/s)
                        !
                        !
        ids(006)=-1     ! geopotential (m2/s2)
        ids(007)=0      ! geopotential height (m)
        ids(008)=0      ! geometric height (m)
        ids(009)=0      ! standard deviation of height (m)
                        !
        ids(011)=1      ! temperature (k)
        ids(012)=1      ! virtual temperature (k)
        ids(013)=1      ! potential temperature (k)
        ids(014)=1      ! pseudo-adiabatic potential temperature (k)
        ids(015)=1      ! maximum temperature (k)
        ids(016)=1      ! minimum temperature (k)
        ids(017)=1      ! dewpoint temperature (k)
        ids(018)=1      ! dewpoint depression (k)
        ids(019)=4      ! temperature lapse rate (k/m)
        ids(020)=0      ! visibility (m)
                        ! radar spectra 1 ()
                        ! radar spectra 2 ()
                        ! radar spectra 3 ()
                        !
        ids(025)=1      ! temperature anomaly (k)
        ids(026)=-1     ! pressure anomaly (pa)
        ids(027)=0      ! geopotential height anomaly (m)
                        ! wave spectra 1 ()
                        ! wave spectra 2 ()
                        ! wave spectra 3 ()
        ids(031)=0      ! wind direction (degrees)
        ids(032)=1      ! wind speed (m/s)
        ids(033)=1      ! zonal wind (m/s)
        ids(034)=1      ! meridional wind (m/s)
        ids(035)=-4     ! streamfunction (m2/s)
        ids(036)=-4     ! velocity potential (m2/s)
        ids(037)=-1     ! montgomery stream function (m2/s2)
        ids(038)=8      ! sigma vertical velocity (1/s)
        ids(039)=3      ! pressure vertical velocity (pa/s)
        ids(040)=4      ! geometric vertical velocity (m/s)
        ids(041)=6      ! absolute vorticity (1/s)
        ids(042)=6      ! absolute divergence (1/s)
        ids(043)=6      ! relative vorticity (1/s)
        ids(044)=6      ! relative divergence (1/s)
        ids(045)=4      ! vertical u shear (1/s)
        ids(046)=4      ! vertical v shear (1/s)
        ids(047)=0      ! direction of current (degrees)
                        ! speed of current (m/s)
                        ! u of current (m/s)
                        ! v of current (m/s)
        ids(051)=4      ! specific humidity (kg/kg)
        ids(052)=0      ! relative humidity (percent)
        ids(053)=4      ! humidity mixing ratio (kg/kg)
        ids(054)=1      ! precipitable water (kg/m2)
        ids(055)=-1     ! vapor pressure (pa)
        ids(056)=-1     ! saturation deficit (pa)
        ids(057)=1      ! evaporation (kg/m2)
        ids(058)=1      ! cloud ice (kg/m2)
        ids(059)=6      ! precipitation rate (kg/m2/s)
        ids(060)=0      ! thunderstorm probability (percent)
        ids(061)=1      ! total precipitation (kg/m2)
        ids(062)=1      ! large-scale precipitation (kg/m2)
        ids(063)=1      ! convective precipitation (kg/m2)
        ids(064)=6      ! water equivalent snowfall rate (kg/m2/s)
        ids(065)=0      ! water equivalent of snow depth (kg/m2)
        ids(066)=2      ! snow depth (m)
                        ! mixed-layer depth (m)
                        ! transient thermocline depth (m)
                        ! main thermocline depth (m)
                        ! main thermocline anomaly (m)
        ids(071)=0      ! total cloud cover (percent)
        ids(072)=0      ! convective cloud cover (percent)
        ids(073)=0      ! low cloud cover (percent)
        ids(074)=0      ! middle cloud cover (percent)
        ids(075)=0      ! high cloud cover (percent)
        ids(076)=1      ! cloud water (kg/m2)
                        !
        ids(078)=1      ! convective snow (kg/m2)
        ids(079)=1      ! large scale snow (kg/m2)
        ids(080)=1      ! water temperature (k)
        ids(081)=0      ! sea-land mask ()
                        ! deviation of sea level from mean (m)
        ids(083)=5      ! roughness (m)
        ids(084)=0      ! albedo (percent)
        ids(085)=1      ! soil temperature (k)
        ids(086)=0      ! soil wetness (kg/m2)
        ids(087)=0      ! vegetation (percent)
                        ! salinity (kg/kg)
        ids(089)=4      ! density (kg/m3)
        ids(090)=1      ! runoff (kg/m2)
        ids(091)=0      ! ice concentration ()
                        ! ice thickness (m)
        ids(093)=0      ! direction of ice drift (degrees)
                        ! speed of ice drift (m/s)
                        ! u of ice drift (m/s)
                        ! v of ice drift (m/s)
                        ! ice growth (m)
                        ! ice divergence (1/s)
        ids(099)=1      ! snow melt (kg/m2)
                        ! sig height of waves and swell (m)
        ids(101)=0      ! direction of wind waves (degrees)
                        ! sig height of wind waves (m)
                        ! mean period of wind waves (s)
        ids(104)=0      ! direction of swell waves (degrees)
                        ! sig height of swell waves (m)
                        ! mean period of swell waves (s)
        ids(107)=0      ! primary wave direction (degrees)
                        ! primary wave mean period (s)
        ids(109)=0      ! secondary wave direction (degrees)
                        ! secondary wave mean period (s)
        ids(111)=0      ! net solar radiative flux at surface (w/m2)
        ids(112)=0      ! net longwave radiative flux at surface (w/m2)
        ids(113)=0      ! net solar radiative flux at top (w/m2)
        ids(114)=0      ! net longwave radiative flux at top (w/m2)
        ids(115)=0      ! net longwave radiative flux (w/m2)
        ids(116)=0      ! net solar radiative flux (w/m2)
        ids(117)=0      ! total radiative flux (w/m2)
                        !
                        !
                        !
        ids(121)=0      ! latent heat flux (w/m2)
        ids(122)=0      ! sensible heat flux (w/m2)
        ids(123)=0      ! boundary layer dissipation (w/m2)
        ids(124)=3      ! u wind stress (n/m2)
        ids(125)=3      ! v wind stress (n/m2)
                        ! wind mixing energy (j)
                        ! image data ()
        ids(128)=-1     ! mean sea-level pressure (stdatm) (pa)
        ids(129)=-1     ! mean sea-level pressure (maps) (pa)
        ids(130)=-1     ! mean sea-level pressure (eta) (pa)
        ids(131)=1      ! surface lifted index (k)
        ids(132)=1      ! best lifted index (k)
        ids(133)=1      ! k index (k)
        ids(134)=1      ! sweat index (k)
        ids(135)=10     ! horizontal moisture divergence (kg/kg/s)
        ids(136)=4      ! speed shear (1/s)
        ids(137)=5      ! 3-hr pressure tendency (pa/s)
        ids(138)=6      ! brunt-vaisala frequency squared (1/s2)
        ids(139)=11     ! potential vorticity (mass-weighted) (1/s/m)
        ids(140)=0      ! rain mask ()
        ids(141)=0      ! freezing rain mask ()
        ids(142)=0      ! ice pellets mask ()
        ids(143)=0      ! snow mask ()
                        !
                        !
                        !
                        !
                        !
                        !
                        ! covariance between v and u (m2/s2)
                        ! covariance between u and t (k*m/s)
                        ! covariance between v and t (k*m/s)
                        !
                        !
        ids(155)=0      ! ground heat flux (w/m2)
        ids(156)=0      ! convective inhibition (w/m2)
                        ! convective ape (j/kg)
                        ! turbulent ke (j/kg)
                        ! condensation pressure of lifted parcel (pa)
        ids(160)=0      ! clear sky upward solar flux (w/m2)
        ids(161)=0      ! clear sky downward solar flux (w/m2)
        ids(162)=0      ! clear sky upward longwave flux (w/m2)
        ids(163)=0      ! clear sky downward longwave flux (w/m2)
        ids(164)=0      ! cloud forcing net solar flux (w/m2)
        ids(165)=0      ! cloud forcing net longwave flux (w/m2)
        ids(166)=0      ! visible beam downward solar flux (w/m2)
        ids(167)=0      ! visible diffuse downward solar flux (w/m2)
        ids(168)=0      ! near ir beam downward solar flux (w/m2)
        ids(169)=0      ! near ir diffuse downward solar flux (w/m2)
                        !
                        !
        ids(172)=3      ! momentum flux (n/m2)
        ids(173)=0      ! mass point model surface ()
        ids(174)=0      ! velocity point model surface ()
        ids(175)=0      ! sigma layer number ()
        ids(176)=2      ! latitude (degrees)
        ids(177)=2      ! east longitude (degrees)
                        !
                        !
                        !
        ids(181)=9      ! x-gradient log pressure (1/m)
        ids(182)=9      ! y-gradient log pressure (1/m)
        ids(183)=5      ! x-gradient height (m/m)
        ids(184)=5      ! y-gradient height (m/m)
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
                        !
        ids(201)=0      ! ice-free water surcace (percent)
                        !
                        !
        ids(204)=0      ! downward solar radiative flux (w/m2)
        ids(205)=0      ! downward longwave radiative flux (w/m2)
                        !
        ids(207)=0      ! moisture availability (percent)
                        ! exchange coefficient (kg/m2/s)
        ids(209)=0      ! number of mixed layer next to sfc ()
                        !
        ids(211)=0      ! upward solar radiative flux (w/m2)
        ids(212)=0      ! upward longwave radiative flux (w/m2)
        ids(213)=0      ! non-convective cloud cover (percent)
        ids(214)=6      ! convective precipitation rate (kg/m2/s)
        ids(215)=7      ! total diabatic heating rate (k/s)
        ids(216)=7      ! total radiative heating rate (k/s)
        ids(217)=7      ! total diabatic nonradiative heating rate (k/s)
        ids(218)=2      ! precipitation index (fraction)
        ids(219)=1      ! std dev of ir t over 1x1 deg area (k)
        ids(220)=4      ! natural log of surface pressure over 1 kpa ()
                        !
        ids(222)=0      ! 5-wave geopotential height (m)
        ids(223)=1      ! plant canopy surface water (kg/m2)
                        !
                        !
                        ! blackadars mixing length (m)
                        ! asymptotic mixing length (m)
        ids(228)=1      ! potential evaporation (kg/m2)
        ids(229)=0      ! snow phase-change heat flux (w/m2)
                        !
        ids(231)=3      ! convective cloud mass flux (pa/s)
        ids(232)=0      ! downward total radiation flux (w/m2)
        ids(233)=0      ! upward total radiation flux (w/m2)
        ids(224)=1      ! baseflow-groundwater runoff (kg/m2)
        ids(225)=1      ! storm surface runoff (kg/m2)
                        !
                        !
        ids(238)=0      ! snow cover (percent)
        ids(239)=1      ! snow temperature (k)
                        !
        ids(241)=7      ! large scale condensation heating rate (k/s)
        ids(242)=7      ! deep convective heating rate (k/s)
        ids(243)=10     ! deep convective moistening rate (kg/kg/s)
        ids(244)=7      ! shallow convective heating rate (k/s)
        ids(245)=10     ! shallow convective moistening rate (kg/kg/s)
        ids(246)=7      ! vertical diffusion heating rate (kg/kg/s)
        ids(247)=7      ! vertical diffusion zonal acceleration (m/s/s)
        ids(248)=7      ! vertical diffusion merid acceleration (m/s/s)
        ids(249)=10     ! vertical diffusion moistening rate (kg/kg/s)
        ids(250)=7      ! solar radiative heating rate (k/s)
        ids(251)=7      ! longwave radiative heating rate (k/s)
                        ! drag coefficient ()
                        ! friction velocity (m/s)
                        ! richardson number ()
                        !
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
      subroutine wryte(lu,lc,c)
c$$$  subprogram documentation block
c
c subprogram:    wryte       write data out by bytes
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31
c
c abstract: efficiently write unformatted a characeter array.
c
c program history log:
c   91-10-31  mark iredell
c
c usage:    call wryte(lu,lc,c)
c
c   input argument list:
c     lu       - integer unit to which to write
c     lc       - integer number of characters or bytes to write
c     c        - characeter (lc) data to write
c
c attributes:
c   language: cray fortran
c
c$$$
      character c(lc)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      write(lu) c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
      function isrchne(n,x,incx,target)
      integer x(*), target
      j=1
      isrchne=0
      if(n.le.0) return
      if(incx.lt.0) j=1-(n-1)*incx
      do i=1,n
        if(x(j).eq.target) then
          isrchne=i
          return
        endif
        j=j+incx
      enddo
      return
      end

