      subroutine sunprm(ko,po,lpcl,lppr,
     &                  ipusun,itlsun,ip1sun,ip2sun,kslp,kli)
c$$$  subprogram documentation block
c
c subprogram:    sunprm      set parameters for sundry fields
c   prgmmr: iredell          org: w/nmc23    date: 92-10-31
c
c abstract: sets parameters for the sundry fields.
c   parameters returned are parameter indicator, level type indicator
c   and two level numbers.
c   the current nsun=29 sundry fields are:
c     1) surface pressure
c     2) surface pressure tendency
c     3) column precipitable water
c     4) column relative humidity
c     5) tropopause temperature
c     6) tropopause pressure
c     7) tropopause zonal wind
c     8) tropopause meridional wind
c     9) tropopause vertical wind speed shear
c    10) surface lifted index
c    11) best lifted index
c    12) maximum wind level temperature
c    13) maximum wind level pressure
c    14) maximum wind level zonal wind
c    15) maximum wind level meridional wind
c    16) surface orography
c    17) sea level pressure
c    18) relative humidity in sigma range (0.44,1.00)
c    19) relative humidity in sigma range (0.72,0.94)
c    20) relative humidity in sigma range (0.44,0.72)
c    21) potential temperature at sigma 0.9950
c    22) temperature at sigma 0.9950
c    23) pressure vertical velocity at sigma 0.9950
c    24) relative humidity at sigma 0.9950
c    25) zonal wind at sigma 0.9950
c    26) meridional wind at sigma 0.9950
c    27) specific humidity at sigma 0.9950
c    28) total cloud water
c    29) total precipitation water
c
c program history log:
c   92-10-31  iredell
c   94-08-01  iredell
c
c usage:    call sunprm(ko,po,ipusun,itlsun,ip1sun,ip2sun,kslp,kli)
c
c   input argument list:
c     ko       - integer number of pressure levels
c     po       - real (ko) pressure in millibars
c
c   output argument list:
c     ipusun   - integer (nsun) parameter indicators
c     itlsun   - integer (nsun) level type indicators
c     ip1sun   - integer (nsun) first level numbers
c     ip2sun   - integer (nsun) second level numbers
c     kslp     - integer (2) relevant pressure levels for slp
c     kli      - integer relevant pressure level for lifted index
c
c subprograms called:
c   isrcheq  - find first value in an array equal to target value
c
c attributes:
c   language: cray fortran
c
c$$$
      dimension po(ko)
      logical lpcl,lppr,lpo3
      parameter(nsun=29)
      parameter(nps  = 1,npst = 2,ntpw = 3,ntrh = 4,ntpt = 5,
     &          ntpp = 6,ntpu = 7,ntpv = 8,ntpsh= 9,nsli =10,
     &          nbli =11,nmwt =12,nmwp =13,nmwu =14,nmwv =15,
     &          nzs  =16,nslp =17,nrh1 =18,nrh2 =19,nrh3 =20,
     &          ns1th=21,ns1t =22,ns1o =23,ns1r =24,ns1u =25,
     &          ns1v =26,ns1q =27,ntcl =28,ntpr =29)
      dimension ipusun(nsun),itlsun(nsun)
      dimension ip1sun(nsun),ip2sun(nsun)
      dimension kslp(2)
      dimension ipudef(nsun),itldef(nsun)
      dimension ip1def(nsun),ip2def(nsun)
      dimension pslp(2)
      data ipudef/001,003,054,052,011,001,033,034,136,131,
     &            132,011,001,033,034,007,002,052,052,052,
     &            013,011,039,052,033,034,51,76,77/
      data itldef/001,001,200,200,007,007,007,007,007,001,
     &            001,006,006,006,006,001,102,108,108,108,
     &            107,107,107,107,107,107,107,200,200/
      data ip1def/000,000,000,000,000,000,000,000,000,000,
     &            000,000,000,000,000,000,000,044,072,044,
     &            00000,00000,00000,00000,00000,00000,00000,
     &            000,000/
      data ip2def/000,000,000,000,000,000,000,000,000,000,
     &            000,000,000,000,000,000,000,100,094,072,
     &            09950,09950,09950,09950,09950,09950,09950,
     &            000,000/
      data pslp/1000.,500./,pli/500./
c 
      ipusun=ipudef
      itlsun=itldef
      ip1sun=ip1def
      ip2sun=ip2def
      kslp(1)=mod(isrcheq(ko,po,1,pslp(1)),ko+1)
      kslp(2)=mod(isrcheq(ko,po,1,pslp(2)),ko+1)
      kli=mod(isrcheq(ko,po,1,pli),ko+1)
      if(kslp(1).eq.0.or.kslp(2).eq.0) ipusun(nslp)=0
      if(kli.eq.0) ipusun(nsli)=0
      if(kli.eq.0) ipusun(nbli)=0
      if(.not.lpcl) ipusun(ntcl)=0
      if(.not.lppr) ipusun(ntpr)=0
c
c  do not produce low-layer rh because of a difficulty in grib 
c
      ipusun(nrh1)=0
      ipusun(nrh2)=0
      ipusun(nrh3)=0
c
      return
      end
