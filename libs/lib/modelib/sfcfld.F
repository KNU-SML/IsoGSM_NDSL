      subroutine sfcfld(sfcftyp,iflag,
     1                  nrecs,lev,svar,maxlev)
c
c  returns number of 2-D fields per record in the surface file.
c
c  input:
c   sfcftyp: character ['osu1','osu2','noa1','ncp1','ncp2']
c   iflag  : integer.  =0 only returns nrecs and maxlev (for allocate)
c                      =l regular return
c      call with iflag=0, allocate the array and call with iflag=1
c  output:
c   nrecs:  integer  .. number of records
c   lev  :  int  array (nrecs)  .. number of 2-D fields in the record
c   svar :  char array (nrecs)  .. short variable name
c   maxlev: int      .. max level/type in this file
c   
c  note that mountain variance is removed from regional model sfc file
c  to make it consistent with global model   
c   
c 
c                                                 (ncep T382 "gfs" sfc files)
c                                                 (200506)    (200608)  (200509 T574)
c      osu1       osu2       noa1       vic1        ncp1        ncp2     ncp3
c---------------------------------------------------------------------------------
c  1   ts         ts         ts         ts          sli         sli      sli
c  2   smc(ksoil) smc(ksoil) smc(ksoil) smc(ksoil)  orog        orog     orog
c  3   sno        sno        sno        sno         ts          ts       ts
c  4   stc(ksoil) stc(ksoil) stc(ksoil) stc(nsoil)  sno         sno      sno
c  5   tg3        tg3        tg3        tg3         tg3         tg3      tg3
c  6   z0         z0         z0         z0          z0          z0       z0
c  7   cv         cv         cv         cv          albvs       albvs    albvs
c  8   cvb        cvb        cvb        cvb         albvw       albvw    albvw
c  9   cvt        cvt        cvt        cvt         albns       albns    albns
c 10   alb(kalbd) alb(kalbd) alb(kalbd) alb(kalbd)  albnw       albnw    albnw
c 11   sli        sli        sli        sli         vfrac       vfrac    vfrac
c 12   plantr     vfrac      vfrac      vfrac       canop       canop    canop
c 13   canop      canop      canop      canop       f10m        f10m     f10m
c 14   f10m       f10m       f10m       f10m        vegtyp      t2m      t2m
c 15              vegtyp     vegtyp     vegtyo      soiltyp     q2m      q2m
c 16              soiltyp    soiltyp    rt(ksoil)   albfs       vegtyp   vegtyp
c 17              albf(2)    albf(2)    albf(2)     albfw       soiltyp  soiltyp
c 18              ustar      ustar      ustar       ustar       albfs    albfs
c 19              fm         fm         fm          fm          albfw    albfw
c 20              fh         fh         fh          fh          ustart   ustart
c 21                         prcp       prcp        hice        fm       fm
c 22                         srflag     srflag      fice        fh       fh
c 23                         snodph     binf        prcp        hice     hice
c 24                         slc(ksoil) ds          srflag      fice     fice
c 25                         shdmin     dsm         snodph      tisfc    tisfc
c 26                         shdmax     ws          shdmin      prcp     prcp
c 27                         slope      cef         shdmax      srflag   srflag
c 28                         snoalb     expt(ksoil) slope       snodph   snodph
c 29                                    kst(ksoil)  snoalb      shdmin   shdmin
c 30                                    dph(ksoil)  stc010      shdmax   shdmax
c 31                                    bub(ksoil)  stc040      slope    slope
c 32                                    qrt(ksoil)  stc100      snoalb   snoalb
c 33                                    bkd(ksoil)  stc200      stc010   stc010
c 34                                    sld(ksoil)  smc010      stc040   stc040
c 35                                    wcr(ksoil)  smc040      stc100   stc100
c 36                                    wpw(ksoil)  smc100      stc200   stc200
c 37                                    smr(ksoil)  smc200      smc010   smc010
c 38                                    rmx(ksoil)  slc010      smc040   smc040
c 39                                    dphn(ksoil) slc040      smc100   smc100
c 40                                    smxn(ksoil) slc100      smc200   smc200
c 41                                    expn(ksoil) slc200      slc010   slc010
c 42                                    bubn(ksoil)             slc040   slc040
c 43                                    alpn(ksoil)             slc100   slc100
c 44                                    betn(ksoil)             slc200   slc200
c 45                                    gamn(ksoil)
c 46                                    flai
c 47                                    silz
c 48                                    snwz
c 49                                    sic(ksoil)
c 50                                    csno
c 51                                    rsno
c 52                                    tsf
c 53                                    tpk
c 54                                    sfw
c 55                                    pkw
c 56                                    lstsn
c
      implicit none
c
      character*4 sfcftyp
      integer iflag
      integer ksoil(100)
c
      integer lev(*)
      character*8 svar(*)
c
      integer nrecs,kalbd,maxlev
c
      if(sfcftyp(1:4).eq.'osu1') then
         nrecs=14 
         ksoil(1)=2
         kalbd=1
         maxlev=2
         if(iflag.eq.0) return
         lev( 1)=1
         lev( 2)=ksoil(1)
         lev( 3)=1
         lev( 4)=ksoil(1)
         lev( 5)=1
         lev( 6)=1
         lev( 7)=1
         lev( 8)=1
         lev( 9)=1
         lev(10)=kalbd
         lev(11)=1
         lev(12)=1
         lev(13)=1
         lev(14)=1
         svar( 1)='ts'
         svar( 2)='smc'
         svar( 3)='sno'
         svar( 4)='stc'
         svar( 5)='tg3'
         svar( 6)='z0'
         svar( 7)='cv'
         svar( 8)='cvb'
         svar( 9)='cvt'
         svar(10)='alb'
         svar(11)='sli'
         svar(12)='plantr'
         svar(13)='canop'
         svar(14)='f10m'
      elseif(sfcftyp(1:4).eq.'osu2') then
         nrecs=20 
         ksoil(1)=2
         kalbd=4
         maxlev=4
         if(iflag.eq.0) return
         lev( 1)=1
         lev( 2)=ksoil(1)
         lev( 3)=1
         lev( 4)=ksoil(1)
         lev( 5)=1
         lev( 6)=1
         lev( 7)=1
         lev( 8)=1
         lev( 9)=1
         lev(10)=kalbd
         lev(11)=1
         lev(12)=1
         lev(13)=1
         lev(14)=1
         lev(15)=1
         lev(16)=1
         lev(17)=2
         lev(18)=1
         lev(19)=1
         lev(20)=1
         svar( 1)='ts'
         svar( 2)='smc'
         svar( 3)='sno'
         svar( 4)='stc'
         svar( 5)='tg3'
         svar( 6)='z0'
         svar( 7)='cv'
         svar( 8)='cvb'
         svar( 9)='cvt'
         svar(10)='alb'
         svar(11)='sli'
         svar(12)='vegcov'
         svar(13)='canop'
         svar(14)='f10m'
         svar(15)='vegtyp'
         svar(16)='soiltyp'
         svar(17)='albf'
         svar(18)='ustar'
         svar(19)='fm'
         svar(20)='fh'
      elseif(sfcftyp(1:4).eq.'noa1') then
         nrecs=28 
         ksoil(1)=4
         kalbd=4
         maxlev=4
         if(iflag.eq.0) return
         lev( 1)=1
         lev( 2)=ksoil(1)
         lev( 3)=1
         lev( 4)=ksoil(1)
         lev( 5)=1
         lev( 6)=1
         lev( 7)=1
         lev( 8)=1
         lev( 9)=1
         lev(10)=kalbd
         lev(11)=1
         lev(12)=1
         lev(13)=1
         lev(14)=1
         lev(15)=1
         lev(16)=1
         lev(17)=2
         lev(18)=1
         lev(19)=1
         lev(20)=1
         lev(21)=1
         lev(22)=1
         lev(23)=1
         lev(24)=ksoil(1)
         lev(25)=1
         lev(26)=1
         lev(27)=1
         lev(28)=1
         svar( 1)='ts'
         svar( 2)='smc'
         svar( 3)='sno'
         svar( 4)='stc'
         svar( 5)='tg3'
         svar( 6)='z0'
         svar( 7)='cv'
         svar( 8)='cvb'
         svar( 9)='cvt'
         svar(10)='alb'
         svar(11)='sli'
         svar(12)='vegcov'
         svar(13)='canop'
         svar(14)='f10m'
         svar(15)='vegtyp'
         svar(16)='soiltyp'
         svar(17)='albf'
         svar(18)='ustar'
         svar(19)='fm'
         svar(20)='fh'
         svar(21)='prcp'
         svar(22)='srflag'
         svar(23)='snodph'
         svar(24)='slc'
         svar(25)='shdmin'
         svar(26)='shdmax'
         svar(27)='slope'
         svar(28)='snoalb'
      elseif(sfcftyp(1:4).eq.'vic1') then
         nrecs=56
         ksoil(1)=3
         ksoil(2)=5
         kalbd=1
         maxlev=5
         if(iflag.eq.0) return
         lev( 1)=1
         lev( 2)=ksoil(1)
         lev( 3)=1
         lev( 4)=ksoil(2)
         lev( 5)=1
         lev( 6)=1
         lev( 7)=1
         lev( 8)=1
         lev( 9)=1
         lev(10)=kalbd
         lev(11)=1
         lev(12)=1
         lev(13)=1
         lev(14)=1
         lev(15)=1
         lev(16)=ksoil(1)
         lev(17)=2
         lev(18)=1
         lev(19)=1
         lev(20)=1
         lev(21)=1
         lev(22)=1
         lev(23)=1
         lev(24)=1
         lev(25)=1
         lev(26)=1
         lev(27)=1
         lev(28)=ksoil(1)
         lev(29)=ksoil(1)
         lev(30)=ksoil(1)
         lev(31)=ksoil(1)
         lev(32)=ksoil(1)
         lev(33)=ksoil(1)
         lev(34)=ksoil(1)
         lev(35)=ksoil(1)
         lev(36)=ksoil(1)
         lev(37)=ksoil(1)
         lev(38)=ksoil(1)
         lev(39)=ksoil(2)
         lev(40)=ksoil(2)
         lev(41)=ksoil(2)
         lev(42)=ksoil(2)
         lev(43)=ksoil(2)
         lev(44)=ksoil(2)
         lev(45)=ksoil(2)
         lev(46)=1
         lev(47)=1
         lev(48)=1
         lev(49)=ksoil(1)
         lev(50)=1
         lev(51)=1
         lev(52)=1
         lev(53)=1
         lev(54)=1
         lev(55)=1
         lev(56)=1
         svar( 1)='ts'
         svar( 2)='smc'
         svar( 3)='sno'
         svar( 4)='stc'
         svar( 5)='tg3'
         svar( 6)='z0'
         svar( 7)='cv'
         svar( 8)='cvb'
         svar( 9)='cvt'
         svar(10)='alb'
         svar(11)='sli'
         svar(12)='vegcov'
         svar(13)='canop'
         svar(14)='f10m'
         svar(15)='vegtyp'
         svar(16)='root'
         svar(17)='albf'
         svar(18)='ustar'
         svar(19)='fm'
         svar(20)='fh'
         svar(21)='prcp'
         svar(22)='srflag'
         svar(23)='binf'
         svar(24)='ds'
         svar(25)='dsm'
         svar(26)='ws'
         svar(27)='cef'
         svar(28)='expt'
         svar(29)='kst'
         svar(30)='dph'
         svar(31)='bub'
         svar(32)='qrt'
         svar(33)='bkd'
         svar(34)='sld'
         svar(35)='wcr'
         svar(36)='wpw'
         svar(37)='smr'
         svar(38)='smx'
         svar(39)='dphn'
         svar(40)='smxn'
         svar(41)='expn'
         svar(42)='bubn'
         svar(43)='alpn'
         svar(44)='betn'
         svar(45)='gamn'
         svar(46)='flai'
         svar(47)='silz'
         svar(48)='snwz'
         svar(49)='sic'
         svar(50)='csno'
         svar(51)='rsno'
         svar(52)='tsf'
         svar(53)='tpk'
         svar(54)='sfw'
         svar(55)='pkw'
         svar(56)='lstsn'
      elseif(sfcftyp(1:4).eq.'ncp1') then
         nrecs=41
         ksoil(1)=4
         kalbd=4
         maxlev=4
         if(iflag.eq.0) return
         lev( 1)=1
         lev( 2)=1
         lev( 3)=1
         lev( 4)=1
         lev( 5)=1
         lev( 6)=1
         lev( 7)=1
         lev( 8)=1
         lev( 9)=1
         lev(10)=1
         lev(11)=1
         lev(12)=1
         lev(13)=1
         lev(14)=1
         lev(15)=1
         lev(16)=1
         lev(17)=1
         lev(18)=1
         lev(19)=1
         lev(20)=1
         lev(21)=1
         lev(22)=1
         lev(23)=1
         lev(24)=1
         lev(25)=1
         lev(26)=1
         lev(27)=1
         lev(28)=1
         lev(29)=1
         lev(30)=1
         lev(31)=1
         lev(32)=1
         lev(33)=1
         lev(34)=1
         lev(35)=1
         lev(36)=1
         lev(37)=1
         lev(38)=1
         lev(39)=1
         lev(40)=1
         lev(41)=1
         svar( 1)='sli'
         svar( 2)='orog'
         svar( 3)='ts'
         svar( 4)='sno'
         svar( 5)='tg3'
         svar( 6)='z0'
         svar( 7)='albvs'
         svar( 8)='albvw'
         svar( 9)='albns'
         svar(10)='albnw'
         svar(11)='vegcov'
         svar(12)='canop'
         svar(13)='f10m'
         svar(14)='vegtyp'
         svar(15)='soiltyp'
         svar(16)='albfs'
         svar(17)='albfw'
         svar(18)='ustar'
         svar(19)='fm'
         svar(20)='fh'
         svar(21)='hice'
         svar(22)='fice'
         svar(23)='prcp'
         svar(24)='srflag'
         svar(25)='snodph'
         svar(26)='shdmin'
         svar(27)='shdmax'
         svar(28)='slope'
         svar(29)='snoalb'
         svar(30)='stc010'
         svar(31)='stc040'
         svar(32)='stc100'
         svar(33)='stc200'
         svar(34)='smc010'
         svar(35)='smc040'
         svar(36)='smc100'
         svar(37)='smc200'
         svar(38)='slc010'
         svar(39)='slc040'
         svar(40)='slc100'
         svar(41)='slc200'
      elseif((sfcftyp(1:4).eq.'ncp2').or.(sfcftyp(1:4).eq.'ncp3')) then
         nrecs=44
         ksoil(1)=4
         kalbd=4
         maxlev=4
         if(iflag.eq.0) return
         lev( 1)=1
         lev( 2)=1
         lev( 3)=1
         lev( 4)=1
         lev( 5)=1
         lev( 6)=1
         lev( 7)=1
         lev( 8)=1
         lev( 9)=1
         lev(10)=1
         lev(11)=1
         lev(12)=1
         lev(13)=1
         lev(14)=1
         lev(15)=1
         lev(16)=1
         lev(17)=1
         lev(18)=1
         lev(19)=1
         lev(20)=1
         lev(21)=1
         lev(22)=1
         lev(23)=1
         lev(24)=1
         lev(25)=1
         lev(26)=1
         lev(27)=1
         lev(28)=1
         lev(29)=1
         lev(30)=1
         lev(31)=1
         lev(32)=1
         lev(33)=1
         lev(34)=1
         lev(35)=1
         lev(36)=1
         lev(37)=1
         lev(38)=1
         lev(39)=1
         lev(40)=1
         lev(41)=1
         lev(42)=1
         lev(43)=1
         lev(44)=1
         svar( 1)='sli'
         svar( 2)='orog'
         svar( 3)='ts'
         svar( 4)='sno'
         svar( 5)='tg3'
         svar( 6)='zorl'
         svar( 7)='albvs'
         svar( 8)='albvw'
         svar( 9)='albns'
         svar(10)='albnw'
         svar(11)='vegcov'
         svar(12)='canop'
         svar(13)='f10m'
         svar(14)='t2m'
         svar(15)='q2m'
         svar(16)='vegtyp'
         svar(17)='soiltyp'
         svar(18)='albfs'
         svar(19)='albfw'
         svar(20)='ustar'
         svar(21)='fm'
         svar(22)='fh'
         svar(23)='hice'
         svar(24)='fice'
         svar(25)='tisfc'
         svar(26)='prcp'
         svar(27)='srflag'
         svar(28)='snodph'
         svar(29)='shdmin'
         svar(30)='shdmax'
         svar(31)='slope'
         svar(32)='snoalb'
         svar(33)='stc010'
         svar(34)='stc040'
         svar(35)='stc100'
         svar(36)='stc200'
         svar(37)='smc010'
         svar(38)='smc040'
         svar(39)='smc100'
         svar(40)='smc200'
         svar(41)='slc010'
         svar(42)='slc040'
         svar(43)='slc100'
         svar(44)='slc200'
      else
         print *,'no such sfc file type'
         call abort
      endif
c
      return
      end
