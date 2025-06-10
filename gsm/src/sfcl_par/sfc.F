#define PI pi
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define JWAV1 jwav1
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVHP levhp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LEVSP levsp
#define LLWAVP llwavp
#define LNGRD lngrd
#define LNGRDP lngrdp
#define LNWAV lnwav
#define LNWAVP lnwavp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NSOIL nsoil
#define NTOTAL ntotal
#define NTRAC ntrac
#include <define.h>
#ifdef ADJSOILW
      subroutine sfc(lugb,iy,im,id,ih,fh,ifcst,
     1             intobsp,mdlpcp,mdlrnf)
#else
      subroutine sfc(lugb,iy,im,id,ih,fh,ifcst)
#endif

c
c  surface program version 3 
c    to be called inside forecast program
c    complete clean-up.
c    supposed to be easier to add additional variables for new 
c    land models
c
c  this program reads climatology and analysis from grib files
c  interpolates in time and space, merge the two, and then
c  merges with predicted sfc field array in common block
c  the program minimizes the grib file read
c
c  lugb  integer unit number
c  iy .. integer year of initial time
c  im .. integer month of initial time
c  id .. integer day of initial time
c  ih .. integer hour  of initial time
c  fh .. real forecast hour
c
      implicit none
#include <paramodel.h>
#ifndef RSM
#define LONF2F lonf2_
#define LATG2F latg2_
#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif /* MP */
#else
#define LONF2F igrd12_
#define LATG2F jgrd12_
#ifdef MP
#define LONF2S igrd12p_
#define LATG2S jgrd12p_
#else
#define LONF2S igrd12_
#define LATG2S jgrd12_
#endif
#endif

#ifdef MP
#ifndef RSM
#include <commpi.h>
#else
#include <comrmpi.h>
#endif
#endif
      integer idim,jdim,ijdim
      parameter(idim=LONF2F/2,jdim=LATG2F*2)
      parameter(ijdim=idim*jdim)
      real, allocatable :: obsst0(:),clsst0(:)
      real factor
#ifdef NFDR
#ifdef RSM
      real, allocatable :: flat(:),flon(:)
      real dlamda0
#endif
#include <comnfdr.h>
#ifdef RSM
#include <rscomloc.h>
#endif
#endif
#ifdef DBGNFDR
      real*4 st2m2pm(ijdim)
      real*4 srh2m2pm(ijdim)
      real*4 sppt24(ijdim)
      real*4 scloud2pm(ijdim)
      real*4 swind2pm(ijdim)
      real*4 sprcpdur(ijdim)
      real*4 stmax24(ijdim)
      real*4 stmin24(ijdim)
      real*4 srhmax(ijdim)
      real*4 srhmin(ijdim)
#endif
c
#include <vargrb.h>
#include <varsfc.h>
c
      real orog(ijdim), slmask(ijdim)
      save orog,slmask
      real glacier(ijdim),amxice(ijdim)
      save glacier,amxice
      real grbfld(ijdim,numsfcs)
c
c  ifcst=0 when called from sfc0
c  ifcst=1 when called from fcst or rfcst
c
      integer ifcst
#ifdef MP
c
      real grbfldp(LONF2S*LATG2S,numsfcs)
      real orogp(LONF2S*LATG2S),slmaskp(LONF2S*LATG2S)
      real glacierp(LONF2S*LATG2S),amxicep(LONF2S*LATG2S)
#endif
c
      real rad
      parameter(rad=pi_/180.)
c
      integer nf
c
c  namelist variables
c
c  fnmskg .. character file namr of grib high resolution global 
c            lat/lon land ocean mask
c  fnorog .. character file name of binary orography (mtn output) 
c  fnmask .. character file name of binary land sea mask (mtn output) 
c  fngt0  .. character array, file name of grib sfc file for ft=0.
c  fngfc  .. character array, file name of grib sfc file for ft>0.
c  lanomp .. logical.  Controlls SST anomaly persistency.
c  
c
c default file name for initial and forecast
c
      character*128 fnmskg,fnorog,fnmask
      data fnmskg/'clim.maskh.grib          '/
      data fnorog/'orogrd.smth              '/
      data fnmask/'slmsk                    '/
c
c  climatological sst file for anomaly persistency
c
      character*128 fnobsst0
      data fnobsst0/'          '/
c
c   initial sfc fields.  blank indicates not available or simply use
c   what is available in initial sfc file.
c
      character*128 fnt0(numgrbs)
      data fnt0(igla)/'clim.glacier.grib    '/
      data fnt0(imxi)/'clim.maxice.grib     '/
      data fnt0(imsk)/'                     '/ ! fnmask used
      data fnt0(iais)/'clim.ice.grib        '/
      data fnt0(iscv)/'                     '/ ! not used
      data fnt0(isno)/'clim.snow.grib       '/ ! to be overwritten by anl
      data fnt0(itsf)/'clim.sst.grib        '/ ! to be overwritten by anl
      data fnt0(iab4)/'clim.yhalbedo.grib   '/
      data fnt0(iab1)/'clim.sibalbedo.grib  '/
      data fnt0(ialf)/'clim.yhalbedo.grib   '/
      data fnt0(iso2)/'                     '/ ! not available
      data fnt0(isn1)/'                     '/ ! not available
      data fnt0(izor)/'clim.sibrough.grib   '/
      data fnt0(iplr)/'clim.sibresis.grib   '/
      data fnt0(itg3)/'clim.tg3.grib        '/
      data fnt0(ito2)/'                     '/ ! not available
      data fnt0(itn1)/'                     '/ ! not available
      data fnt0(iveg)/'clim.vegfrac40m.grib '/
      data fnt0(ivet)/'clim.vegtype40m.grib '/
      data fnt0(isot)/'clim.soiltype40m.grib'/
      data fnt0(islo)/'clim.slptyp.grib     '/
      data fnt0(isna)/'clim.snoalb.grib     '/
      data fnt0(ismx)/'clim.shdmax.grib     '/
      data fnt0(ismn)/'clim.shdmin.grib     '/
#ifdef VIC
      data fnt0(ivgv)/'glob.veg1typ.vegcov.grib'/
      data fnt0(ivtv)/'glob.veg1typ.vegtyp.grib'/
      data fnt0(ivrt)/'glob.veg1typ.vegrt.grib'/
      data fnt0(ibif)/'glob.soil.binf.grib  '/
      data fnt0(ids )/'glob.soil.ds.grib    '/
      data fnt0(idsm)/'glob.soil.dsmax.grib '/
      data fnt0(iws )/'glob.soil.ws.grib    '/
      data fnt0(icef)/'glob.soil.cef.grib   '/
      data fnt0(iexp)/'glob.soil.expt.grib  '/
      data fnt0(ikst)/'glob.soil.ksat.grib  '/
      data fnt0(idph)/'glob.soil.dph.grib   '/
      data fnt0(ibub)/'glob.soil.bble.grib  '/
      data fnt0(iqrt)/'glob.soil.qurtz.grib '/
      data fnt0(ibkd)/'glob.soil.blkd.grib  '/
      data fnt0(isld)/'glob.soil.sden.grib  '/
      data fnt0(iwcr)/'glob.soil.wcr.grib   '/
      data fnt0(iwpw)/'glob.soil.wpwp.grib  '/
      data fnt0(ismr)/'glob.soil.smr.grib   '/
      data fnt0(islz)/'glob.soil.silz.grib  '/
      data fnt0(isnz)/'glob.soil.snwz.grib  '/
      data fnt0(ilai)/'glob.veg.lai.grib    '/
      data fnt0(ivg2)/'glob.veg.vegcov.grib '/
      data fnt0(irt1)/'glob.veg.rt1.grib    '/
      data fnt0(irt2)/'glob.veg.rt2.grib    '/
      data fnt0(irt3)/'glob.veg.rt3.grib    '/
#else
      data fnt0(ivgv)/'                     '/
      data fnt0(ivtv)/'                     '/
      data fnt0(ivrt)/'                     '/
      data fnt0(ibif)/'                     '/
      data fnt0(ids )/'                     '/
      data fnt0(idsm)/'                     '/
      data fnt0(iws )/'                     '/
      data fnt0(icef)/'                     '/
      data fnt0(iexp)/'                     '/
      data fnt0(ikst)/'                     '/
      data fnt0(idph)/'                     '/
      data fnt0(ibub)/'                     '/
      data fnt0(iqrt)/'                     '/
      data fnt0(ibkd)/'                     '/
      data fnt0(isld)/'                     '/
      data fnt0(iwcr)/'                     '/
      data fnt0(iwpw)/'                     '/
      data fnt0(ismr)/'                     '/
      data fnt0(islz)/'                     '/
      data fnt0(isnz)/'                     '/
      data fnt0(ilai)/'                     '/
      data fnt0(ivg2)/'                     '/
      data fnt0(irt1)/'                     '/
      data fnt0(irt2)/'                     '/
      data fnt0(irt3)/'                     '/
#endif
      data fnt0(isv1)/'                     '/
      data fnt0(itv1)/'                     '/
      data fnt0(isv2)/'                     '/
      data fnt0(itv2)/'                     '/
#ifdef ADJSOILW
      data fnt0(iobp)/'                     '/
#endif
c
c  forecast sfc fields.  blank indicate time invariant or not available.
c  for time invariant climatology, initial sfc fields needs to be persisted,
c  thus blank to avoid unnecessary io.
c
      character*128 fnfc(numgrbs)
      data fnfc(igla)/'                     '/ ! time invariant
      data fnfc(imxi)/'                     '/ ! time invariant
      data fnfc(imsk)/'                     '/ ! fnmask used 
      data fnfc(iais)/'clim.ice.grib        '/ ! monthly varying
      data fnfc(iscv)/'                     '/ ! predicted
      data fnfc(isno)/'                     '/ ! predicted
      data fnfc(itsf)/'clim.sst.grib        '/ ! monthly varying
      data fnfc(iab4)/'clim.yhalbedo.grib   '/ ! monthly varying
      data fnfc(iab1)/'clim.sibalbedo.grib  '/ ! monthly varying
      data fnfc(ialf)/'clim.yhalbedo.grib   '/ ! monthly varying
      data fnfc(iso2)/'                     '/ ! predicted
      data fnfc(isn1)/'                     '/ ! predicted
      data fnfc(izor)/'clim.sibrough.grib   '/ ! monthly varying
      data fnfc(iplr)/'clim.sibresis.grib   '/ ! monthly varying
      data fnfc(itg3)/'                     '/ ! time invariant
      data fnfc(ito2)/'                     '/ ! predicted
      data fnfc(itn1)/'                     '/ ! predicted
      data fnfc(iveg)/'clim.vegfrac40m.grib '/ ! monthly varying
      data fnfc(ivet)/'                     '/ ! time invariant
      data fnfc(isot)/'                     '/ ! time invariant
      data fnfc(islo)/'                     '/ ! time invariant
      data fnfc(isna)/'                     '/ ! time invariant
      data fnfc(ismx)/'                     '/ ! time invariant
      data fnfc(ismn)/'                     '/ ! time invariant
      data fnfc(ivgv)/'                     '/ ! time invariant
      data fnfc(ivtv)/'                     '/ ! time invariant
      data fnfc(ivrt)/'                     '/ ! time invariant
      data fnfc(ibif)/'                     '/ ! time invariant
      data fnfc(ids )/'                     '/ ! time invariant
      data fnfc(idsm)/'                     '/ ! time invariant
      data fnfc(iws )/'                     '/ ! time invariant
      data fnfc(icef)/'                     '/ ! time invariant
      data fnfc(iexp)/'                     '/ ! time invariant
      data fnfc(ikst)/'                     '/ ! time invariant
      data fnfc(idph)/'                     '/ ! time invariant
      data fnfc(ibub)/'                     '/ ! time invariant
      data fnfc(iqrt)/'                     '/ ! time invariant
      data fnfc(ibkd)/'                     '/ ! time invariant
      data fnfc(isld)/'                     '/ ! time invariant
      data fnfc(iwcr)/'                     '/ ! time invariant
      data fnfc(iwpw)/'                     '/ ! time invariant
      data fnfc(ismr)/'                     '/ ! time invariant
      data fnfc(islz)/'                     '/ ! time invariant
      data fnfc(isnz)/'                     '/ ! time invariant
#ifdef VIC
      data fnfc(ilai)/'glob.veg.lai.grib    '/ ! monthly varying
#else
      data fnfc(ilai)/'                     '/ ! monthly varying
#endif
      data fnfc(ivg2)/'                     '/ ! time invariant
      data fnfc(irt1)/'                     '/ ! time invariant
      data fnfc(irt2)/'                     '/ ! time invariant
      data fnfc(irt3)/'                     '/ ! time invariant
      data fnfc(isv1)/'                     '/ ! predicted 
      data fnfc(itv1)/'                     '/ ! predicted
      data fnfc(isv2)/'                     '/ ! predicted
      data fnfc(itv2)/'                     '/ ! predicted
#ifdef ADJSOILW
      data fnfc(iobp)/'                     '/ ! daily varying
#endif
c
      integer lugb
      integer iy,im,id,ih
      real fh
c
      integer is2g(numsfcs)
      logical lsf(numsfcs)
c
      character*128 fngrb(numgrbs)
      character*128 condir,bindir
c
      integer i,j,ij,k,n
#ifdef ADJSOILW
      integer intobsp
      real mdlpcp(LONF2S*LATG2S),mdlrnf(LONF2S*LATG2S)
#endif
c
      data condir/'    '/,bindir/'    '/
      namelist/namsfc/fnmskg,fnorog,fnmask,fnobsst0,
     &                fngrb,condir,bindir
c
      save fnmskg,fnorog,fnmask,fnobsst0,fngrb,condir,bindir
c
      integer ifp,jfp
      data ifp/0/,jfp/0/
      save ifp,jfp
c
#ifdef MP
      if( mype.eq.master ) then
#endif
c
      if(ifp.eq.0) then
        if(fh.eq.0.) then
           do n=1,numgrbs
             fngrb(n)=fnt0(n)
           enddo
        else
           do n=1,numgrbs
             fngrb(n)=fnfc(n)
           enddo
        endif
#ifdef RSM
        if(ifcst.eq.0) then
          nf=5
          close(nf)
          open(unit=nf,file='rsfc0.parm',form='formatted',status='old')
        else
          nf=2
          close(nf)
          open(unit=nf,file='rfcst.parm',form='formatted',status='old')
        endif
#else
        if(ifcst.eq.0) then
          nf=5
          close(nf)
          open(unit=nf,file='sfc0.parm',form='formatted',status='old')
        else
          nf=1
          close(nf)
          open(unit=nf,file='fcst.parm',form='formatted',status='old')
        endif
#endif
        read(nf,namsfc)
        ifp=1
      endif
#ifdef ADJSOILW
      if(fh.gt.0..and.fngrb(iobp)(1:4).eq.'    ') then
         print *,'adjsoilw run but obs precip file is empty'
         call abort
      endif
#endif
c
c  check input sfc file names for consistency
c  if condir and bindir are given, process them
c
c  note that if filenames starts from '/', it is assumed that
c  full directory is already given and condir, bindir are not added.
c
        call checkfn(condir,bindir,fngrb,numgrbs,
     1               fnmskg,fnorog,fnmask)
c
c  translate unpacked grib climatology/analysis file field record 
c  index and sfc file record index.
c     is2g(index of sfc file)=index of grb file
c
      call inds2g(is2g)
#ifdef DBG
#ifdef MP
      if(mype.eq.master) then
#endif
      print *,'is2g'
      do n=1,numsfcs
        print *,n,is2g(n)
      enddo
      print *,'ksfc'
      do n=1,numsfcs
        print *,n,ksfc(n)
      enddo
#ifdef MP
      endif
#endif
#endif
c
c  initialize file exist flag
c
      do n=1,numsfcs
        lsf(n)=.false.
      enddo
c
c  initialize grib read in array
c
      do k=1,numsfcs
        do ij=1,idim*jdim
          grbfld(ij,k)=0.
        enddo
      enddo
c
c
c
c  read sfc fields from grib (climatology as well as analysis)
c  in the order of forecast sfc file record
c
      call rdgrbsfc(lugb,idim,jdim,
     &              numsfcs,numsfcv,ksfc,
     &              is2g,
     &              iy,im,id,ih,fh,
     &              fnmask,fnorog,fnmskg,fngrb,numgrbs,
     &              orog,slmask,
     &              grbfld,lsf)
c
c  process multiple grib files that correspond to a single
c  forecast sfc record.  currently, snow/snow-cover belongs
c  to this category.
c
c  read snow cover field (snow depth is read by default)
c  if snow depth is missing and snow cover is given,
c  do special processing in sfcmrg.  set snow depth to -999.,
c  where the snow cover is observed.
c
c
c
c
      if(.not.lsf(jsno).and.fngrb(iscv)(1:4).ne.'    ') then
#ifdef DBG
#ifdef MP
        if(mype.eq.master) then
#endif
        print *,'iscv superead, iscv=',iscv
#ifdef MP
        endif
#endif
#endif
        call superead(lugb,fngrb(iscv),idim,jdim,slmask,
     1                iscv,iy,im,id,ih,fh,
     2                fnmskg,grbfld(1,jsno))
        lsf(jsno)=.true.
c
c  in order to distinguish from actual snow depth, set snow mask
c  as -999.
c

        do ij=1,idim*jdim
          grbfld(ij,jsno)=grbfld(ij,jsno)*(-999.)
        enddo
      endif
#ifdef GCA
c
      call gulfofca(grbfld(1,jtsf),slmask)
#endif
c
c
c
c
c  obtain glacier points and maximum possible ice extent from
c  grib file
c
c  get glacier points from grib input file
c
      if(jfp.eq.0) then
#ifdef DBG
        print *,'igla superead, igla=',igla
#endif
        if(fngrb(igla)(1:4).eq.'    ') then
           print *,'glacier grib file name empty'
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
        call superead(lugb,fngrb(igla),idim,jdim,slmask,
     1                igla,iy,im,id,ih,fh,
     2                fnmskg,glacier)
c
c  get maximum possible sea ice extent from grib file
c
#ifdef DBG
#ifdef MP
        if(mype.eq.master) then
#endif
        print *,'imxi superead, imxi=',imxi
#ifdef MP
        endif
#endif
#endif
        if(fngrb(imxi)(1:4).eq.'    ') then
           print *,'max ice grib file name empty'
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
        call superead(lugb,fngrb(imxi),idim,jdim,slmask,
     1                imxi,iy,im,id,ih,fh,
     2                fnmskg,amxice)
        jfp=1
      endif
c
c  sst anomaly persistency
c  Note that superead is called with ft=0. to compute initial anomaly.
c  If you wish, make factor as function of fhour.
c  For this to work, fnobsst0 should be non blank and fnfc must be
c  a climatological sst.
c  
      if(fnobsst0(1:4).ne.'    '.and.fngrb(itsf)(1:4).ne.'    ') then
        print *,'this is sst anomaly persistency run.'
        allocate(obsst0(ijdim),clsst0(ijdim))
        call superead(lugb,fnobsst0,idim,jdim,slmask,
     1                itsf,iy,im,id,ih,0.,
     2                fnmskg,obsst0)
        call superead(lugb,fngrb(itsf),idim,jdim,slmask,
     1                itsf,iy,im,id,ih,0.,
     2                fnmskg,clsst0)
        factor=1.
        do ij=1,ijdim
          grbfld(ij,jtsf)=grbfld(ij,jtsf)+
     1                    (obsst0(ij)-clsst0(ij))*factor
        enddo
        deallocate(obsst0,clsst0)
      endif
c
c  rearrange all the arrays to transform format before merge
c
#ifndef RSM
        call row1ns(orog)
        call row1ns(slmask)
        call row1ns(glacier)
        call row1ns(amxice)
        do k=1,numsfcs
          call row1ns(grbfld(1,k))
        enddo
#else
        call shalfi(orog,1)
        call shalfi(slmask,1)
        call shalfi(glacier,1)
        call shalfi(amxice,1)
        call shalfi(grbfld,numsfcs)
#endif
c
#ifdef MP
      endif
#ifndef RSM
#define MPBCASTL mpbcastl
#define MPBCASTI mpbcasti
#define MPBCASTR mpbcastr
#define MPGF2P   mpgf2p
#define MPSYNALL mpsynall
#else
#define MPBCASTL rmpbcastl
#define MPBCASTI rmpbcasti
#define MPBCASTR rmpbcastr
#define MPGF2P   rmpgf2p
#define MPSYNALL rmpsynall
#endif
c
      call MPBCASTL(lsf    ,numsfcs)
      call MPBCASTI(is2g   ,numsfcs)
      call MPBCASTR(orog   ,  ijdim)
      call MPBCASTR(slmask ,  ijdim)
      call MPBCASTR(glacier,  ijdim)
      call MPBCASTR(amxice ,  ijdim)
      call MPBCASTR(grbfld ,ijdim*numsfcs)
c
      call MPSYNALL
      call MPGF2P(orog  ,LONF2F,LATG2F,  orogp,LONF2S,LATG2S,1)
      call MPGF2P(slmask,LONF2F,LATG2F,slmaskp,LONF2S,LATG2S,1)
      call MPGF2P(glacier,LONF2F,LATG2F,glacierp,
     1            LONF2S,LATG2S,1)
      call MPGF2P(amxice,LONF2F,LATG2F,amxicep,LONF2S,LATG2S,1)
      call MPGF2P(grbfld,LONF2F,LATG2F,grbfldp,
     1            LONF2S,LATG2S,numsfcs)
#endif
c
c  merge with forecast sfc array
c
#ifdef MP
#define GRBFLD grbfldp
#define OROG   orogp
#define SLMASK slmaskp
#define GLACIER glacierp
#define AMAXICE amxicep
#else
#define GRBFLD grbfld
#define OROG   orog
#define SLMASK slmask
#define GLACIER glacier
#define AMAXICE amxice
#endif
      call sfcmrg(GRBFLD,LONF2S,LATG2S,iy,im,id,ih,fh,
     &            OROG,SLMASK,GLACIER,AMAXICE,
#ifdef ADJSOILW
     &            intobsp,mdlpcp,mdlrnf,
#endif
     &            lsf,is2g)
c
c  rearrange back to input format 
c
#ifndef RSM
#ifdef MP
        if(mype.eq.master) then
#endif
        call rowsep(orog)
        call rowsep(slmask)
        call rowsep(glacier)
        call rowsep(amxice)
#ifdef MP
        endif
#endif
#else
#ifdef MP
        if(mype.eq.master) then
#endif
        call shalfo(orog,1)
        call shalfo(slmask,1)
        call shalfo(glacier,1)
        call shalfo(amxice,1)
#ifdef MP
        endif
#endif
#endif
c
#ifdef NFDR
c
c  nfdr
c
#ifdef RSM
      allocate (flat(ijdim),flon(ijdim))
      call setgrd(flat,flon,rdelx,rdely,dlamda0)
      rlat1=flat(1)/rad
      rlon1=flon(1)/rad
      rlat2=flat(idim*2)/rad
      rlon2=flon(idim*2)/rad
      deallocate (flat,flon)
#endif
#ifdef DBGNFDR
#ifdef MP
      if(mype.eq.master) then
#endif
      fh=24
      open(unit=33,file='./nfdr.forcing.ft24',form='unformatted',
     1     status='old')
      read(33) st2m2pm
      read(33) srh2m2pm
      read(33) scloud2pm
      read(33) swind2pm
      read(33) sppt24
      read(33) sprcpdur
      read(33) stmax24
      read(33) stmin24
      read(33) srhmax
      read(33) srhmin
      do ij=1,ijdim
        t2m2pm(ij,1)=st2m2pm(ij)        ! degree kelvin
        rh2m2pm(ij,1)=srh2m2pm(ij)      ! percent
        cloud2pm(ij,1)=scloud2pm(ij)    ! fraction
        wind2pm(ij,1)=abs(swind2pm(ij)) ! m/s
        ppt24(ij,1)=sppt24(ij)          ! meter unit
        prcpdur(ij,1)=sprcpdur(ij)      ! in second
        tmax24(ij,1)=stmax24(ij)        ! kelvin
        tmin24(ij,1)=stmin24(ij)        ! kelvin
        rhmax(ij,1)=srhmax(ij)          ! percent
        rhmin(ij,1)=srhmin(ij)          ! percent
      enddo
#ifdef RSM
      call shalfi(t2m2pm,1)
      call shalfi(rh2m2pm,1)
      call shalfi(cloud2pm,1)
      call shalfi(wind2pm,1)
      call shalfi(ppt24,1)
      call shalfi(prcpdur,1)
      call shalfi(tmax24,1)
      call shalfi(tmin24,1)
      call shalfi(rhmax,1)
      call shalfi(rhmin,1)
#else
      call row1ns(t2m2pm)
      call row1ns(rh2m2pm)
      call row1ns(cloud2pm)
      call row1ns(wind2pm)
      call row1ns(ppt24)
      call row1ns(prcpdur)
      call row1ns(tmax24)
      call row1ns(tmin24)
      call row1ns(rhmax)
      call row1ns(rhmin)
#ifdef MP
      endif
#endif
#endif
#endif
c
c  substitute grbfld(ij,jobp) into ppt24 with proper unit conversion should be
c  able to use observed precipitation for nfdr.
c
      if(fh.ne.0.) then
#ifdef ADJSOILW
        ij=0
        do j=1,LATG2S
          do i=1,LONF2S
            ij=ij+1
            ppt24(ij,1)=GRBFLD(ij,jobp)
          enddo
        enddo
c
#endif
        call nfdr_pred(iy,im,id,ih,fh,condir)
c
      endif
#endif
c
      return
      end
