#define G g
#define IGRD1 igrd1
#define IGRD igrd
#define IWAV1 iwav1
#define JCAP jcap
#define JGRD1 jgrd1
#define JGRD jgrd
#define JWAV1 jwav1
#define NCLDG ncldg
#define NTOTAL ntotal
#define NTRAC ntrac
#define LATG latg
#define LEVP1 levp1
#define LEVS levs
#define LONF lonf
#define PI pi
#define RERTH
#define RD rd
#define RV rv
      program p2sig
      implicit none
#include <paramodel.h>
#define G2G
#undef DBG
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    p2sig      vertically interpolate standard pressure level
c                           data to sigma level and convert to either
c                           spectral coefficients or RSM grid suitable
c                           for G-RSM input.
c
c                           This program assumes that all the input standard
c                           pressure level fields are given in one grib format
c                           file (file name fn given as a namelist input).
c
c   prgmmr: masao kanamitsu date: 07-03-28
c           succeeded by kei yoshimura
c           haiqin li added regional to regional: 09-06-04  
c
c abstract: this routine spatially interpolates standard pressure level data
c           to forecast model gaussian grid, or rsm grid, interpolates surface
c           pressure to the new orography, vertically interpolates wind,
c           temperature and moisture, fixes moisture in the stratosphere if
c           necessary, then transform the data either to spectral coeffs
c           or rsm grid suitable for G-RSM input.
c           if this program is used to create lateral forcing for rsm
c           on rsm grid, you need to run this program using larger domain
c           rsm as a model configuration.
c           If defined G2G, p2sig will implement Global interpolation, else the
c           regional pressure to regional sigma interpolation will be implemented.
c
c           The defination of G2G in: p2sig.F setgrd.F shalfi.F shalfo.F prsgrdbrd.F
c  
c attributes:
c   language: fortran 90.
c
      integer idim,jdim,idimt,jdimhf,ijdim
#ifdef G2G
      parameter(idim=lonf_,jdim=latg_,idimt=idim*2,jdimhf=jdim/2,
     1          ijdim=idim*jdim)
#else
      parameter(idim=igrd1_,jdim=jgrd1_,idimt=idim,jdimhf=jdim,
     1          ijdim=idim*jdim)
      integer levmax
      parameter(levmax=100)
      real  dummy(2*levmax),ext(512-(6+2*levmax))
#endif
      integer kdim,kdimp,kdimh
      parameter(kdim=levs_,kdimp=kdim+1,kdimh=levs_*ntotal_)
c
      integer nfld
      parameter(nfld=5+ntotal_-1)
      integer mkd
      parameter(mkd =(kdim*nfld+2)*2)
c
      real xlat(idimt)
      real ps(idimt),zs(idimt),expps(idimt)
      real fsig(idimt,kdim*nfld+2)
!     for regional to reginal output 
      real rfsig(idimt,jdimhf,kdim*nfld+2)
      real pi,twopi,hfpi,qtpi,rad
      parameter(pi=pi_,twopi=2.0*pi,hfpi=0.5*pi,qtpi=0.5*hfpi)
      parameter(rad=pi/180.)
 
!
      real qe(idimt,kdim)
c
      character*32 lab
      character*8 sfcftyp
#ifdef G2G
      real waves,xlayers,trun,order,realform,gencode,
     1     rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers,
     2     cloud,pdryini,subcen,ppid,slid,vcid,vmid,vtid,
     3     runid,usrid
      integer kdum,kdum2,kens
      parameter(kdum=201-levs_-1-levs_,kdum2=21,kens=2)
      real dummy(kdum),dummy2(kdum2),ensemble(kens)
#else
#include <rscomloc.h>
      real  truth,orient,a2,x00,y1,y0,y00,cenlon,cenlat,rsoa2
      real  fm2(idim*jdim),fm2x(idim*jdim),fm2y(idim*jdim)
      real  flat(idim*jdim),flon(idim*jdim)
      real  rgln,rcsln(idim*jdim),rsnln(idim*jdim)
      real  dlamda0
      real  gz2(idim,jdim), ps2(idim,jdim)
      real  te2(idim,jdim,kdim),uu2(idim,jdim,kdim)
      real  vv2(idim,jdim,kdim),rq2(idim,jdim,kdim)
      integer idate2(4)
      real  fhour2
#endif
      integer idate(4)
      real fhour
c
#ifdef G2G
c  arrays for global spectral conversion
c
      integer mdim,mdimv
      parameter(mdim =(jcap_+1)*(jcap_+2))
      parameter(mdimv=(jcap_+1)*(jcap_+4))
      real gzo,qo,teo,dio,zeo,rqo
      common gzo(mdim),qo (mdim),teo(mdim,kdim),dio(mdim,kdim),
     1       zeo(mdim,kdim),rqo(mdim,kdimh)
c
      real zss,pss,tts,qqs,uus,vvs
      common zss(idimt),pss(idimt),tts(idimt,kdim),
     1       uus(idimt,kdim),vvs(idimt,kdim),qqs(idimt,kdimh)
      real zsa,psa,tta,qqa,uua,vva
      common zsa(idimt),psa(idimt),tta(idimt,kdim),
     1       uua(idimt,kdim),vva(idimt,kdim),qqa(idimt,kdimh)
c
      real eps,colrad,wgt,wgtcs,rcs2
      common/plncom/eps   (mdimv ),colrad(jdimhf),wgt(jdimhf ),
     1              wgtcs (jdimhf),rcs2(jdimhf)
c
      real qlnt(mdim),qlnv(mdimv),qdert (mdim),qlnwct(mdim)
c
      real b(mdim)
      integer indxnn,indxmm
      common/comind/ indxnn(mdim),indxmm(mdim)
#endif
c
c  model orography
c
      real zso(idimt,jdimhf)
c
c  model sigma level
c
      real si(kdimp),sl(kdim),del(kdim)
c
#ifdef DBG
      real zout(idimt,jdimhf),pout(idimt,jdimhf)
      real tout(idimt,jdimhf,kdim),qout(idimt,jdimhf,kdimh)
      real uout(idimt,jdimhf,kdim),vout(idimt,jdimhf,kdim)
      real*4 outss(idim,jdim)
      real*4 outs(idim,jdim,kdim,nfld),orog(idim,jdim)
#endif
      integer ij,tt,nn
c
      real tensn
      parameter(tensn=10.0)
      real gamma,ps1,rh1,psk
      parameter(gamma=6.5e-3)
      integer lugb
      data lugb/10/
c
      integer iy,im,id,ih,indmois
      real    fh
      integer maxstdprs
      parameter(maxstdprs=100)
      real    stdprs(maxstdprs)
      character*128 fn
      data indmois/0/
      integer isfcp
      data isfcp/1/
      integer inst
      data inst/0/ !! 0: ncep rean grib id, 1: ecmwf rean grib id
      integer iintp
      data iintp/0/
      integer irhex
      data irhex/0/
      integer ihyb
      data ihyb/0/ !! 0: pressure level, 1: s-p hybrid level, 2: sigma level
      real hyal(0:26),hybl(0:26) !! top to bottom (ccsm coordinate)
      data hyal/1.e-10,0.003544638,0.0073888135,0.013967214,0.023944625,
     $     0.03723029,0.053114605,0.07005915,0.07791257,0.07660701,
     $     0.075071085,0.07326415,0.071138385,0.068637535,0.065695415,
     $     0.062234155,0.058162165,0.05337168,0.047735925,0.041105755, 
     $     0.0333057,0.02496844,0.01709591,0.01021471,0.004803175,
     $     0.00126068,0./
      data hybl/
     $     0.,0.,0.,0.,0.,0.,0.,0.,0.007526545,0.023907685,0.04317925,
     $     0.065851245,0.092523685,0.1239024,0.16081785,0.204247,
     $     0.2553391,0.3154463,0.3861593,0.4693495,0.5672185,0.67182785,
     $     0.77060615,0.85694605,0.9248457,0.96929415,0.9925561/
      integer inorhv
      data inorhv/0/
      integer izrtr
      data izrtr/0/
c
      integer kdimprs
c
      real es,fpvs0,qs,xeps,xepsm1
      parameter(xeps=rd_/rv_,xepsm1=rd_/rv_-1.)
c
      integer i,j,in,in2,k,n,m,lat,kdimi,kdimqi,ie,iw
c
      real, allocatable :: out(:,:,:,:),out2(:,:,:)
      real, allocatable :: slin(:,:)
      real, allocatable :: siin(:,:)
c
c   lugb : unit number for grib read file
c   fn   : input standard pressure level grid file name
c   iy,im,id,ih,fh   : input year,month,day,hour,forecast hour
c   stdprs : standard pressure level to be read in mb unit
c   kdimprs : number of standard pressure levels
c   indmois : =0 specific humidity ; =1 relative humidity
c   iorsmgsm : choice of output. =0 global spectral ; =1 regional
c   isfcp : =0 use p-level geoheight ; =1 use sfcP and topo
c   isfcp : =2 use sea level pressure and zero-topo
c   irhex : =0 use extrapolation for RH (default)
c   irhex : =1 do not extrapolate for upper than top level
c   inorhv: =0 use RH for vertical interpolation
c   inorhv: =1 use SH for vertical interpolation
c   izrtr:  =0 requires real input for tracers
c   izrtr:  =1 assume zero tracers (input data not required)
c
      namelist/namp2s/fn,iy,im,id,ih,fh,stdprs,indmois,isfcp,inst,iintp
     $     ,irhex,ihyb,inorhv,izrtr
c
      do k=1,maxstdprs
        stdprs(k)=0.
      enddo
c
      read (5,namp2s)
      write(6,namp2s)
c
      if (ihyb.eq.1 .and. isfcp.ne.1) then
        print*,'if ihyb=1, isfcp should be 1'
        call abort
      endif
      if (ihyb.eq.2 .and. isfcp.ne.1) then
        print*,'if ihyb=2, isfcp should be 1'
        call abort
      endif
      kdimprs=1
      dowhile(stdprs(kdimprs).ne.0.)
        kdimprs=kdimprs+1
      enddo
      kdimprs=kdimprs-1
      print *,'kdimprs=',kdimprs
c
c  define constants for global spherical transform
c
#ifdef G2G
      call glats (jdimhf, colrad, wgt, wgtcs, rcs2)
      call epslon(eps,jcap_)
      call cmpind
#endif
c
c  define model sigma levels
c
      call newsig2(si,sl,del)
c
c define constants for saturation vapor pressure
c
      call gpvs
c
c  read model surface topography (grid space).
c
      read(11) zso
c
      in=0
      in2=0
c
c interpolate/extrapolate from slin to sl
c
      print *,' this program converts from pres to sl '
      print *,' pres=',(stdprs(k),k=1,kdimprs)
c
c  zero out spectral array for accumulation in lat loop.
c
#ifdef G2G
      do i = 1, mdim
        gzo(i)=0.0
        qo(i)=0.0
      enddo
      do k=1,kdim
        do i = 1, mdim
           teo(i,k)=0.0
           dio(i,k)=0.0
           zeo(i,k)=0.0
        enddo
      enddo
      do k=1,kdimh
        do i = 1, mdim
           rqo(i,k)=0.0
        enddo
      enddo
#endif
c
      allocate (out(idimt,jdimhf,kdimprs,nfld))
      allocate (out2(idimt,kdimprs,nfld))
      allocate (slin(idimt,kdimprs))
      allocate (siin(idimt,kdimprs+1))
c
c  fsig(1,1)=zg
c  fsig(1,2)=ps
c  fsig(1,3)=t
c  fsig(1,3+kdim*1)=u
c  fsig(1,3+kdim*2)=v
c  fsig(1,3+kdim*3)=q
c  fsig(1,3+kdim*4)=o3
c  fsig(1,3+kdim*5)=cl
c  fsig(1,3+kdim*6)=pr
c
      call prsgrbrd(lugb,fn,iy,im,id,ih,fh,stdprs,kdimprs,indmois,nfld,
     1              out,idim,jdim,isfcp,inst,iintp,ihyb,izrtr)
#ifndef G2G
      print*,'rproj,rorient,rtruth,rdelx,rdely=',
     1        rproj,rorient,rtruth,rdelx,rdely
      print*,'rlat1,rlon1,rlat2,rlon2=',rlat1,rlon1,rlat2,rlon2
#endif
c
#ifdef G2G
#ifdef DBG
      open(52,recl=4*idim*jdim,access='direct')
      do j=1,jdimhf
        do i=1,idim
          orog(i,j*2-1)=zso(i,j)
        enddo
        do i=idim+1,idimt
          orog(i-idim,j*2)=zso(i,j)
        enddo
      enddo
      write(52,rec=1) orog
      do n=1,nfld
        do k=1,kdimprs
          do j=1,jdimhf
            do i=1,idim
              outs(i,j*2-1,k,n)=out(i,j,k,n)
            enddo
            do i=idim+1,idimt
              outs(i-idim,j*2,k,n)=out(i,j,k,n)
            enddo
          enddo
        enddo
      enddo
      tt=1
      do n=1,nfld
        do k=1,kdimprs
          tt=tt+1
          write(52,rec=tt) ((outs(i,j,k,n),i=1,idim),j=1,jdim)
        enddo
      enddo
#endif
      do n=1,nfld
        do k=1,kdimprs
           call row1ns2(out(1,1,k,n))
        enddo
      enddo
      call row1ns2(zso)
#endif
c
c  at this point grid, values are available
c      out(1,1,1,1) ... wind x-direction
c      out(1,1,1,2) ... wind y-direction
c      out(1,1,1,3) ... temperature
c      out(1,1,1,4) ... moisture
c      if isfcp=0
c       out(1,1,1,5) ... geopotential height 
c      if isfcp=1
c       out(1,1,1,5) ... surface pressure (Pa)
c       out(1,1,2,5) ... surface height 
c      if isfcp=2
c       out(1,1,1,5) ... sealevel pressure (Pa)
c       out(1,1,2,5) ... zero topography
c      out(1,1,1,6) ... o3
c      out(1,1,1,7) ... cloud
c      out(1,1,1,8) ... precipitation
c
c  note that idimt=idim and jdimhf=jdim for rsm
c
!  begine the long loop
      do lat=1,jdimhf
         print*,'kei',lat,isfcp
c
        if (isfcp.eq.0) then
c  ps=1000mb=100cb (stdprs(1))
c
        do i=1,idimt
          ps(i)=log(stdprs(1)/10.)
c          ps(i)=log(1000.)
          zs(i)=out(i,lat,1,5)
          fsig(i,1)=zso(i,lat)
        enddo
c
        else
c
        do i=1,idimt
          if (out(i,lat,1,5).eq.0) then
            ie=i-1
            iw=i+1
            if (i.eq.1) then
              ie=idimt
            elseif (i.eq.idimt) then
              iw=1
            endif
            if ((out(ie,lat,1,5).ne.0.).and.
     $           (out(iw,lat,1,5).ne.0.)) then
              ps(i)=log((out(ie,lat,1,5)+out(iw,lat,1,5))*0.5/1000.)
              zs(i)=(out(ie,lat,2,5)+out(iw,lat,2,5))*0.5
              print*,'zero ps. averaging.',i,lat,exp(ps(i))*10.
            else
              print*,'invalid ps',i,lat
              call abort
            endif
          else
            ps(i)=log(out(i,lat,1,5)/1000.)
            zs(i)=out(i,lat,2,5)
          endif
          fsig(i,1)=zso(i,lat)
        enddo
c
        endif   
c
cc sigma level coordination cc
c
      if (isfcp.eq.0) then !! surface = lowest p-level
        do i=1,idimt
        do k=1,kdimprs
          slin(i,k)=stdprs(k)/stdprs(1)
        enddo
        siin(i,1)=1.D0
        siin(i,kdimprs+1)=0.D0
        do k=2,kdimprs
          siin(i,k)=(slin(i,k-1)+slin(i,k))/2.D0
        enddo
        enddo
      elseif (isfcp.ge.1) then !! surface = new surface (out(i,j,1,5))
                               !! unit of input ps is Pa
        if (ihyb.eq.0) then
          do i=1,idimt
            do k=1,kdimprs
              slin(i,k)=stdprs(k)/(exp(ps(i))*10.)
            enddo
          enddo
        elseif (ihyb.eq.1) then
          do i=1,idimt
            do k=1,kdimprs
              slin(i,k)=hyal(int(stdprs(k)))*1000./(exp(ps(i))*10.)
     $            +hybl(int(stdprs(k)))
            enddo
          enddo
        elseif (ihyb.eq.2) then
          do i=1,idimt
            do k=1,kdimprs
              slin(i,k)=stdprs(k)/10000.
            enddo
          enddo
        endif
        do i=1,idimt
        siin(i,1)=1.D0
        siin(i,kdimprs+1)=0.D0
c        siin(i,1)=min(stdprs(1)/out(i,lat,1,5)*100.,1)
c        siin(i,1)=stdprs(1)/out(i,lat,1,5)*100.
        do k=2,kdimprs
c          siin(i,k)=min((stdprs(k-1)+stdprs(k))/2.D0/
c     $          out(i,lat,1,5)*100.,1)
c          siin(i,k)=(stdprs(k-1)+stdprs(k))/2.D0/
c     $          (exp(ps(i))*10.)
          if (slin(i,k-1).lt.slin(i,k)) then
             print*,'strange vertical coord.',i,k,slin(i,k-1),slin(i,k)
          endif
          siin(i,k)=(slin(i,k-1)+slin(i,k))/2.D0
        enddo
        enddo
      endif
c
c      print *,' slin=',(slin(1,k),k=1,kdimprs)
c      print *,' siin=',(siin(1,k),k=1,kdimprs)
c      print *,' sl=  ',sl      
c
#ifdef G2G
        call pln2t(qlnt,qlnv,colrad,lat)
c
c  compute surface pressure from model topography
c
c        print'(i5,4f15.7)',
c     $       (i,exp(ps(i)),out(i,lat,1,5),zs(i),fsig(i,1),
c     $       i=1,idimt)
c
c  u,v multiplied by coslat
c
        do i=1,idimt
          do k=1,kdimprs
            out(i,lat,k,1)=out(i,lat,k,1)*cos(pi_/2.0-colrad(lat))
            out(i,lat,k,2)=out(i,lat,k,2)*cos(pi_/2.0-colrad(lat))
          enddo
        enddo
#endif
c
c  convert rh to q for newps
c
        if(indmois.eq.1) then
c
          if (ihyb.eq.0) then
c
          do k=1,kdimprs
            do i=1,idimt
c#ifdef ICE
c             es=fpvs(out2(i,k,3))
c#else
              es=fpvs0(out(i,lat,k,3)) !! cb
c#endif
              qs=xeps*es/(stdprs(k)/10.+xepsm1*es)
              out(i,lat,k,4)=out(i,lat,k,4)/100.*qs
            enddo
          enddo
c
          elseif (ihyb.eq.1) then
c
          do k=1,kdimprs
            do i=1,idimt
c#ifdef ICE
c             es=fpvs(out2(i,k,3))
c#else
              es=fpvs0(out(i,lat,k,3)) !! cb
c#endif
              ps1=hyal(int(stdprs(k)))*1000.+
     $             hybl(int(stdprs(k)))*(exp(ps(i))*10.)
              qs=xeps*es/(ps1/10.+xepsm1*es)
              out(i,lat,k,4)=out(i,lat,k,4)/100.*qs
            enddo
          enddo
c
          elseif (ihyb.eq.2) then
c
          do k=1,kdimprs
            do i=1,idimt
c#ifdef ICE
c             es=fpvs(out2(i,k,3))
c#else
              es=fpvs0(out(i,lat,k,3)) !! cb
c#endif
              ps1=stdprs(k)/10000.*(exp(ps(i))*10.)
              qs=xeps*es/(ps1/10.+xepsm1*es)
              out(i,lat,k,4)=out(i,lat,k,4)/100.*qs
            enddo
          enddo
c
          endif
c
        endif
c
c  Temporarily, simple extrapolation for undefined (T=0) values for T
c
        if (isfcp.ge.1) then
        do i=1,idimt
          do k=kdimprs,1,-1
            if (ihyb.eq.0) then
              ps1=log(stdprs(k+1)/10.)
              psk=log(stdprs(k)/10.)
            elseif (ihyb.eq.1) then
              ps1=log(hyal(int(stdprs(k+1)))*100.+
     $            hybl(int(stdprs(k+1)))*(exp(ps(i))))
              psk=log(hyal(int(stdprs(k)))*100.+
     $            hybl(int(stdprs(k)))*(exp(ps(i))))
            elseif (ihyb.eq.2) then
              ps1=log(stdprs(k+1)/10000.*exp(ps(i)))
              psk=log(stdprs(k)/10000.*exp(ps(i)))
            endif
            if (out(i,lat,k,3).eq.0) then
              out2(i,k,1)=out2(i,k+1,1)
              out2(i,k,2)=out2(i,k+1,2)
c              if (isfcp.eq.1) then
              out2(i,k,3)=out2(i,k+1,3)*exp(gamma*rd_/g_*(psk-ps1))
c              elseif (isfcp.eq.2) then
c              out2(i,k,3)=out2(i,k+1,3)
c                 print*,'Kei Warning',i,k,out2(i,k,3)
c              endif
              out2(i,k,4)=out2(i,k+1,4)
            else
              out2(i,k,1)=out(i,lat,k,1)
              out2(i,k,2)=out(i,lat,k,2)
              out2(i,k,3)=out(i,lat,k,3)
              out2(i,k,4)=out(i,lat,k,4)
            endif
            if (ntotal_.gt.1) then
              do nn=2,ntotal_
                out2(i,k,4+nn)=out(i,lat,k,4+nn)
              enddo
            endif
          enddo
          out2(i,1,5)=out(i,lat,1,5)
          out2(i,2,5)=out(i,lat,2,5)
        enddo
c
        elseif (isfcp.eq.0) then
        do n=1,nfld
          do k=1,kdimprs
            do i=1,idimt
              out2(i,k,n)=out(i,lat,k,n)
            enddo
          enddo
        enddo
        endif
c
c  compute surface pressure for new topography
c  newps2: uses input ps or Z1000 only to estimate sfc pressure. (isfcp=1,2)
c  newps3: uses all Z to estimate sfc pressure. (better if isfcp=0)
c
        if (isfcp.ge.1) then
        call newps2(out2(1,1,3),out2(1,1,4),
     $       kdimprs,ps,zs,fsig(1,1),fsig(1,2),siin)
c                                !! unit of ps and fsig(i,2) : ln(cb)
c                                !! unit of slin : [ND]
        elseif (isfcp.eq.0) then
        call newps3(out2(1,1,3),out2(1,1,4),
     $       kdimprs,ps,out2(1,1,5),fsig(1,1),fsig(1,2),siin,slin)
                                !! unit of ps and fsig(i,2) : ln(cb)
                                !! unit of slin : [ND]
        endif
c        print*,'kei'
c        print'(i5,5f15.7)',
c     $       (i,exp(ps(i))*10.,exp(fsig(i,2))*10.,
c     $       zs(i),fsig(i,1),siin(i,1),
c     $       i=1,idimt)
c
        if (inorhv.eq.0) then
c
c  convert q to rh before vertical interpolation
c        
        if (ihyb.eq.0) then
c
        do k=1,kdimprs
          do i=1,idimt
c#ifdef ICE
c            es=fpvs(out2(i,k,3))
c#else
            es=fpvs0(out2(i,k,3)) !! cb
c#endif
            qs=xeps*es/(stdprs(k)/10.+xepsm1*es)
            out2(i,k,4)=min(max(out2(i,k,4)/qs,0.),1.)*100.
          enddo
        enddo
c
        elseif (ihyb.eq.1) then
c
        do k=1,kdimprs
          do i=1,idimt
c#ifdef ICE
c            es=fpvs(out2(i,k,3))
c#else
            es=fpvs0(out2(i,k,3)) !! cb
c#endif
            ps1=hyal(int(stdprs(k)))*1000.+
     $           hybl(int(stdprs(k)))*exp(ps(i))*10.
            qs=xeps*es/(ps1/10.+xepsm1*es)
            out2(i,k,4)=min(max(out2(i,k,4)/qs,0.),1.)*100.
          enddo
        enddo
c
        elseif (ihyb.eq.2) then
c
        do k=1,kdimprs
          do i=1,idimt
c#ifdef ICE
c            es=fpvs(out2(i,k,3))
c#else
            es=fpvs0(out2(i,k,3)) !! cb
c#endif
            ps1=stdprs(k)/10000.*exp(ps(i))*10.
            qs=xeps*es/(ps1/10.+xepsm1*es)
            out2(i,k,4)=min(max(out2(i,k,4)/qs,0.),1.)*100.
          enddo
        enddo
c
        endif
c
        endif
c
c  vertical interpolation in the order of u, v, t, q
c
        call sg2sg(ps,slin,out2(1,1,1),fsig(1,2),sl,fsig(1,1*kdim+3),
     1             idimt,1,kdimprs,kdim,in,tensn,1)
        in=0
        call sg2sg(ps,slin,out2(1,1,2),fsig(1,2),sl,fsig(1,2*kdim+3),
     1             idimt,1,kdimprs,kdim,in,tensn,1)
        in=0
        call sg2sg(ps,slin,out2(1,1,3),fsig(1,2),sl,fsig(1,0*kdim+3),
     1             idimt,1,kdimprs,kdim,in,tensn,1)
        in=0
        if (inorhv.eq.0) then
           call sg2sgrh(ps,slin,out2(1,1,4),fsig(1,2),
     $          sl,fsig(1,3*kdim+3),
     1          idimt,1,kdimprs,kdim,in2,tensn,1,irhex)
           in2=0
        elseif (inorhv.eq.1) then
           call sg2sg(ps,slin,out2(1,1,4),fsig(1,2),
     $          sl,fsig(1,3*kdim+3),
     1          idimt,1,kdimprs,kdim,in2,tensn,1)
           in2=0
        endif
        if (ntotal_.gt.1) then
          do nn=2,ntotal_
            call sg2sg(ps,slin,out2(1,1,4+nn),fsig(1,2)
     $            ,sl,fsig(1,(2+nn)*kdim+3),
     1            idimt,1,kdimprs,kdim,in,tensn,1)
            in=0
          enddo
        endif
c
c sg2sg sets values to constant outside of input domain.
c fix lapse rate and relative humidity below input surface.
c
      do 299 i=1,idimt
      ps1=log(slin(i,1))+ps(i)
      do 298 k=1,kdim
      psk=log(sl(k))+fsig(i,2)
      if(psk.lt.ps1) goto 299
      fsig(i,k+2)=out2(i,1,3)*exp(gamma*rd_/g_*(psk-ps1))
 298  continue
 299  continue
c
#ifdef G2G
#ifdef DBG
c
c        print *,'t :(out2(i,10,3),i=1,10)=',(out2(i,10,3),i=1,10)
c        print *,'rh:(out2(i,10,4),i=1,10)=',(out2(i,10,4),i=1,10)
c        print *,'u :(out2(i,10,1),i=1,10)=',(out2(i,10,1),i=1,10)
c        print *,'v :(out2(i,10,2),i=1,10)=',(out2(i,10,2),i=1,10)
c
        print *,'zs:(out2(1,2,5))=',out2(1,2,5)
        print *,'zs:(fsig(1,1))=',fsig(1,1)
        print *,'ps:(ps(1))=',exp(ps(1))*10.
        print *,'ps:(fsig(1,2))=',exp(fsig(1,2))*10.
        print *,'slin:(slin(1,k)=',
     $       (slin(1,k),k=1,kdimprs)
        print *,'sl:(sl(k)=',
     $       (sl(k),k=1,kdim)
        print *,'t :(out2(1,k,3),k=1,kdimprs)=',
     $       (out2(1,k,3),k=1,kdimprs)
        print *,'t :(fsig(1,2+k),k=1,kdim)=',
     $       (fsig(1,2+k),k=1,kdim)
        print *,'rh:(out2(1,k,4),k=1,kdimprs)=',
     $       (out2(1,k,4),k=1,kdimprs)
        print *,'rh:(fsig(1,2+3*kdim+k),k=1,kdim)=',
     $       (fsig(1,2+3*kdim+k),k=1,kdim)
        print *,'u :(out2(1,k,1),k=1,kdimprs)=',
     $       (out2(1,k,1),k=1,kdimprs)
        print *,'u :(fsig(1,2+1*kdim+k),k=1,kdim)=',
     $       (fsig(1,2+1*kdim+k),k=1,kdim)
        print *,'v :(out2(1,k,2),k=1,kdimprs)=',
     $       (out2(1,k,2),k=1,kdimprs)
        print *,'v :(fsig(1,2+2*kdim+k),k=1,kdim)=',
     $       (fsig(1,2+2*kdim+k),k=1,kdim)
        
c        do i=1,idimt
c           if (out2(i,2,5).gt.4000.) then
c              print *,'lat and i:',lat,i
c              print *,'zsin and zsout:',
c     $             out2(i,2,5),fsig(i,1)
c              print *,'psin and psout:',
c     $             exp(ps(i))*10.,exp(fsig(i,2))*10.
c              do k=1,kdimprs
c                 print *,'in :',slin(i,k),out2(i,k,3)
c              enddo
c              do k=1,kdim
c                 print *,'out:',sl(k),fsig(i,2+k)
c              enddo
c           endif
c        enddo
#endif
#endif
        if (isfcp.eq.0) then
ccc extrapolation for below old surface (constant lapse rate)
        do i=1,idimt
          ps1=log(slin(i,1))+ps(i)
          do k=1,kdim
            psk=log(sl(k))+fsig(i,2)
            if(psk.ge.ps1) then !! below 1000hPa level
              fsig(i,k+2)=out2(i,1,3)*exp(gamma*rd_/g_*(psk-ps1)) !! t
              fsig(i,3*kdim+k+2)=out2(i,1,4) !! rh
            endif
          enddo
        enddo
        endif
c
c  convert rh to q and temp to virtual temp
c
        do k=1,kdim
          do i=1,idimt
c#ifdef ICE
c           es=fpvs(fsig(i,k+2))
c#else
            es=fpvs0(fsig(i,k+2)) !! es is in cb
c#endif
            qs=xeps*es/((exp(fsig(i,2))*sl(k))+xepsm1*es)
c            if (fsig(i,3*kdim+k+2).lt.0.) then
c            if ((k.eq.kdim).and.(i.eq.1)) then
c            if ((i.eq.1)) then
c            print*,'kei',qs,xeps,es,xepsm1,
c     $           fsig(i,3*kdim+k+2),exp(fsig(i,2))*sl(k)
c            endif
            if (inorhv.eq.0) then
               fsig(i,3*kdim+k+2)=fsig(i,3*kdim+k+2)/100.*qs
            endif
            fsig(i,k+2)=fsig(i,k+2)*(1+0.61*fsig(i,3*kdim+k+2))
          enddo
        enddo
c
        if (ihyb.eq.0) then
c
        do k=1,kdimprs
          do i=1,idimt
c#ifdef ICE
c           es=fpvs(fsig(i,k+2))
c#else
            es=fpvs0(out2(i,k,3)) !! es is in cb
c#endif
            qs=xeps*es/(stdprs(k)/10.+xepsm1*es)
            if (inorhv.eq.0) then
               out2(i,k,4)=out2(i,k,4)/100.*qs
            endif
            out2(i,k,3)=out2(i,k,3)*(1+0.61*out2(i,k,4))
          enddo
        enddo
c
        elseif (ihyb.eq.1) then
c
        do k=1,kdimprs
          do i=1,idimt
c#ifdef ICE
c           es=fpvs(fsig(i,k+2))
c#else
            es=fpvs0(out2(i,k,3)) !! es is in cb
c#endif
            ps1=hyal(int(stdprs(k)))*1000.+
     $           hybl(int(stdprs(k)))*exp(ps(i))*10.
            qs=xeps*es/(ps1/10.+xepsm1*es)
            if (inorhv.eq.0) then
               out2(i,k,4)=out2(i,k,4)/100.*qs
            endif
            out2(i,k,3)=out2(i,k,3)*(1+0.61*out2(i,k,4))
          enddo
        enddo
c
        elseif (ihyb.eq.2) then
c
        do k=1,kdimprs
          do i=1,idimt
c#ifdef ICE
c           es=fpvs(fsig(i,k+2))
c#else
            es=fpvs0(out2(i,k,3)) !! es is in cb
c#endif
            ps1=stdprs(k)/10000.*exp(ps(i))*10.
            qs=xeps*es/(ps1/10.+xepsm1*es)
            if (inorhv.eq.0) then
               out2(i,k,4)=out2(i,k,4)/100.*qs
            endif
            out2(i,k,3)=out2(i,k,3)*(1+0.61*out2(i,k,4))
          enddo
        enddo
c
        endif
c
        print *,'Tv:(out2(1,k,3),k=1,kdimprs)=',
     $       (out2(1,k,3),k=1,kdimprs)
        print *,'Tv:(fsig(1,2+k),k=1,kdim)=',
     $       (fsig(1,2+k),k=1,kdim)
        print *,'q :(out2(1,k,4),k=1,kdimprs)=',
     $       (out2(1,k,4),k=1,kdimprs)
        print *,'q :(fsig(1,2+3*kdim+k),k=1,kdim)=',
     $       (fsig(1,2+3*kdim+k),k=1,kdim)
c
c
c  extrap q to stratosphere if needed 
c
#ifdef G2G
        do i=1,idimt
          xlat(i)=pi_/2.0-colrad(lat)
        enddo
        do i=1,idimt
          expps(i)=exp(fsig(i,2))
c          expps(i)=exp(ps(i))
        enddo
c        call exmois2(sl,expps,fsig(1,3+3*kdim),fsig(1,3),
c     1               fsig(1,3+3*kdim),idimt,kdim,kdim,xlat)
cc     1               qe(1,1),idimt,kdim,kdim,xlat)
cc        call exmois2(sl,expps,out2(1,1,4),out2(1,1,3),
cc     1               fsig(1,3+3*kdim),idimt,kdim,kdimprs,xlat)
c        print *,'qe:(fsig(1,2+3*kdim+k),k=1,kdim)=',
c     $       (fsig(1,2+3*kdim+k),k=1,kdim)
cc     $       (qe(1,k),k=1,kdim)
c
#ifdef DBG
c
c        print *,'gz:(fsig(i,1),i=1,10)=',(fsig(i,1),i=1,10)
c        print *,'lnps:(fsig(i,2),i=1,10)=',(fsig(i,2),i=1,10)
c        print *,'t:(fsig(i,3+20),i=1,10)=',
c     $       (fsig(i,3+20),i=1,10)
c        print *,'q:(fsig(i,3+3*kdim+20),i=1,10)=',
c     1       (fsig(i,3+3*kdim+20),i=1,10)
c        print *,'u:(fsig(i,3+kdim*1+20),i=1,10)=',
c     $       (fsig(i,3+kdim*2+20),i=1,10)
c        print *,'v:(fsig(i,3+kdim*2+20),i=1,10)=',
c     $       (fsig(i,3+kdim*3+20),i=1,10)
c
        do i=1,idimt
         zout(i,lat)=fsig(i,1)
         pout(i,lat)=fsig(i,2)
        enddo
        do k=1,kdim
          do i=1,idimt
            tout(i,lat,k)=fsig(i,2+kdim*0+k)
            uout(i,lat,k)=fsig(i,2+kdim*1+k)
            vout(i,lat,k)=fsig(i,2+kdim*2+k)
            qout(i,lat,k)=fsig(i,2+kdim*3+k)
          enddo
        enddo
        if (ntotal_.gt.1) then
          do k=kdim+1,kdimh
            do i=1,idimt
              qout(i,lat,k)=fsig(i,2+kdim*3+k)
            enddo
          enddo
        endif
#endif
#endif
c
#ifdef G2G
c  gsm output
c  spherical transform
c
        call fftidim(fsig,zss,mkd,-1)
c
        do i=1,mdim
          qdert(i) = qlnt(i) * wgt(lat)
        enddo
c        print*,'qdert',qdert        
c
        call symasa(zss,zsa,1)
        call fl22a(zss,zsa,gzo,qdert,1)
c
        call symasa(pss,psa,1)
        call fl22a(pss,psa,qo ,qdert,1)
c
        call symasa(tts,tta,kdim)
        call fl22a(tts,tta,teo,qdert,kdim)
c
        call symasa(qqs,qqa,kdimh)
        call fl22a(qqs,qqa,rqo,qdert,kdimh)
c
        call gozrmt(qlnt,qlnv,qdert,eps,lat,qlnwct,rcs2,wgt(lat))
c        print*,'qlnwct',qlnwct
c
        call symasa(uus,uua,kdim)
        call symasa(vvs,vva,kdim)
        call psum2a(uua,uus,vva,vvs,zeo,qlnwct,qdert,kdim)
        call msum2a(uua,uus,vva,vvs,dio,qlnwct,qdert,kdim)
c
c        print *,'u :(fsig(i,2+1*kdim+1),i=1,idimt)=',
c     $       (fsig(i,2+1*kdim+1),i=1,idimt)
c        print *,'us:(uus(i,1),i=1,idimt)=',
c     $       (uus(i,1),i=1,idimt)
c        print *,'ua:(uua(i,1),i=1,idimt)=',
c     $       (uua(i,1),i=1,idimt)
#else
!     for regional output   haiqinli
        do k=1,kdim*nfld+2
         do i=1,idim
          rfsig(i,lat,k)=fsig(i,k)
         enddo
        enddo 
#endif
      enddo
!     end of the long lat loop
#ifdef G2G
      do k=1,kdim
        do m=1,mdim
          dio(m,k)=-dio(m,k)
        enddo
      enddo
c
#ifdef DBG
c
      open(53,form='unformatted',access='direct',
     $     recl=4*idim*jdim)
      tt=0
      call rmaxmin(zout,idim*jdim,1,1,1,'zout')
      call rmaxmin(pout,idim*jdim,1,1,1,'pout')
      call rmaxmin(tout,idim*jdim,kdim,1,kdim,'tout')
      call rmaxmin(qout,idim*jdim,kdim,1,kdimh,'qout')
      call rmaxmin(uout,idim*jdim,kdim,1,kdim,'uout')
      call rmaxmin(vout,idim*jdim,kdim,1,kdim,'vout')
c
      call rowsep2(zout)
      call rowsep2(pout)
      do k=1,kdim
        call rowsep2(tout(1,1,k))
        call rowsep2(uout(1,1,k))
        call rowsep2(vout(1,1,k))
      enddo
      do k=1,kdimh
        call rowsep2(qout(1,1,k))
      enddo
      do ij=1,idim*jdim
        outss(ij,1)=zout(ij,1)
        if ((outss(ij,1).gt.-1.E32).and.(outss(ij,1).lt.1.E32)) then
           continue
        else
           print*,'strange data z',ij
        endif
      enddo
      tt=tt+1
      write(53,rec=tt) (outss(ij,1),ij=1,idim*jdim)
      do ij=1,idim*jdim
        outss(ij,1)=pout(ij,1)
        if ((outss(ij,1).gt.-1.E32).and.(outss(ij,1).lt.1.E32)) then
           continue
        else
           print*,'strange data p',ij
        endif
      enddo
      tt=tt+1
      write(53,rec=tt) (outss(ij,1),ij=1,idim*jdim)
      do k=1,kdim
        do ij=1,idim*jdim
          outss(ij,1)=tout(ij,1,k)
        if ((outss(ij,1).gt.-1.E32).and.(outss(ij,1).lt.1.E32)) then
           continue
        else
           print*,'strange data t',ij
        endif
        enddo
      tt=tt+1
      write(53,rec=tt) (outss(ij,1),ij=1,idim*jdim)
      enddo
      do k=1,kdimh
        do ij=1,idim*jdim
          outss(ij,1)=qout(ij,1,k)
        if ((outss(ij,1).gt.-1.E32).and.(outss(ij,1).lt.1.E32)) then
           continue
        else
           print*,'strange data q',ij
        endif
        enddo
      tt=tt+1
      write(53,rec=tt) (outss(ij,1),ij=1,idim*jdim)
      enddo
      do k=1,kdim
        do ij=1,idim*jdim
          outss(ij,1)=uout(ij,1,k)
        if ((outss(ij,1).gt.-1.E32).and.(outss(ij,1).lt.1.E32)) then
           continue
        else
           print*,'strange data u',ij
        endif
        enddo
      tt=tt+1
      write(53,rec=tt) (outss(ij,1),ij=1,idim*jdim)
      enddo
      do k=1,kdim
        do ij=1,idim*jdim
          outss(ij,1)=vout(ij,1,k)
        if ((outss(ij,1).gt.-1.E32).and.(outss(ij,1).lt.1.E32)) then
           continue
        else
           print*,'strange data v',ij
        endif
        enddo
      tt=tt+1
      write(53,rec=tt) (outss(ij,1),ij=1,idim*jdim)
      enddo
#endif
#endif
c
      idate(4)=iy
      idate(2)=im
      idate(3)=id
      idate(1)=ih
      fhour=fh
      lab='        '
c
#ifdef G2G
c      call rmsgto(qo,dio,teo,dio,del,rqo)
      call rmsgto(qo,dio,teo,zeo,del,rqo)
      do k=1,kdum
        dummy(k)=0.
      enddo
      waves=jcap_
      xlayers=levs_
      trun=1.
      order=2.
      realform=1.
      gencode=7.
      rlond=lonf_
      rlatd=latg_
      rlonp=lonf_
      rlatp=latg_
      rlonr=lonf_
      rlatr=latg_
c      tracers=1.
c      cloud=0.
      tracers=ntrac_
      cloud=ncldg_
      pdryini=1.
      subcen=0.
      ensemble(1)=0.
      ensemble(2)=0.
      ppid=0.
      slid=0.
      vcid=0.
      vmid=0.
      vtid=0.
      runid=0.
      usrid=0.
      do k=1,kdum2
        dummy2(k)=0.
      enddo
c      do k=1,kdim
c        do m=1,mdim
c          dio(m,k)=-dio(m,k)
c        enddo
c      enddo
      write(51) lab
      write(51) fhour,idate,si,sl
     &       ,dummy,waves,xlayers,trun,order,realform,gencode
     &       ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers
     &       ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid
     &       ,pdryini,dummy2,cloud
      do i = 1, mdim
        b(indxmm(i)) = gzo(i)
      enddo
      write(51) b
      do i = 1, mdim
        b(indxmm(i)) = qo(i)
      enddo
      write(51) b
      do k=1,kdim
        do i = 1, mdim
          b(indxmm(i)) = teo(i,k)
        enddo
        write(51) b
      enddo
      do k=1,kdim
        do i = 1, mdim
          b(indxmm(i)) = dio(i,k)
        enddo
        write(51) b
        do i = 1, mdim
          b(indxmm(i)) = zeo(i,k)
        enddo
        write(51) b
      enddo
      do k=1,kdimh
        do i = 1, mdim
          b(indxmm(i)) = rqo(i,k)
        enddo
        write(51) b
      enddo
c
#else
      do i=1,2*levmax+1-levs_-1-levs_
      dummy(i)=0.0
      enddo
      do i=1,512-(6+2*levmax)
      ext(i)=0.0
      enddo
c
c     to calculate rcenlat,rcenlat,rlftgrd,rbtmgrd (haiqin li)
c     now only support Mecator & Polar Stereographic projection
      if(rproj.eq.0) then ! Mecator projection  
        rlftgrd=idim/2
        rlftgrd=jdim/2
        truth=rtruth*rad
        rlon1=rlon1*rad
        rlon1=mod(rlon1,twopi)
        rlat1=rlat1*rad
        a2=rerth_*cos(truth)
        x00=-(rlftgrd-1.)*rdelx - 0.5*rdelx
        cenlon=rlon1-x00/a2
        rcenlon=cenlon/rad
        if(rcenlon.gt.180) rcenlon=rcenlon-360.

        y1=(log(tan((rlat1/2.)+qtpi)))*a2
        y00=-(rbtmgrd-1.)*rdely-0.5*rdely
        y0=y1-y00
        cenlat=(atan(exp(y0/a2))-qtpi)*2.
        cenlat=mod(cenlat,pi_)
        rcenlat=cenlat/rad
        print*,'Mecator rcenlat,rcenlon',rcenlat,rcenlon
      elseif(abs(rproj).eq.1.) then ! Polar stereographic projection
        if(rproj.eq.1) then
           rcenlat=90.
           rtruth=60.
        elseif(rproj.eq.-1.) then
           rcenlat=-90.
           rtruth=-60.
        endif
        rcenlon=0.
        print*,'rcenlat,rcenlon,rorient',rcenlat,rcenlon,rorient
        rlat1=rlat1*rad
        rlon1=rlon1*rad
        truth  = rtruth * rad
        truth  = rproj * truth
        orient  = rorient * rad
        dlamda0 = orient + hfpi
        print*,'dlamda0=',dlamda0
        a2 =  rerth_ * ( 1.0 + sin(truth) )
        cenlat = rcenlat
        rlat1=rlat1/rproj
        rsoa2=tan((hfpi-rlat1)/2.0e0)
        rlon1=(rlon1-dlamda0)/rproj
        x00=((rsoa2*a2)**2/(1+(tan(rlon1-pi))**2))
        x00=-sqrt(x00)
!   (add a small number for integer truncation)
        y00=x00*tan(rlon1)-100.
        print*,'x00=',x00
        print*,'y00=',y00
        rlftgrd=1-x00/rdelx
        rbtmgrd=abs(1-y00/rdely)
        print*,'polar rlftgrd=',rlftgrd
        print*,'polar rbtmgrd=',rbtmgrd
      endif

      ext(1)=iwav1_
      ext(2)=jwav1_
      ext(3)=igrd1_
      ext(4)=jgrd1_
      ext(5)=levs_
      ext(6)=2+levs_*4+5
      ext(7)=rproj             !projection flag
      ext(8)=rtruth
      ext(9)=rorient
      ext(10)=rcenlat
      ext(11)=rcenlon
      ext(12)=rlftgrd
      ext(13)=rbtmgrd
      ext(14)=rdelx
      ext(15)=rdely
      ext(19)=ntrac_
      ext(20)=ncldg_
c     begin to calculate flat & flon  (haiqin li)
      call setgrd(flat,flon,rdelx,rdely,dlamda0)     

!     begin to calculate map factor (haiqin li)
c
c define regional map factors
c
c      rad = pi_ / 180.0
c      hfpi= 0.5 * pi_
      truth = rtruth
c      nproj = rproj
      if( rproj.eq.1 .or. rproj.eq.-1 ) then
        orient = rorient
        dlamda0 = orient * rad + hfpi
        a2 = ( 1.0 + sin( rproj * truth * rad ) )
        do n=1,idim*jdim
          fm2(n)=(1.+sin( rproj * flat(n) ))/a2  ! 1/(map factor)
          rgln = flon(n) - dlamda0
          rgln = rproj * rgln
          rcsln(n) = rproj * cos( rgln )
          rsnln(n) = rproj * sin( rgln )
        enddo
      else if( rproj.eq.0 ) then
        a2 = cos( truth * rad )
        do n=1,idim*jdim
          fm2(n)=cos( flat(n) )/a2               ! 1/(map factor)
          rcsln(n) =  0.0
          rsnln(n) = -1.0
        enddo
      endif
c
      call shalfo(fm2,1)
      call todxdy(fm2,fm2x,fm2y,igrd1_,jgrd1_,rdelx,rdely)
      call shalfi(fm2,1)
      call shalfi(fm2x,1)
      call shalfi(fm2y,1)
c
      do ij=1,idim*jdim
        fm2(ij) = 1.0/ ( fm2(ij) * fm2(ij) )
        fm2x(ij) = -2.0 * fm2(ij)**1.5 * fm2x(ij)
        fm2y(ij) = -2.0 * fm2(ij)**1.5 * fm2y(ij)
      enddo
      call rmaxmin(fm2,idim*jdim,1,1,1,'fm2')
      call rmaxmin(fm2x,idim*jdim,1,1,1,'fm2x')
      call rmaxmin(fm2y,idim*jdim,1,1,1,'fm2y')
c
      call shalfo(flat, 1)
      call shalfo(flon, 1)
      call shalfo(fm2,1)
      call shalfo(fm2x,1)
      call shalfo(fm2y,1)
c
      print*,'fm2(1)',fm2(1),fm2(2),fm2(3)
      print*,'fm2x(1)',fm2x(1)
      print*,'fm2y(1)',fm2y(1)
      print*,'flat(1)/rad',flat(1)/rad
      print*,'flon(1)/rad',flon(1)/rad
      print*,'flat(idim*jdim)/rad',flat(idim*jdim)/rad
      print*,'flon(idim*jdim)/rad',flon(idim*jdim)/rad
c     begin to write regional out (haiqin li)

      write(51) lab
      write(51) fhour,idate,si,sl,
     &    (dummy(i),i=1,2*levmax+1-levs_-1-levs_),ext
      write(51)((rfsig(i,j,1),i=1,idim),j=1,jdim)  ! gz
      write(51)((rfsig(i,j,2),i=1,idim),j=1,jdim)  ! ps
      do k=3,2+kdim
      write(51)((rfsig(i,j,k),i=1,idim),j=1,jdim)  ! t
      enddo
      do k=1,kdim
       write(51)((rfsig(i,j,2+kdim*1+k),i=1,idim),j=1,jdim) ! u
       write(51)((rfsig(i,j,2+kdim*2+k),i=1,idim),j=1,jdim) ! v
      enddo
      do k=3+3*kdim,2+kdim*4
      write(51)((rfsig(i,j,k),i=1,idim),j=1,jdim)  ! q
      enddo
      write(51)fm2
      write(51)fm2x
      write(51)fm2y
      write(51)flat
      write(51)flon
#endif

      end
