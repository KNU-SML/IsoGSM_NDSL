#define IO2 io2
#define JO2 jo2
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JCAPP jcapp
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LCAPP lcapp
#define LEVH levh
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LNUV lnuv
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LPNT lpnt
#define LSOIL lsoil
#define LTSTP ltstp
#define MLVARK mlvark
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NPES npes
#define NROW nrow
#define NSOIL nsoil
#define NTRAC ntrac
#define SLVARK slvark
#define TWOJ1 twoj1


        subroutine gsmstep
c$$$  subroutine program documentation block
c
c subroutine: gsmstep        make global forecast with spectral model
c   prgmmr: juang            org: w/np5       date: 00-06-01
c
c abstract: make global forecast step with spectral model.
c           this routine combines forward-step (stepone) and
c           full-time leap-frog time step together.
c
c program history log:
c   81-01-01  sela
c   99-05-21  hann-ming henry juang   add mpi in one loop version
c   99-06-22  hann-ming henry juang   add mpi in all loop version
c   00-06-01  hann-ming henry juang   change gsm into gsmini, gsmstep, and gsmend
c
c attributes:
c   language: fortran
c
c$$$
c
#include <paramodel.h>
c include comfibm inludes other includes ... i.e. geshem runoff
#include <comfibm.h>
#include <comcon.h>
#include <comgpd.h>
#include <comfgsm.h>
#ifdef LFM
#include <comlfm.h>
#endif
c
#ifdef SCM
      common /comscmbdy/dtbdy, curtime, vvel(lnt22_,levs_)
#ifdef CLM_CWF
     +                , wdiv(lnt22_,levs_), hadq(lnt22_,levs_)
#endif
#endif
c
#ifdef MP
#include <commpi.h>
#endif
#include <sfcfld.h>
c
#ifdef ISOTOPE
#include <comfiso.h>
      dimension rresv(lonf2_,latg2_,lsoil_,2)
#endif
#ifdef NISLQ
#include <comnislq.h>
#endif
c
      character*128 fno
      integer ncho
c
#ifdef MP
#define LNT22S lnt22p_
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LNT22S lnt22_
#define LONF2S lonf2_
#define LATG2S latg2_
#endif

#ifdef SAS_DIAG
      common/comsdiag/dcu(LONF2S,levs_,latg2_),dcv(LONF2S,levs_,latg2_),
     +                dct(LONF2S,levs_,latg2_),dcq(LONF2S,levs_,latg2_),
     +                dch(LONF2S,levs_,latg2_),
     +                fcu(LONF2S,levs_,latg2_),fcd(LONF2S,levs_,latg2_),
     +                deltb(LONF2S,latg2_),delqb(LONF2S,latg2_),
     +                delhb(LONF2S,latg2_),cbmf(LONF2S,latg2_),
     +                dlt(LONF2S,levs_,latg2_),dlq(LONF2S,levs_,latg2_),
     +                dlh(LONF2S,levs_,latg2_)
      dimension cvbp(LONF2S,latg2_), cvtp(LONF2S,latg2_)
#endif
#ifdef CONDENS
c ...............................................................
      common/comcond/
     +       qcicnv(LONF2S,levs_,latg2_), qrscnv(LONF2S,levs_,latg2_),
     +       qcilrg(LONF2S,levs_,latg2_), qrslrg(LONF2S,levs_,latg2_),
     +       taucld(LONF2S,levs_,latg2_),
     +       cldwp(LONF2S,levs_,latg2_),  cldip(LONF2S,levs_,latg2_)
#endif
#ifdef CLM_CWF
      parameter (npmx=1000)
      common /comfcwrk/cgs(LONF2S,levs_,latg2_),
     +                 cwgt(npmx,lonf2_,latg2_)
      common /comfcwrki/inx(npmx,lonf2_,latg2_),jnx(npmx,lonf2_,latg2_),
     +                 npnt(lonf2_,latg2_), lwconv
      logical lwconv
#ifdef DBG
      print *, 'CWRK-lwconv', lwconv
#endif
#endif
#ifdef CO2PROG
** co2 **
      dimension
     *     out(lonf_*latg_),dummy2(lonf_,latg_),
     *     data(720*360),dummy1(720,360)
      common /co2/co2fl(LONF2S,LATG2S)
      integer  kpd(25),imdata,jmdata,imaxgrb,jmaxgrb,
     $     kpds5,kgds1,iret,kgau,idim,jdim,ijdim
      character*128 fngrib
#endif

       real , dimension(LONF2S,levh_,LATG2S) :: q2_rev
c
#ifdef RIVER
#include <comfriv.h>
#endif
c
c  start one time step
c
      if( stepone ) then
        if (jdt.gt.limlow) inistp=0
#ifdef KEN
        if (jdt.gt.limlow) isave = 0
#endif
        kdt = jdt - limlow + 1
        lastep=jdt.eq.1
      else
        kdt=jdt
        lastep=jdt.eq.maxstp
#ifdef KEN
        if(npoint.gt.0.and.itnum.gt.nstken.and.iope) then
#ifndef NOPRINT
          if(iope)
     &     print *,'ken points disabled - time levels exceed ',nstken
#endif
          npoint=0
        endif
#endif
      endif
c 
c  The following print is kept
c
#ifndef NOPRINT
      if(iope) then
        if(stepone) then
          print 81,kdt
        else
          print 82,kdt
        endif
      endif
#else
      if(iope) then
        if(mod(kdt-1,10).eq.0) then
          if(stepone) then
            print 81,kdt
          else
            print 82,kdt
          endif
        endif
      endif
#endif
81      format(1h ,'forward step: kdt in gsmstep=',i3)
82      format(1h ,'full step   : kdt in gsmstep=',i3)
c 
      if( .not. stepone ) then
c
c   set switch for saving kuo data (for interactive clouds)..
c
        cvmod= mod (solhr+dthr,dtcvav)
        if(cvmod.lt.hdthr.or.cvmod.ge.dtcvav-hdthr) then
          clstp=min(dtcvav,(shour+deltim)/3600.)
        elseif(clstp.gt.0.) then
          clstp=0.
        else
          clstp=-10.
        endif
      endif
c
#ifdef ONELOOP
c
c  one loop includes gloopr, gloopa, and gloopb.
c
      call gloop
#ifdef DBG
      if(iope) then
        if(stepone) then
          print *,' forward time step:  after gloop '
        else
          print *,' full time step: after gloop '
        endif
      endif
      call rmsgt(z ,x ,y ,w ,del,rt)
#endif
c
#else
c
#ifndef NO_PHYSICS
c
c  radiation loop
c
      if( jdt.eq.limlow .or. .not. stepone ) then
        call gloopr
#ifdef DBG
        if(iope) then
          if(stepone) then
            print *,' forward time step:  after gloopr '
          else
            print *,' full time step: after gloopr '
          endif
        endif
        call rmsgt(q,di,te,ze,del,rq)
#endif
      endif
#ifdef PSPLIT
c
c process splitting physics tendency initialization
c
      if( jdt.eq.limlow) then
        call gloopb
      endif
#endif
#endif
c
c  dynamic loop
c
      call gloopa
#ifdef DBG
      if(iope) then
        if(stepone) then
          print *,' forward time step:  after gloopa '
        else
          print *,' full time step: after gloopa '
        endif
      endif
      call rmsgt(z ,x  ,y  ,w  ,del,rt)
#endif
c
#endif
c
c advance shour (second) and thour (total hour)
c
      if(stepone) then
        shour=deltim
        dtpost=deltim
      endif
      shour=shour+deltim
      dtpost=dtpost+deltim
      if(.not.stepone) thour=fhour+ifix(shour/3600.+0.5e0)
#ifdef SCM
      curtime = fhour + shour/3600.
#endif
c
#ifndef SCM
c
c  semi-implicit
c
#ifndef HYBRID
      call sicdif(dim,tem,qm,x,y,z,uln,vln,lnts2,lnoffset)
#else
      call sicdif_hyb(dim,tem,qm,x,y,z,di,te,q,uln,vln,lnts2,lnoffset)
#endif

#ifdef DBG
      if(iope) then
        if(stepone) then
          print *,' forward time step:  after sicdif '
        else
          print *,' full time step: after sicdif '
        endif
      endif
      call rmsgt(z ,x  ,y  ,w  ,del,rt)
#endif
c
c  horizontal diffusion
c
      call deldif(rt,w,deltim,qm,sl,x,y,lnts2,lnoffset)
#ifdef DBG
      if(iope) then
        if(stepone) then
          print *,' forward time step:  after deldif '
        else
          print *,' full time step: after deldif '
        endif
      endif
      call rmsgt(z ,x  ,y  ,w  ,del,rt)
#endif

#ifdef NISLQ
c
c forward semi-Largrangian advection for specific humidity (t+dt)
c
      do j = 1, LATG2S
        do k = 1, levh_
          do i = 1, LONF2S
            q2_rev(i,k,j) = slq_q1(i,levh_+1-k,j)
          enddo
        enddo
      enddo

#ifdef NISLQ_MASS
      call nislq_mass_advect(deltim,slq_psfc2,slq_u2,slq_v2,
     &                       slq_w2,q2_rev,slq_q3)
#else
      call nislq_mono_advect(deltim,slq_psfc2,slq_u2,slq_v2,
     &                       slq_w2,slq_q1,slq_q3)
#endif

#ifdef NISLQ_MONOMASS
c
c 3D mass conservation for nislq
c
      if(kdt.ne.1) call nislq_mass_adjustment(slq_q1,slq_q3,qm,z)
#endif /* NISLQ_MONOMASS */

#endif /* NISLQ end */

#endif /* ~SCM end */
c
c  partial time filter to n and save in n-1
c

#ifdef NISLQ
#define RM slq_q1
#define RQ slq_q2
#else
#define RM rm
#define RQ rq
#endif
      if( .not. stepone ) then
        call filtr1(tem,te,dim,di,zem,ze,RM,RQ,qm,q,
     1                       filta,LNT22S,lnts2,levs_,levh_)
#ifdef DBG
        if(iope) then
          print *,' full time step: after filtr1 '
        endif
        call rmsgt(z ,x  ,y  ,w  ,del,rt)
#endif
      endif
      do j=1,lnts2
        q(j)=z(j)
      enddo
c
#ifndef ONELOOP
#ifndef PSPLIT
c
c  physics loop
c

#ifndef NO_PHYSICS

        call gloopb
c      slq_qtemp = slq_q3 
#ifdef DBG
      if(iope) then
        if(stepone) then
          print *,' forward time step:  after gloopb '
        else
          print *,' full time step: after gloopb '
        endif
      endif
      call rmsgt(q ,x  ,y  ,w  ,del,rt)
#endif

#endif
c
#endif
#endif
#ifdef RIVER
c
c  river discharge (kei)
c
      call pgriv
#endif
c
#ifndef SCM
c  zonal wind damping
c
      call dampux(x,w,y,rt,deltim,uln,vln,spdmax,
     1                       snnp1,lnts2,lnoffset)
#ifdef DBG
      if(iope) then
        if(stepone) then
          print *,' forward time step:  after dampux '
        else
          print *,' full time step: after dampux '
        endif
      endif
      call rmsgt(q ,x  ,y  ,w  ,del,rt)
#endif
#endif
c
c  update n by n+1
c
      do k=1,levs_
        do j=1,lnts2
          di(j,k)=x(j,k)
          ze(j,k)=w(j,k)
          te(j,k)=y(j,k)
        enddo
      enddo
#ifdef NISLQ 
      do j = 1, LATG2S
        do k = 1, levh_
          do i = 1, LONF2S
            slq_q2(i,k,j)=slq_q3(i,k,j)
          enddo
        enddo
      enddo
#else
      do k=1,levh_
        do j=1,lnts2
          rq(j,k)=rt(j,k)
        enddo
      enddo
#endif 
c
c  complete time filter for n and save to n-1
c
      if( .not. stepone ) then
#ifdef NISLQ
#define RM slq_q1
#define RQ slq_q2
#endif
        call filtr2(tem,te,dim,di,zem,ze,RM,RQ,
     1                       filta,LNT22S,lnts2,levs_,levh_)
c
#ifdef DBG
        if(iope) then
          print *,' full time step: after filtr2 '
        endif
        call rmsgt(q ,di  ,te  ,ze  ,del,rq)
#endif
      endif
#ifdef LFM
c
c  lfm filtering
c
      if( .not. stepone ) then
        ifstep=kdt+ipstep
        call lfmfilt(ifstep,thour)
      endif
#endif
c
#ifdef MP
#ifdef GSMDAMP
      call gloopn
#endif
#endif
c
c update the matrix in case of changing deltim within stepone
c
      if( stepone ) then
        deltim=deltim*2.e0
#ifndef SCM
        call gsicdf(deltim,am,bm,gv,sv,cm)
#ifdef DBG
        if( iope ) print *,' forward time: gsicdf for deltim= ',deltim
#endif
#endif
      endif
c
#ifdef DFI
c
c dfini : call digital filter initialization every step if con(3).gt.0.0
c
      if(con(3).ne.0.0) call dfini(n1,con(3),shour/3600.,solsec,1)
#endif
c
c  advance solhr
c
      if( .not. stepone ) then
        solsec=solsec+deltim
        solhr=solsec/3600.e0
        iday=solhr/24.e0
        solhr=solhr-iday*24.e0
      endif
#ifdef KEN
c
c   for grid point diag advance itnum, if proper time, and set isave
c
      isave = 0
#endif
#ifdef DFI
      if(numsum.lt.0) then
#endif
#ifdef KEN
        if (ikfreq.gt.1) then
          imodk = mod(jdt,ikfreq)
          if (imodk.eq.0) then
            isave = 1
            itnum = itnum + 1
          end if
        else
          isave = 1
          itnum = itnum + 1
        end if
#endif
#ifdef DFI
      end if
#endif
#ifdef RSM
c
c rsm : call main routine of regional forecast
c
      call rsmsmf(fhour,shour,gz,q,te,di,ze,rq)
#endif
#ifndef SCM
c
c  mass conservation, surface merge and output files
c
      if((mod(kdt,max(kpfix,1)).eq.0.or.lastep).and.(kpfix.ne.0)) then
        if( .not. stepone ) then
          call gloopp(q,avprs0)
          call gloopp(qm,avprs0)
        endif
          avprs0=0.
          call gloopp(q,avprs0)
      endif
#endif
#ifdef DBG
        if(iope) then
          print *,' full time step: after gloopp '
        endif
        call rmsgt(q ,di  ,te  ,ze  ,del,rq)
#endif
c 
c write sigma and surface files
c
      if( .not. stepone ) then
c
c write sigit for diagnostics and rsm boundary conditions
c
        if(mod(kdt,max(ksig,1)).eq.0.or.lastep) then
#ifndef HYBRID
          call wrisig(n1,thour,idate,qm,tem,dim,zem,rm,sl,si,gz,z00,
#else
          call wrisig(n1,thour,idate,qm,tem,dim,zem,rm,ak5,bk5,gz,z00,
#endif
#ifdef SCM
     1            vvel,
#endif
     1            1)
#ifdef CLM_CWF
          print *, 'RDWRCGS write kdt=',kdt,' thour=',thour
          call rdwrcgs(lugi+10,thour,idate,cgs,cwgt,inx,jnx,npnt,1)
#endif
#ifndef NOPRINT
c
c  monitor sigma file
c
          if( iope )
     &    print *,' predicted full values at the end of forecst segment'
          call rmsgt(q,di,te,ze,del,rq)
#endif
        endif
#ifdef ADJSOILW
c
c  for soil adjustment, sfc execution and post processing intervals need
c  to be the same (to match the accumulation period).
c
        kpost=ksfcx
#endif
c
c  change surface field by merging 
c      (not executed at the beginning of time loop).
c
        ifcst=1
#ifndef SKIPSFCMRG
        if(mod(kdt,max(ksfcx,1)).eq.0) then
#ifdef ADJSOILW
#ifdef DBG
          CALL stats(geshem,LONF2S,LATG2S,'gesh')
#endif
          intobsp=ksfcx*con(1)/3600.+0.5
          call sfc(n1,idate(4),idate(2),idate(3),idate(1),thour,ifcst,
     1             intobsp,geshem,runoff)
#else
          call sfc(n1,idate(4),idate(2),idate(3),idate(1),thour,ifcst)
#endif
c
#ifdef CO2PROG
** co2 **
        imdata=720 !! temporary
        jmdata=360 
        ijmdim=imdata*jmdata
        idim=lonf_
        jdim=latg_
        ijdim=idim*jdim
#ifdef MP
      if(iope) then
#endif
        print *,'========== start reading co2 flux file ========='
        do k=1,25
          kpd(k)=-1
        enddo
        kpd(5)=1
        kpd(6)=1
        kpd(7)=0
cc please specify your co2 flux file location cc
        fngrib='/home/kei/data/FLUX_data/transcom_nep_clim.grb'
        call fixrd_clim(n1,fngrib,
     $             kpd,
     &             idate(4),idate(2),idate(3),idate(1),thour,
     &             imdata,jmdata,
     &             imaxgrb,jmaxgrb,
     &             kpds5,kgds1,
     &             proj,orient,truth,cotru,
     &             delx,dely,rlat1,rlat2,rlon1,rlon2,
     &             dlon,dlat,wlon,rnlat,
     &             data,iret)
        if(iret.gt.0) then
          print *,' read failed '
#ifdef MP
          call mpabort
#else
          call abort
#endif
        endif
        do ij=1,imaxgrb*jmaxgrb
          dummy1(ij,1)=1.
        enddo
        do ij=1,idim*jdim
          dummy2(ij,1)=1.
        enddo
        kgau=0
        if(kgds1.eq.4) kgau=1
        call la2ga(data,imaxgrb,jmaxgrb,
     &             abs(dlon),abs(dlat),wlon,rnlat,
     &             out,idim,jdim,.false.,dummy1,dummy2,kgau,
     &             0)
c
#ifdef MP
        do ij=1,idim*jdim
          out(ij)=out(ij)*44./12. !! kgC->kgCO2
        enddo
        call row1ns(out)
      endif 
        call mpgf2p(out,lonf2_,latg2_,co2fl,
     1            LONF2S,LATG2S,1)
#else
        do i=1,LONF2S
          do j=1,LATG2S
            ij=(j-1)*LONF2S+i
            co2fl(i,j)=out(ij)*44./12. !! kgC->kgCO2
          enddo
        enddo
        call row1ns(co2fl)
#endif
#ifdef MP
      if(iope) then
#endif
        print *,'========== end reading co2 flux file ========='
#ifdef MP
      endif
#endif
#endif
        endif
#endif
c
c  write surface file
c
        if(mod(kdt,max(ksfc,1)).eq.0.or.lastep) then
          call fnam('sfc',3,thour,fno,ncho)
          call fixio(n1,fno,sfcftyp,
     1               labs,idate(4),idate(2),idate(3),idate(1),thour,
     2               tsea,LONF2S,LATG2S,1)
#ifdef ISOTOPE
#ifdef MP
          call mpgp2f(rresv1,LONF2S,LATG2S,
     $         rresv(1,1,1,1),lonf2_,latg2_,lsoil_)
          call mpgp2f(rresv2,LONF2S,LATG2S,
     $         rresv(1,1,1,2),lonf2_,latg2_,lsoil_)
          if (iope) then
#endif
            close(95)
            call fnam('rresvi',6,thour,fno,ncho)
            open(95,file=fno,form='unformatted')
            print*,'writing ',fno
            write(95) rresv
            close(95)
#ifdef MP
          endif
#endif          
#endif
#ifdef RIVER
c
c  river file
c
#ifdef MP
          if (iope) then
#endif
            close(94)
            call fnam('rivi',4,thour,fno,ncho)
            open(94,file=fno,form='unformatted')
            print*,'writing ',fno
            write(94)gdriv
#ifdef ISOTOPE
            write(94)gdriv1
            write(94)gdriv2
#endif
            close(94)
#ifdef MP
          endif
#endif
#endif
c
        endif
      endif
c
c  write restart files
c
        if(mod(kdt,max(krestart,1)).eq.0.or.lastep) then
#ifndef HYBRID           
          call wrisig(n1,thour,idate,qm,tem,dim,zem,rm,sl,si,gz,z00,
#else
          call wrisig(n1,thour,idate,qm,tem,dim,zem,rm,ak5,bk5,gz,z00,
#endif
#ifdef SCM
     1            vvel,
#endif
     1            4)
#ifndef HYBRID
          call wrisig(n1,thour,idate,q ,te ,di ,ze ,rq,sl,si,gz,z00,
#else
          call wrisig(n1,thour,idate,q ,te ,di ,ze ,rq,ak5,bk5,gz,z00,
#endif
#ifdef SCM
     1            vvel,
#endif
     1            5)
          fno='sfci '
          call fixio(n1,fno,sfcftyp,
     1               labs,idate(4),idate(2),idate(3),idate(1),thour,
     2               tsea,LONF2S,LATG2S,1)
#ifdef ISOTOPE
#ifdef MP
          call mpgp2f(rresv1,LONF2S,LATG2S,
     $         rresv(1,1,1,1),lonf2_,latg2_,lsoil_)
          call mpgp2f(rresv2,LONF2S,LATG2S,
     $         rresv(1,1,1,2),lonf2_,latg2_,lsoil_)
          if (iope) then
#endif
            close(95)
            fno='rresvi'
            open(95,file=fno,form='unformatted')
            print*,'writing rresvi'
            write(95) rresv
            close(95)
#ifdef MP
          endif
#endif
#endif
c
#ifdef RIVER
c
c  river file
c
#ifdef MP
          if (iope) then
#endif
            close(94)
            open(94,file='rivi',form='unformatted')
            print*,'writing rivi'
            write(94)gdriv
#ifdef ISOTOPE
            write(94)gdriv1
            write(94)gdriv2
#endif
            close(94)
#ifdef MP
          endif
#endif
#endif
c
        endif
#ifdef COUPLE_GSM_AND_RSM
c
c  parallel coupling to rsm
c
        if(mod(kdt,max(kcouple,1)).eq.0) then
          call send_sig(qm,tem,dim,zem,rm)
        endif
#endif
c
c  write diagnostic files
c

      if(mod(kdt,max(kpost,1)).eq.0.or.lastep) then
        if( stepone ) then
          xhour=shour/3600.
        else
          xhour=thour
        endif
#ifndef MRG_POST
        call wriflx(n1)
#else
c
#ifndef HYBRID
        call subpgb(thour,idate,sl,si,
#else
        call subpgb(thour,idate,ak5,bk5,
#endif
     &              gz,q,te,di,ze,rq,colrad,dummy,n1,.false.)
#endif
#ifdef MP
#ifdef COUPLE_GSM_AND_ROMS
        call send_flx_recv_sst
#endif
#endif
c
        if(iope) then
          print *,' full time step: after subpgb '
        endif
        call rmsgt(q ,di  ,te  ,ze  ,del,rq)
c
#ifdef SAS_DIAG
        do i = 1, lonf2_
        do j = 1, latg2_
          if (cvb(i,j).ne.100.) then
            kbcnv = int(cvb(i,j))
            cvbp(i,j) = sl(kbcnv) * psurf(i,j) * 1e+3
          else
            cvbp(i,j) =  -999.0
          endif
          if (cvt(i,j).ne.0.) then
            ktcnv = int(cvt(i,j))
            cvtp(i,j) = sl(ktcnv) * psurf(i,j) * 1e+3
          else
            cvtp(i,j) =  -999.0
          endif
        enddo
        enddo
        call wri_diag(thour,lonf2_,latg2_,levs_,
     +                dcu,dcv,dct,dcq,dch,fcu,fcd,deltb,delqb,delhb,
     +                cbmf,cvbp,cvtp,dlt,dlq,dlh)
#endif
#ifdef DGZ
c
c  create zonal diagnostic fields
c
        call wriznl(n1)
        call zerznl
#endif
#ifdef DG3
c
c  write dg3 file
c
#ifdef CONDENS
        call wri_cldrad (fhour,thour,idate,sl,colrad,
     +                   fluxr,cvavg,qcicnv,qrscnv,
     +                   qcilrg,qrslrg,taucld,cldwp,cldip,
     +                   n1)
#else
        call wridia (slmsk,n1)
#endif
        if( .not. stepone ) call zerdia
#endif
#ifdef KEN
c
c write ken
c
#ifdef MP
        call mpgetken
#endif
        call wriken(npoint,ikfreq,imodk,itnum,svdata,
     1              lab,thour,idate,si,sl,nvrken,nptken,nstken,n1)
#endif
        if( .not. stepone ) then
        call zerflx(0)
        dtpost=0.
        endif
        call zertmm
      endif
c
#ifdef LFM
      if(mod(kdt,max(klfm,1)).eq.0.or.lastep) then
        call lfmfout(ifstep,thour,n1)
      endif
#endif
#ifdef RSM
c
c  save rsm out
c
      if(mod(kdt,max(krsm,1)).eq.0.or.lastep) then
        call rsmsav(fhour)
      endif
#endif
c
c  reset timers
c
      if(.not.stepone ) then
        if(mod(kdt,max(min(ksfcx,kpost),1)).eq.0) then
          fhour=thour
          solhr=fhour+idate(1)
          iday=solhr/24.e0
          solhr=solhr-iday*24.e0
          solsec=solhr*3600.
          shour=0.
#ifndef NOPRINT
          if(iope) print *,'reset fhour=',fhour
#endif
        endif
      endif
c
c end of stepone
c
!hoon       if(stepone) then
!hoon         call zerflx(1)
!hoon         dtpost=0.
!hoon       endif
      if( lastep .and. stepone ) then
        stepone=.false.
        solsec=solsec+deltim
        solhr=solsec/3600.0
        deltim=con(1)
#ifdef DG3
        call zerdia
#endif
#ifdef KEN
        call zertmm
#endif
#ifdef DGZ
        call zerznl
#endif
#ifdef LFM
        ifstep=2
        call lfmfilt(ifstep,fhour)
        ipstep=1
#endif
#ifndef NOPRINT
        if(iope) print *,' forward step done. '
        call rmsgt(q, di, te, ze, del, rq)
#endif
      endif
c
c  end  time  step
c
      return
      end
