#define NTOTAL ntotal
#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP1 jcap1
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

        subroutine gsmini
c$$$  subroutine documentation block
c
c subroutine: gsmini         make global initilization with spectral model
c   prgmmr: juang             org: w/np5       date: 00-06-01
c
c abstract: make global forecast with spectral model.
c
c program history log:
c   81-01-01  sela
c   99-05-21  hann-ming henry juang   add mpi in one loop version
c   99-06-22  hann-ming henry juang   add mpi in all loop version
c   00-06-01  hann-ming henry juang   change gsm into gsmini and gsmstep gsmend
c
c attributes:
c   language: fortran
c
c$$$
c
#include <paramodel.h>
#include <comfibm.h>
c
#include <comcon.h>
#include <comgpd.h>
#include <comfgsm.h>
#include <sfcfld.h>
c
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
#ifdef ISOTOPE
#include <comfiso.h>
      dimension rresv(lonf2_,latg2_,lsoil_,2)
      parameter (niso=lonf2_*latg2_*lsoil_*2)
      data      rresv /niso*1./
#endif
#ifdef NISLQ
#include <comnislq.h>
#endif
c
#ifdef MP
#ifdef RSM
#define MPABORT rmpabort
#else
#define MPABORT mpabort
#endif
#else
#define MPABORT abort
#endif
c
#ifdef CONDENS
c ...............................................................
      common/comcond/
     +       qcicnv(LONF2S,levs_,latg2_), qrscnv(LONF2S,levs_,latg2_),
     +       qcilrg(LONF2S,levs_,latg2_), qrslrg(LONF2S,levs_,latg2_),
     +       taucld(LONF2S,levs_,latg2_),
     +       cldwp(LONF2S,levs_,latg2_),  cldip(LONF2S,levs_,latg2_)
#endif
#ifdef SAS_DIAG
c ...............................................................
      common/comsdiag/dcu(LONF2S,levs_,latg2_),dcv(LONF2S,levs_,latg2_),
     +                dct(LONF2S,levs_,latg2_),dcq(LONF2S,levs_,latg2_),
     +                dch(LONF2S,levs_,latg2_),
     +                fcu(LONF2S,levs_,latg2_),fcd(LONF2S,levs_,latg2_),
     +                deltb(LONF2S,latg2_),delqb(LONF2S,latg2_),
     +                delhb(LONF2S,latg2_),cbmf(LONF2S,latg2_),
     +                dlt(LONF2S,levs_,latg2_),dlq(LONF2S,levs_,latg2_),
     +                dlh(LONF2S,levs_,latg2_)
#endif
c
      character*128 fno
c
#ifdef CLM_CWF
      parameter (npmx=1000)
      common /comfcwrk/cgs(LONF2S,levs_,latg2_),
     +                 cwgt(npmx,lonf2_,latg2_)
      common /comfcwrki/inx(npmx,lonf2_,latg2_),jnx(npmx,lonf2_,latg2_),
     +                 npnt(lonf2_,latg2_), lwconv
      logical lwconv
      lwconv = .true.
#endif
c
#ifndef OMP
       itread=0
#endif
#ifdef MP
      call mpdimset(jcap_,levs_,lonf_,latg_)
      iope=mype.eq.master
      lnts2=lntlen(mype)*2
      lons2=lonlen(mype)*2
      lats2=latlen(mype)
      lnoffset=lntstr(mype)*2
#else
      do l=1,jcap1_
        lwvdef(l)=l-1
      enddo
      do l=1,latg2_
        latdef(l)=l
      enddo
      iope=.true.
      lnts2=lnt2_
      lons2=lonf2_
      lats2=latg2_
      lnoffset=0
#endif
c
      call getcon(n1,
#ifdef RSM
     3 nrsmi1,nrsmi2,nrflip,
     & nrsmo1,nrsmo2,nrflop,nrsfli,nrsflx,nrinit,nrpken,
#endif
#ifdef LFM
     & nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo,klenp,weight,filtwin,
#endif
#ifdef COUPLE
     5 kpfix,ksfcx,ksig,ksfc,kpost,krestart,klfm,krsm,kcouple)
#else
     5 kpfix,ksfcx,ksig,ksfc,kpost,krestart,klfm,krsm)
#endif
c
#ifdef LFM
#ifndef NOPRINT
      if(iope)
     1 write(6,*) 'nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo,klenp,wght,fltwin=',
     1        nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo,klenp,(weight(i),i=1,10),
     2        filtwin
#endif
#endif
c
c  fhour is forecast hour of the start of the segment
c  thour is forecast hour rounded to nearest full hour
c  shour is forecast time in seconds measured from the start of the segment
c
c  shour is initialized to zero in step1
c
c  thour=fhour+int(shour/3600.+0.5e0)
c
      thour=fhour
#ifdef DFI
c
c dfini:  do digital filter initialization setup
c         nummax and numsum saved and passed in comver
      if( con(3).ne.0.0 ) then
         nummax=nint(con(3)*3600./con(1)/2.)
         numsum=-nummax-1
#ifndef NOPRINT
         if(iope) print *,' do gsm digital filter initialization '
#endif
         call dfini(0,con(3),shour/3600.,solsec,lugi)
      else
        nummax=0
        numsum=-1
      endif
#endif
c
c  zero-out arrays
c
      call zerflx(0)
      call zertmm
#ifdef NFDR
      call nfdr_zer
#endif
c
#ifdef DG3
      call inddia
      call zerdia
#endif
#ifdef DGZ
      call zerznl
#endif
c
c   read the t-dt or t=0 sigma file
c
      if(iope) then
      open (unit=n1,file='sigit ',form='unformatted',err=999)
      go to 998
  999 continue
      print *,'error opening sigit in rsmini'
#ifdef MP
#ifdef RSM
      call rmpabort
#else
      call mpabort
#endif
#else
      call abort
#endif
  998 continue
      endif
#ifndef HYBRID
      call rdsig(n1,fhour,idate,gz,qm,tem,dim,zem,rm,sl,si,z00)
#else
      call rdsig(n1,fhour,idate,gz,qm,tem,dim,zem,rm,ak5,bk5,z00)
#endif

#ifdef NISLQ

      call sph_fft_driver(-1,slq_q1,rm,ntotal_)
c   remove negative values (qmin=0.)
      slq_q1=max(slq_q1,0.)
#endif


      if( iope ) then
        rewind n1
#ifndef NOPRINT
        print 9877,n1,itread,fhour
9877    format(1h ,'n1,itread,fhour after tread',2(i4,1x),f10.1)
        print *,' input t=t0 full values'
#endif
      endif
#ifndef NOPRINT
      call rmsgt(qm,dim,tem,zem,del,rm)
#endif
c
      if(fhour.eq.0.) then
        do i=1,lnts2
          q(i)=qm(i)
        enddo
        do k=1,levs_
          do i=1,lnts2
            te(i,k)=tem(i,k)
            di(i,k)=dim(i,k)
            ze(i,k)=zem(i,k)
          enddo
        enddo
#ifndef NISLQ
        do k=1,levh_
          do i=1,lnts2
            rq(i,k)=rm(i,k)
          enddo
        enddo
#endif
#ifdef NISLQ
!
! nislq moisture at n time (n=n-1)
!
        slq_q2(:,:,:)=slq_q1(:,:,:)
#endif



        nfstep=2
        deltim=con(1)/2.e0**nfstep
        limlow=2-nfstep
        inistp=1
#ifdef KEN
        isave=1
#endif
        stepone=.true.
      else
        if(iope) then
        close(n1)
        open (unit=n1,file='sigitdt ',form='unformatted',err=889)
        go to 888
  889   continue
        print *,'error opening sigitdt in gsmini'
#ifdef MP
#ifdef RSM
        call rmpabort
#else
        call mpabort
#endif
#else
        call abort
#endif
  888   continue
        endif
#ifndef HYBRID
        call rdsig(n1,fhour,idate,gz,q,te,di,ze,rq,sl,si,z00)
#else
        call rdsig(n1,fhour,idate,gz,q,te,di,ze,rq,ak5,bk5,z00)
#endif


#ifdef NISLQ
c nislq moisture at n time (wave to grid)
        call sph_fft_driver(-1,slq_q2,rq,ntotal_)

c remove negative values (qmin=0.)
        slq_q2=max(slq_q2,0.)
#endif



        if( iope ) then
          rewind n1
#ifndef NOPRINT
          print 9878,n1,itread,fhour
9878      format(1h ,'n1,itread,fhour after tread',2(i4,1x),f10.1)
          print *,' input t=t0+dt full values'
#endif
        endif
        deltim=con(1)
        limlow=1
        inistp=0
#ifdef KEN
        isave=0
#endif
        stepone=.false.
#ifdef CLM_CWF
#ifndef SCM
        print *, 'RDWRCGS read  thour=',thour,' fhour =',fhour
        call rdwrcgs(lugi+10,fhour,idate,cgs,cwgt,inx,jnx,npnt,0)
        print *, fhour, idate
        lwconv = .false.
#endif
#endif
      endif
c
#ifndef NOPRINT
      call rmsgt(q,di,te,ze,del,rq)
#endif
c
      shour=0.0
      dtpost=0.0
      maxstp=num(7)
      dthr = con(1)/3600.e0
      hdthr = 0.5 * dthr
c
c  set initial solhr
c
      solhr=fhour+idate(1)
      iday=solhr/24.e0
      solhr=solhr-iday*24.e0
      solsec=solhr*3600.
      if( iope ) then
        print *,' initial solhr = ',solhr
      endif
c
c   read fixed fields from fixfld prog
c
      fno='sfci '
      call fixio(n1,fno,sfcftyp,
     1           labs,idate(4),idate(2),idate(3),idate(1),fhour,
     2           tsea,LONF2S,LATG2S,0)
#ifdef DBG
      print *,'just after fixio'
      call maxmin(tsea,lons2_,lats2_,1,1,'tsea')
      call maxmin(albedo,lons2_,lats2_,1,1,'albedo')
      call maxmin(slmsk,lons2_,lats2_,1,1,'slmsk')
#endif
#ifdef ISOTOPE
      if (iope) then
        close(95)
        open(95,file='rresvi',form='unformatted',status='old',err=99)
        print*,'rresvi opened'
        read(95,err=99) rresv
        close(95)
        goto 98
 99     print*,'no rresvi'
 98     continue
      endif
#ifdef MP      
      call mpgf2p(rresv(1,1,1,1),lonf2_,latg2_,
     $     rresv1,LONF2S,LATG2S,lsoil_)
      call mpgf2p(rresv(1,1,1,2),lonf2_,latg2_,
     $     rresv2,LONF2S,LATG2S,lsoil_)
#else
      do i=1,LONF2S
        do j=1,LATG2S
          do k=1,lsoil_
             rresv1(i,j,k)=rresv(i,j,k,1)
             rresv2(i,j,k)=rresv(i,j,k,2)
          enddo
        enddo
      enddo
#endif
#endif
c
#ifdef DFI
c      limlow=limlow-nummax
#endif
c
#ifdef RSM
c
c rsm : initial call to regional spectral model
c
      call rsmini(fhour,filta,ncpus,n1,n2,
     &            nrsmi1,nrsmi2,nrflip,
     &            nrsmo1,nrsmo2,nrflop,nrsfli,nrsflx,nrinit,nrpken)
#endif
c
#ifdef LFM
      call lfmfini(ipstep,fhour)
      icstep=nint(fhour*3600./deltim)+1
#ifndef NOPRINT
      if(iope) then
        write(6,*) 'ipstep from lfmfini =',ipstep
        write(6,*) 'icstep computed from fhour=',icstep
      endif
#endif
#endif
c
      if( .not. stepone ) then
        avprs0=0.0
        call gloopp(q,avprs0)
      endif
#ifndef SCM
      call gsicdf(deltim,am,bm,gv,sv,cm)
#ifdef DBG
      if(iope) print *,' gsmini:  gsicdf for deltim= ',deltim
#endif
#endif
c
#ifdef KEN
      if(npoint.gt.0) then
        isave = 1
        itnum = 1
        if(fhour.eq.0.) then
          if(istep.eq.1) then
           if (ikfreq.gt.1) isave = 0
           if (ikfreq.eq.1) itnum = 2
          end if
        end if
      end if
#endif
c
      return
      end
