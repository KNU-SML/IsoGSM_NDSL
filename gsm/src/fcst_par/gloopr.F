#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP1 jcap1
#define JCAP1P jcap1p
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
#define LCAP22P lcap22p
#define LCAPP lcapp
#define LEVH levh
#define LEVHP levhp
#define LEVM1 levm1
#define LEVP1 levp1
#define LEVS levs
#define LEVSP levsp
#define LLN22P lln22p
#define LLN2P lln2p
#define LLNP llnp
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LNUV lnuv
#define LONF22 lonf22
#define LONF22P lonf22p
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MSUB msub
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NCPUS ncpus
#define NPES npes
#define NROW nrow
#define NSOIL nsoil
#define NTOTAL ntotal
#define NTRAC ntrac
#define TWOJ1 twoj1
#define PI pi


      subroutine gloopr
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    gloop  computes dynamic non-linear tendency terms
c of temp. div. ln(ps)
c computes predicted values of vorticity and moisture
c   prgmmr: joseph sela      org: w/nmc23    date: 88-05-13
c
c abstract:
c   program  starts with spectral  coefficients temp.
c   of vorticity, divergence, specific humidity, and
c   ln((psfc).  converts them to the gaussian grid at each
c   latitude and calls fidi,  for the northern and southern
c   hemispheres at the same time.  after return from fidi
c   sr.  completes calculation of tendencies of temp. div. and lnps.
c   specific humidity, and vorticity are predicted by sr. sigvor
c   all input/output  is via commons.
c
c program history log:
c   91-03-06  joseph sela
c   99-01-04  hann-ming henry juang	merge loopb and loopr in
c   99-05-14  hann-ming henry juang	add mpi
c
c usage:    call gloopa
c
c attributes:
c   language: fortran 77
c
c$$$
c....
#include <paramodel.h>
#include <comfibm.h>
#include <comgrad.h>
#include <comgda.h>
#ifdef NISLQ
#include <comnislq.h>
#endif
#ifndef DYNAMIC_ALLOC
      parameter(ncpus=ncpus_)
#endif
c...sof include..........................................
#ifdef MP
#define JCAP1S jcap1p_
#define LONF2S lonf2p_
#define LEVSS levsp_
#define LEVHS levhp_
#define LCAP22S lcap22p_
#define LONF22S lonf22p_
#else
#define JCAP1S jcap1_
#define LONF2S lonf2_
#define LEVSS levs_
#define LEVHS levh_
#define LCAP22S lonf22_
#define LONF22S lonf22_
#endif




c...
c.................................................................
c syn(1, 0*levs_+1, lan)  di
c syn(1, 1*levs_+1, lan)  te
c syn(1, 2*levs_+1, lan)  rq
c syn(1, 3*levs_+1, lan)  uln
c syn(1, 4*levs_+1, lan)  vln
c syn(1, 5*levs_+1, lan)  dpdphi
c syn(1, 5*levs_+2, lan)  dpdlam
c syn(1, 5*levs_+3, lan)  q
c.................................................................
#ifndef NISLQ
      parameter(lots    = 4*levs_+levh_+3,
#else
      parameter(lots    = 4*levs_+3,
#endif
     &          lotst   =2*levs_+1,
     &          ksd     = 1,
     &          kst     = ksd+levs_,
#ifndef NISLQ
     &          ksr     = kst+levs_,
     &          ksu     = ksr+levh_, kstb=ksu,
#else
     &          ksu     = kst+levs_, kstb=ksu,
#endif
     &          ksv     = ksu+levs_,
     &          kspphi  = ksv+levs_,
     &          ksplam  = kspphi+1,
     &          ksp     = ksplam+1)
c...........................................................
#ifndef NISLQ
      parameter(lotss   = 4*LEVSS+LEVHS+3, 
#else
      parameter(lotss   = 4*LEVSS+3, 
#endif
     &          lotsts=2*LEVSS+1,
     &          ksds    = 1,
     &          ksts    = ksds+LEVSS,
#ifndef NISLQ
     &          ksrs    = ksts+LEVSS,
     &          ksus    = ksrs+LEVHS, kstbs=ksus,
#else
     &          ksus    = ksts+LEVSS, kstbs=ksus,
#endif
     &          ksvs    = ksus+LEVSS,
     &          kspphis = ksvs+LEVSS,
     &          ksplams = kspphis+1,
     &          ksps    = ksplams+1)
c...........................................................
c local array
#ifdef SCM
      common /comscmbdy/dtbdy, curtime, vvel(lnt22_,levs_)
#ifdef CLM_CWF
     +                , wdiv(lnt22_,levs_), hadq(lnt22_,levs_)
#endif
#endif
#ifdef MP
#include <commpi.h>
#include <compspec.h>
c
      dimension
     1 syf(lonf22_,lotss,latg2p_),grs(lonf22p_,lots,latg2p_)
#endif
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
c ...............................................................
#ifdef MP
#define NCPUSS latg2_
#else
#define NCPUSS ncpus
#endif
      dimension scs2(JCAP1S),
     2 syn(LCAP22S,lotss,NCPUSS),syntop(2,JCAP1S,lotsts)
#undef NCPUSS
c ...............................................................
#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus
#endif
#ifdef DG
      dimension cldt(LONF2S,levs_,NCPUSS),clcv(LONF2S,levs_,NCPUSS)
#endif
#ifdef DGZ
      dimension sfcp(LONF2S,NCPUSS)
#endif
#ifdef DG3
      dimension gda(nwgda,kdgda,NCPUSS)
#endif
#undef NCPUSS

#ifdef CONDENS
c ...............................................................
      common/comcond/
     +       qcicnv(LONF2S,levs_,latg2_), qrscnv(LONF2S,levs_,latg2_),
     +       qcilrg(LONF2S,levs_,latg2_), qrslrg(LONF2S,levs_,latg2_),
     +       taucld(LONF2S,levs_,latg2_),
     +       cldwp(LONF2S,levs_,latg2_),  cldip(LONF2S,levs_,latg2_)
#endif
c ...............................................................


#ifdef NISLQ
      real, dimension(LONF22S,levh_,LATG2S) :: slq_in
#endif

c ...............................................................
      logical dorad
c
      call gradpre(idate,dorad)
c
      if(.not.dorad) return
c......................................................................
#ifdef MP
#ifndef NOPRINT
      if(mype.eq.master) 
     1  print *,' run radiation at solar hour=',solhr
#endif
c...................................................................
      llstr=lwvstr(mype)
      llens=lwvlen(mype)
      jstr=latstr(mype)
      jend=latstr(mype)+latlen(mype)-1
      lons2=lonlen(mype)*2
      lats2=latlen(mype)
#else
      llstr=0
      llens=jcap1_
      lons2=lonf2_
      lats2=latg2_
#endif
c=================================================================
#ifndef SCM
#ifdef MP

#ifndef NISLQ
      call mpnn2nk(ze ,lnt22p_,levs_, zea,lln22p_,levsp_,3+ntotal_)
#else
      call mpnn2nk(ze ,lnt22p_,levs_, zea,lln22p_,levsp_,3)
#endif
      call mpnn2n(q ,lnt22p_, qa,lln22p_,1)
#define QS qa
#define DPDPHIS dpdphia
#define DPDLAMS dpdlama
#define DIS dia
#define ZES zea
#define ULNS ulna
#define VLNS vlna
#else
#define QS q
#define DPDPHIS dpdphi
#define DPDLAMS dpdlam
#define DIS di
#define ZES ze
#define ULNS uln
#define VLNS vln
#endif
      call dellnp(QS,DPDPHIS,syntop(1,1,2*LEVSS+1),DPDLAMS,
     1            llstr,llens,lwvdef)
      call dzuvle(DIS,ZES,ULNS,VLNS,syntop(1,1,1),
     1            syntop(1,1,LEVSS+1),llstr,llens,lwvdef)
#undef QS
#undef DPDPHIS
#undef DPDLAMS
#undef DIS
#undef ZES
#undef ULNS
#undef VLNS
#endif                  /* SCM */
c...................................................................
#ifndef MP
c compute latitude band limits
      last=mod(latg2_,ncpus)
      nggs=(latg2_-last)/ncpus
      if(last.ne.0)nggs=nggs+1
      inclat=ncpus
      lat1=1-ncpus
      lat2=0
      latdon=0
cc
      do 10000 ngg=1,nggs
      if((ngg.eq.nggs).and.(last.ne.0)) inclat=last
      lat1=lat1+ncpus
      lat2=lat2+inclat
#endif
c.................................................................
#ifdef DG3
#ifdef MP
      lat1=jstr
      lat2=jend
      latdon=jstr-1
#define LATX lan
#else
#define LATX lat
#endif
      do lat=lat1,lat2
        lan=lat-latdon
        call getdia(LATX,nwgda*kdgda,gda(1,1,lan))
      enddo
#undef LATX
#endif
c...................................................................
#ifndef SCM
#ifdef MP
      lat1=1
      lat2=latg2_
      latdon=0
#define DIS dia
#else
#define DIS di
#endif
c first lat loop
#ifdef ORIGIN_THREAD
c$doacross share(syntop,syn,qtt,qvv,lat1,lat2,latdon,
c$&        colrad,dia,di,llstr,llens,lwvdef,lcapd,lcapdp,mype)
c$&        local(lat,lan,llensd)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(syntop,syn,qtt,qvv,lat1,lat2,latdon)
cmic$1 shared(colrad,dia,di,llstr,llens,lwvdef,lcapd,lcapdp,mype)
cmic$1 private(lat,lan,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,llensd)
#endif
c
      do 1000 lat=lat1,lat2
        lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
        llensd=lcapdp(lat,mype)
#else
        llensd=lcapd(lat)
#endif
#else
        llensd=llens
#endif
        call sums2i(DIS,syn(1,1,lan),qtt(1,lat),
     &            llstr,llensd,lwvdef,lotss)
        call sumtop(syn(1,kstbs,lan),syntop,qvv(1,lat),
     &            llstr,llensd,lwvdef,lotsts)
1000  continue
c
#undef ZES
c......................................................................
#ifdef MP
      call mpnl2ny(syn,lcap22p_,latg2_,
     1             syf,lonf22_,latg2p_,lotss,1,lotss)
c
      lat1=jstr
      lat2=jend
      latdon=jstr-1
#define SYNS syf
#else
#define SYNS syn
#endif
c
#ifdef ORIGIN_THREAD
c$doacross share(syn,syf,lat1,lat2,latdon,lcapd,lonfd,latdef),
c$&        local(lat,lan,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$& shared(syn,syf)
cmic$& shared(lat1,lat2,latdon,lcapd,lonfd,latdef)
cmic$& private(lat,lan,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,lcapf,lonff)
#endif
       do 1300 lat=lat1,lat2
         lan=lat-latdon
#ifdef REDUCE_GRID
         lcapf=lcapd(latdef(lat))
         lonff=lonfd(latdef(lat))
#else
         lcapf=jcap1_
         lonff=lonf_
#endif
         call ftilonf(SYNS(1,1,lan),SYNS(1,1,lan),lotss*2,
     1                lcapf,lonff,latdef(lat),1)
 1300  continue
c
#undef SYNS
#undef DYNS
#endif                            /* SCM */
c....................................................................
#ifdef MP
#ifndef NISLQ
             klev=4+ntotal_
#else
             klev=4
#endif
       call mpnk2nx(syf,lonf22_,lotss,
     &              grs,lonf22p_,lots,latg2p_,levsp_,levs_,1,1,
     &              klev)
       call mpx2nx (syf,lonf22_,lotss,
     &              grs,lonf22p_,lots,latg2p_,kspphis,kspphi,3)
c
      lat1=jstr
      lat2=jend
      latdon=jstr-1
#define SYNS grs
#define LATX lan
#else
#define SYNS syn
#define LATX lat
#endif
c ............................................................
#ifdef ORIGIN_THREAD
c$doacross share(syn,grs,lat1,lat2,latdon,gda,latdef,
c$&              cldt,clcv,sfcp,lons2,lats2,lonfdp,mype),
c$&        local(lat,lan,lonsd2,latrue)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$& shared(syn,grs,lat1,lat2,latdon,gda,latdef)
cmic$& shared(cldt,clcv,sfcp,lons2,lats2,lonfdp,mype)
cmic$& private(lat,lan,lonsd2,latrue)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,lonsd2,latrue)
#endif
c
        do 1400 lat=lat1,lat2
#ifdef SCM
          nxlat = nint(xlat(1,lat)*180./pi_)
          latrue = 90 - nxlat
          if (latrue.gt.90) latrue = 180 - latrue
#else
          latrue=latdef(lat)
#endif
          lan=lat-latdon
#ifdef REDUCE_GRID
#ifdef MP
          lonsd2=lonfdp(lan,mype)*2
#else
          lonsd2=lonfd(latdef(lat))*2
#endif
#else
          lonsd2=lons2
#endif

#ifdef NISLQ
!
!      nislq moisture at n+1 time
!
       do k = 1,levh_
         do i=1,lonsd2
c           if ( stepone ) then
             slq_in(i,k,lan) = slq_q2(i,k,LATX)
c           else
c             slq_in(i,k,lan) = slq_q3(i,k,LATX)
c           endif
         enddo
       enddo
#endif /* NISLQ end */



#ifdef MP
          if( lonsd2.gt.0 ) then
#endif
          call grphys(lonsd2,lats2,
#ifdef SCM
     &                di,ze,q,te,rq,vvel,
#else
     &    SYNS(1,ksplam,lan),SYNS(1,kspphi,lan),
     &    SYNS(1,ksu,lan),SYNS(1,ksv,lan),SYNS(1,ksp,lan),
#ifdef NISLQ
     &    SYNS(1,kst,lan),slq_in(1,1,lan),
#else
     &    SYNS(1,kst,lan),SYNS(1,ksr,lan),
#endif
     &    SYNS(1,ksd,lan),
#endif
#ifdef DG3
     &    gda(1,1,lan),
#endif
#ifdef DG
     &    cldt(1,1,lan),clcv(1,1,lan),
#endif
#ifdef DGZ
     &    sfcp(1,lan),
#endif
#ifdef CONDENS
     +    taucld(1,1,LATX),
     +    cldwp(1,1,LATX), cldip(1,1,LATX),
#endif
     &    LATX,latrue)
#ifdef MP
          endif
#endif
 1400   continue
c
c ................................................................
#ifdef DG3
#ifdef MP
#define LATX lan
#else
#define LATX lat
#endif
        do lat=lat1,lat2
          lan=lat-latdon
          call putdia(LATX,nwgda*kdgda,gda(1,1,lan))
        enddo
#undef LATX
#endif
c ............................................................
#ifndef MP
      latdon=latdon+(lat2-lat1+1)
10000 continue
c..................................................................
#endif
c
      if(inistp.ne.0) then
        runrad=.false.
      else
        runrad=.true.
      endif
c
      return
      end
