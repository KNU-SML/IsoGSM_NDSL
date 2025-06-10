#define LONF22S lonf22s
#define CP cp
#define G g
#define HFUS hfus
#define HVAP hvap
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
#define RD rd
#define TWOJ1 twoj1
      subroutine gloopb
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    gloop  computes model physics tendency terms
c
c   prgmmr: joseph sela      org: w/nmc23    date: 88-05-13
c
c abstract:
c   program  starts with spectral  coefficients temp.
c   of vorticity, divergence, specific humidity, and
c   ln((psfc).  converts them to the gaussian grid at each
c   latitude and calls gbphys,  for the northern and southern
c   hemispheres at the same time.  after return from fidi
c   sr.  completes calculation of tendencies of temp. div. and lnps.
c   specific humidity, and vorticity are predicted by sr. sigvor
c   all input/output  is via commons.
c
c program history log:
c   91-03-06  joseph sela
c   99-05-14  hann-ming henry juang	add mpi one loop
c   99-06-24  hann-ming henry juang	add mpi full loop
c
c usage:    call gloopb
c
c attributes:
c   language: fortran 77
c
c$$$
c
#include <paramodel.h>
#ifdef RASV2
#undef RAS
#endif
#include <comfibm.h>
#ifdef DG
#include <comgda.h>
#endif
#ifdef NISLQ
#include <comnislq.h>
#endif
#ifndef DYNAMIC_ALLOC
      parameter(ncpus=ncpus_)
#endif
c
#ifdef MP
#define JCAP1S jcap1p_
#define LEVSS levsp_
#define LEVHS levhp_
#define LCAPS lcapp_
#define LCAP22S lcap22p_
#else
#define JCAP1S jcap1_
#define LEVSS levs_
#define LEVHS levh_
#define LCAPS lonf_
#define LCAP22S lonf22_
#endif
c
c syn(1, 0*levs_+1, lan)  y
c syn(1, 1*levs_+1, lan)  rt
c syn(1, 2*levs_+1, lan)  x
c syn(1, 3*levs_+1, lan)  uln
c syn(1, 4*levs_+1, lan)  vln
c syn(1, 5*levs_+1, lan)  dpdphi
c syn(1, 5*levs_+2, lan)  dpdlam
c syn(1, 5*levs_+3, lan)  q
c syn(1, 5*levs_+4, lan)  qlap
c
c anl(1, 0*levs_+1, lan)  uu    dvdt
c anl(1, 1*levs_+1, lan)  vv    dtdt
c anl(1, 2*levs_+1, lan)  te    dtdt
c anl(1, 3*levs_+1, lan)  rq    drdt
c

#ifdef NISLQ
      parameter(lots    =4*levs_+4,
#else
      parameter(lots    =4*levs_+levh_+4,
#endif
     &          lotst   =2*levs_+1,
     &          kst     =1,
#ifndef NISLQ
     &          ksr     =kst+levs_,
     &          ksd     =ksr+levh_,
#else
     &          ksd     =kst+levs_,
#endif
     &          ksu     =ksd+levs_,
     &          kstb    =ksd+levs_,
     &          ksv     =ksu+levs_,
     &          kspphi  =ksv+levs_,
     &          ksplam  =kspphi+1,
     &          ksp     =ksplam+1,
     &          ksplap  =ksp   +1)
#ifndef NISLQ
      parameter(lota    =3*levs_+levh_,
#else
      parameter(lota    =3*levs_,
#endif
     &          lotat   =2*levs_,
     &          kau     =1,
     &          kav     =kau+levs_,
     &          kat     =kav+levs_
#ifndef NSILQ
     &,         kar     =kat+levs_)
#else
     &          )
#endif

#ifndef NISLQ
      parameter(lotss   =4*LEVSS+LEVHS+4,
#else
      parameter(lotss   =4*LEVSS+4,
#endif
     &          lotsts  =2*LEVSS+1,
     &          ksts    =1,
#ifndef NISLQ
     &          ksrs    =ksts+LEVSS,
     &          ksds    =ksrs+LEVHS,
#else
     &          ksds    =ksts+LEVSS,
#endif
     &          ksus    =ksds+LEVSS,
     &          kstbs   =ksds+LEVSS,
     &          ksvs    =ksus+LEVSS,
     &          kspphis =ksvs+LEVSS,
     &          ksplams =kspphis+1,
     &          ksps    =ksplams+1,
     &          ksplaps =ksps   +1)
#ifndef NISLQ
      parameter(lotas   =3*LEVSS+LEVHS,
#else
      parameter(lotas   =3*LEVSS,
#endif
     &          lotats  =2*LEVSS,
     &          kaus    =1,
     &          kavs    =kaus+LEVSS,
     &          kats    =kavs+LEVSS
#ifndef NISLQ
     &,         kars    =kats+LEVSS)
#else
     &          )
#endif

c
      parameter(no3p=28,no3l=29)


#ifdef SCM
      common /comscmbdy/dtbdy, curtime, vvel(lnt22_,levs_)
#ifdef CLM_CWF
     +                , wdiv(lnt22_,levs_), hadq(lnt22_,levs_)
#endif
#endif                  /* SCM */
#ifdef MP
#include <commpi.h>
#include <compspec.h>
c
      dimension
     1 syf(lonf22_,lotss,latg2p_),grs(lonf22p_,lots,latg2p_),
     3 anf(lonf22_,lotas,latg2p_),gra(lonf22p_,lota,latg2p_)
#endif
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
c
#ifdef MP
#define NCPUSS latg2_
#else
#define NCPUSS ncpus
#endif
      dimension
     2 syn(LCAP22S,lotss,NCPUSS),syntop(2,JCAP1S,lotsts),
     4 anl(LCAP22S,lotas,NCPUSS),anltop(2,JCAP1S,lotats),
     5 flp(2,JCAP1S,lotas,NCPUSS),flm(2,JCAP1S,lotas,NCPUSS)
#undef NCPUSS
c
#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus
#endif
#ifdef DG
      dimension tgmxl(NCPUSS),igmxl(NCPUSS),kgmxl(NCPUSS)
      dimension tgmnl(NCPUSS),igmnl(NCPUSS),kgmnl(NCPUSS)
#endif
#ifdef DG3
      dimension gda(nwgda,kdgda,NCPUSS)
#endif
#undef NCPUSS
      logical ladj
#ifdef SCM
      parameter(ladj=.false.)
#else
      parameter(ladj=.true.)
#endif
#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif

c#ifdef PSPLIT
#ifdef MP
#define LONF22S lonf22p_
#else
#define LONF22S lonf22_
#endif
#include <compsplit.h>
c#endif
c 

c ...............................................................
#ifdef CLM_CWF
      parameter (npmx=1000)
      common /comfcwrk/cgs(LONF2S,levs_,latg2_),
     +                 cwgt(npmx,lonf2_,latg2_)
      common /comfcwrki/inx(npmx,lonf2_,latg2_),jnx(npmx,lonf2_,latg2_),
     +                 npnt(lonf2_,latg2_), lwconv
      logical lwconv
      dimension qugs(LONF2S*latg2_,levs_), qvgs(LONF2S*latg2_,levs_)
      dimension qu(lonf_,latg_), qv(lonf_,latg_), conv(lonf2_,latg2_),
     +          wconv(lonf2_,latg2_)
      character*80 fno
c
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
#ifdef CONDENS
c ...............................................................
      common/comcond/
     +       qcicnv(LONF2S,levs_,latg2_), qrscnv(LONF2S,levs_,latg2_),
     +       qcilrg(LONF2S,levs_,latg2_), qrslrg(LONF2S,levs_,latg2_),
     +       taucld(LONF2S,levs_,latg2_),
     +       cldwp(LONF2S,levs_,latg2_),  cldip(LONF2S,levs_,latg2_)
#endif

#ifdef NISLQ
        real, dimension(LONF22S,levh_,LATG2S) :: slq_in, slq_out
#endif

#ifdef MUL_CLDTOP
c ...............................................................
c
      dimension xkt2(lonf2_*ncldtop_,latg2_)
      Integer, allocatable ::  NRND(:)
      integer krsize, iseed
      logical first
      real    seed0
      data first/.true./
      save    first, krsize, nrnd, seed0



c
      if (first) then
        CALL RANDOM_SEED(SIZE=krsize)
        allocate (nrnd(krsize))
        seed0 = idate(1) + idate(2) + idate(3) + idate(4)
cyh        write(99,*) ' fisrt=',first,' seed0=',seed0
        first = .false.
      endif
      me = 0.0
      iseed = mod(100.0*sqrt(fhour*3600+100.0*me),1.0E9) + 1 + seed0
      nrnd  = iseed
#ifdef IBM
      call random_seed(generator=2)
#endif
      CALL RANDOM_SEED(PUT=NRND)
      CALL RANDOM_NUMBER(xkt2)
#ifdef SCM
      xkt2m = 0.0
      do jj = 1, latg2_
      do ii = 1, lonf2_*ncldtop_
        xkt2m = xkt2m + xkt2(ii,jj)
      enddo
      enddo
      xkt2m = xkt2m/real(latg2_*lonf2_*ncldtop_)
      do jj = 1, latg2_
      do ii = 1, lonf2_*ncldtop_
        xkt2(ii,jj) = xkt2m
      enddo
      enddo
#endif                  /* SCM */
#endif                  /* MUL_CLDTOP */
#
#ifdef RAS
      logical ras
      parameter (ras=.true.)
      parameter (nsphys=1)
      parameter (cp=cp_, alhl=hvap_, grav=g_, rgas=rd_)
      dimension sig(levs_+1), prj(levs_+1), prh(levs_)
     *,         fpk(levs_),   hpk(levs_),   sgb(levs_)
     *,         ods(levs_), rasal(levs_),  prns(levs_/2)
     *,      rannum(200)
c
       call setras(levs_, si, sl, del, cp, rgas, deltim
     *,            nsphys, thour
     *,            sig, sgb, prh, prj, hpk, fpk, ods, prns
     *,            rasal, lm, krmin, krmax, nstrp
     *,            ncrnd, rannum, afac, ufac)
#endif
#ifdef RASV2
      logical ras
      parameter (ras=.true.)
      parameter (nsphys=1)
      parameter (cp=cp_, alhl=hvap_, grav=g_, rgas=rd_, hfus=hfus_)
      real sig(levs_+1), prj(levs_+1), sgb(levs_),   rasal(levs_)
     *,    rannum(200,5), pdd, dsfc(LONF2S,LATG2S)
      integer mct
c
c  pdd is the lowest pressure above which downdraft is allowed
c  mct is a number of cloud types
c
      parameter (pdd=600.0, mct=15)
      integer krmin, krmax, kfmax, ncrnd, kctop(mct+1)
      common /rassav/sig, sgb, prj, rasal, krmin, krmax, kfmax
     *,              kctop, ncrnd, dsfc
c
        frain = .5
        if(inistp.eq.1) frain = 1.
        dtp= deltim / frain
c
c  check fhour and thour
c
        call setrasv2(levs_,  si,  sl, cp_, rd_, dtp, nsphys, fhour
     *,             sig,   sgb, prj
     *,             rasal, krmin, krmax, kfmax, ncrnd, rannum
     *,             mct,   kctop, deltim, dsfc, LONF2S*LATG2S)
#endif
c
#ifndef O3C
       if(ntrac_.ge.2) call seto3p(idate,fhour,no3p,no3l)
#endif

#ifdef NISLQ
       slq_in=0.
       slq_out=0.
#endif


 
c
#ifdef MP
      llstr=lwvstr(mype)
      llens=lwvlen(mype)
      jstr=latstr(mype)
      jend=latstr(mype)+latlen(mype)-1
      lons2=lonlen(mype)*2
      lats2=latlen(mype)
      lnts2=lntlen(mype)*2
      lnoffset=lntstr(mype)*2
#else
      llstr=0
      llens=jcap1_
      lons2=lonf2_
      lats2=latg2_
      lnts2=lnt2_
      lnoffset=0
#endif
#ifdef PSPLIT
c
      do k=1,levs_
        do i=1,LNT22S
          y (i,k)=tem(i,k)
          x (i,k)=dim(i,k)
          w (i,k)=zem(i,k)
        enddo
      enddo
      do k=1,levh_
        do i=1,LNT22S
          rt(i,k)= rm(i,k)
        end do
      end do
      do i=1,LNT22S
        q(i)=qm(i)
      enddo

#endif
c
#ifndef SCM
#ifdef MP
      call mpnn2n(q ,lnt22p_, qa,lln22p_,1)
#ifndef NISLQ
      call mpnn2nk(y ,lnt22p_,levs_, ya,lln22p_,levsp_,3+ntotal_)
#else
      call mpnn2nk(y ,lnt22p_,levs_, ya,lln22p_,levsp_,3)
#endif
#define QS qa
#define DPDPHIS dpdphia
#define DPDLAMS dpdlama
#define QLAPS qlapa
#define XS xa
#define WS wa
#define ULNS ulna
#define VLNS vlna
#else
#define QS q
#define DPDPHIS dpdphi
#define DPDLAMS dpdlam
#define QLAPS qlap
#define XS x
#define WS w
#define ULNS uln
#define VLNS vln
#endif
      call dellnp(QS,DPDPHIS,syntop(1,1,2*LEVSS+1),DPDLAMS,
     1            llstr,llens,lwvdef)
      call deldfsp(QS,QLAPS,llstr,llens,lwvdef)
      call dzuvle(XS,WS,ULNS,VLNS,syntop(1,1,1),
     1            syntop(1,1,LEVSS+1),llstr,llens,lwvdef)
#undef QS
#undef DPDPHIS
#undef DPDLAMS
#undef QLAPS
#undef XS
#undef WS
#undef ULNS
#undef VLNS
c
#ifndef PSPLIT
      do n=1,2
        do j=1,JCAP1S
          do l=1,lotats
            anltop(n,j,l)=0.0
          enddo
        enddo
      enddo
#endif
#ifdef DG
      tgmx=-1.e20
      tgmn= 1.e20
#endif
#endif                                  /* SCM */
c
#ifdef MP
#define LNT2X lln22p_
#define ZES zea
#define DIS dia
#define TES tea
#define RQS rqa
#else
#define LNT2X lnt2_
#define ZES ze
#define DIS di
#define TES te
#define RQS rq
#endif
c
#ifndef PSPLIT
#ifdef ORIGIN_THREAD
c$doacross share(ze,zea,di,dia,te,tea,rq,rqa),local(j,k)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(ze,zea,di,dia,te,tea,rq,rqa)
cmic$1 private(j,k)
#endif
#ifdef OPENMP
!$omp parallel do private(j,k)
#endif
      do k=1,LEVSS
       do j=1,LNT2X
        ZES(j,k)=0.0
        DIS(j,k)=0.0
        TES(j,k)=0.0
       enddo
      enddo
#ifndef NISLQ
      do k=1,LEVHS
       do j=1,LNT2X
        RQS(j,k)=0.0
       enddo
      enddo
#endif
#undef LNT2X
#undef ZES
#undef DIS
#undef TES
#undef RQS
#endif
c
#ifndef MP
c
c compute latitude band limits
c
      last=mod(latg2_,ncpus)
      nggs=(latg2_-last)/ncpus
      if(last.ne.0)nggs=nggs+1
      inclat=ncpus
      lat1=1-ncpus
      lat2=0
      latdon=0
c
      do 10000 ngg=1,nggs
      if((ngg.eq.nggs).and.(last.ne.0)) inclat=last
      lat1=lat1+ncpus
      lat2=lat2+inclat
#endif
c
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
#ifdef DG
#ifdef MP
#define NCPUSS latg2p_
#else
#define NCPUSS ncpus
#endif
      do k=1,NCPUSS
        tgmxl(k)=tgmx
        tgmnl(k)=tgmn
      enddo
#undef NCPUSS
#endif
c
#ifndef SCM
#ifdef MP
      lat1=1
      lat2=latg2_
      latdon=0
#define YS ya
#define ULNS ulna
#else
#define YS y
#define ULNS uln
#endif
c first lat loop
#ifdef ORIGIN_THREAD
c$doacross share(syntop,syn,qtt,qvv,lat1,lat2,latdon,
c$&        colrad,ya,y,ulna,uln,llstr,llens,lwvdef,lcapd,lcapdp,mype)
c$&        local(lat,lan,llensd)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(syntop,syn,qtt,qvv,lat1,lat2,latdon)
cmic$1 shared(colrad,ya,y,ulna,uln,llstr,llens,lwvdef,lcapd,lcapdp,mype)
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

#ifndef NISLQ
        call sums2i(YS,syn(1,ksts,lan),qtt(1,lat),
     &            llstr,llensd,lwvdef,LEVSS*2+LEVHS)
        call sums2i(ULNS,syn(1,ksus,lan),qtt(1,lat),
     &            llstr,llensd,lwvdef,LEVSS*2+4)
#else
        call sums2i(YS,syn(1,ksts,lan),qtt(1,lat),
     &            llstr,llensd,lwvdef,LEVSS*2)
        call sums2i(ULNS,syn(1,ksus,lan),qtt(1,lat),
     &            llstr,llensd,lwvdef,LEVSS*2+4)
#endif
        call sumtop(syn(1,kstbs,lan),syntop,qvv(1,lat),
     &            llstr,llensd,lwvdef,lotsts)
1000  continue
#undef YS
#undef ULNS
c
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
c$doacross share(syf,syn,lat1,lat2,latdon,latdef,lcapd,lonfd),
c$&        local(lat,lan,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$& shared(syf,syn,latdef,lcapd,lonfd)
cmic$& shared(lat1,lat2,latdon)
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
#endif                          /* SCM */
c
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
     &              grs,lonf22p_,lots,latg2p_,kspphis,kspphi,4)
c
      lat1=jstr
      lat2=jend
      latdon=jstr-1
#define SYNS grs
#define ANLS gra
#define LATX lan
#else
#define SYNS syn
#define ANLS anl
#define LATX lat
#endif
c
#ifdef ORIGIN_THREAD
c$doacross share(syn,grs,anl,gra,lat1,lat2,latdon,gda,rdt,lons2,
c$&              tmgx,igmx,kgmx,tgmn,igmn,kgmn,
c$&              tmgxl,igmxl,kgmxl,tgmnl,igmnl,kgmnl,ras,
c$&              cp,alhl,grav,rgas,sig,sgb,prh,prj,hpk,fpk,ods,
c$&              prns,rasal,lm,krmin,krmax,nstrp,ncrnd,rannum,afac,
c$&              ufac,lonfdp,mype),local(lat,lan,lonsd2,jgmx,jgmn)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$& shared(syn,grs,anl,gra,lat1,lat2,latdon,gda,rdt,lons2)
cmic$& shared(tgmx ,igmx ,kgmx ,tgmn ,igmn ,kgmn )
cmic$& shared(tgmxl,igmxl,kgmxl,tgmnl,igmnl,kgmnl)
cmic$& shared(ras,cp,alhl,grav,rgas)
cmic$& shared(sig, sgb, prh, prj, hpk, fpk, ods, prns)
cmic$& shared(rasal, lm, krmin, krmax, nstrp)
cmic$& shared(ncrnd, rannum, afac, ufac,lonfdp,mype,jgmx,jgmn)
cmic$& private(lat,lan,lonsd2)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,lonsd2)
#endif




c
      do 1500 lat=lat1,lat2
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
           slq_in(i,k,lan) = slq_q3(i,k,LATX)
         enddo
       enddo
#endif /* NISLQ end */


c
#ifdef MP
        if( lonsd2.gt.0 ) then
#endif
c
        call gbphys(lonsd2,
#ifdef SCM
     +  x,w,z,y,rt,vvel,
     +  te,rq,di,ze,
#else
     x  SYNS(1,ksplam,lan),SYNS(1,kspphi,lan),
     x  SYNS(1,ksu,lan),SYNS(1,ksv,lan),SYNS(1,ksp,lan),
#ifdef NISLQ
     x  SYNS(1,kst,lan),slq_in(1,1,lan),
#else
     x  SYNS(1,kst,lan),SYNS(1,ksr,lan),
#endif
     x  SYNS(1,ksd,lan),
     x  SYNS(1,ksplap,lan),
#ifdef PSPLIT
     x  dtdtm(1,1,LATX),drdtm(1,1,LATX),
     x  dudtm(1,1,LATX),dvdtm(1,1,LATX),
#else
#ifdef NISLQ
     x  ANLS(1,kat,lan),slq_out(1,1,lan),
#else
     x  ANLS(1,kat,lan),ANLS(1,kar,lan),
#endif
     x  ANLS(1,kau,lan),ANLS(1,kav,lan),
#endif
#endif                       /* SCM */
#ifdef DG
     x  tgmxl(lan),igmxl(lan),kgmxl(lan),
     x  tgmnl(lan),igmnl(lan),kgmnl(lan),
#endif
#ifdef DG3
     x  gda(1,1,lan),
#endif
#ifdef RAS
     x  ras,levs_,cp,alhl,grav,rgas,
     x  sig, sgb, prh, prj, hpk, fpk, ods, prns,
     x  rasal, lm, krmin, krmax, nstrp,
     x  ncrnd, rannum, afac, ufac,
#endif
#ifdef RASV2
     &  ras,rgas, cp, grav, alhl,
     &  sig,prj,sgb,rasal, rannum, dsfc(1,LATX),
     &  pdd,krmin, krmax, kfmax, ncrnd, mct,kctop,
#endif

#ifdef MUL_CLDTOP
     x  xkt2(1,lan),ncldtop_,
#endif
#ifdef CLM_CWF
     x  qugs, qvgs, cgs(1,1,LATX), lwconv,
#endif
#ifdef SAS_DIAG
     +  dcu(1,1,LATX), dcv(1,1,LATX), dct(1,1,LATX),
     +  dcq(1,1,LATX), dch(1,1,LATX),
     +  fcu(1,1,LATX), fcd(1,1,LATX),
     +  deltb(1,LATX), delqb(1,LATX), delhb(1,LATX), cbmf(1,LATX),
     +  dlt(1,1,LATX), dlq(1,1,LATX), dlh(1,1,LATX),
#endif
#ifdef CONDENS
     +  qcicnv(1,1,LATX), qrscnv(1,1,LATX),
     +  qcilrg(1,1,LATX), qrslrg(1,1,LATX),
#endif


#ifdef RSM
     x  LATX,1.0)
#else
#ifndef VIC
     x  LATX,1.0)
#else
     x  LATX,1.0,idate)
#endif
#endif
c
#ifdef DG
        if(tgmxl(lan).gt.tgmx) then
          tgmx=tgmxl(lan)
          igmx=igmxl(lan)
          kgmx=kgmxl(lan)
          jgmx=lat
        endif
        if(tgmnl(lan).lt.tgmn) then
          tgmn=tgmnl(lan)
          igmn=igmnl(lan)
          kgmn=kgmnl(lan)
          jgmn=lat
        endif
#endif
#ifdef MP
        endif
#endif

#ifdef NISLQ
!
!      nislq moisture update at n+1
!
       do k = 1,levh_
         do i=1,lonsd2
           slq_q3(i,k,LATX)=slq_out(i,k,lan)
         enddo
       enddo
#endif

#undef LATX
c
 1500 continue
#ifndef PSPLIT
c
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
c
#undef SYNS
#undef DYNS
#undef ANLS
c
#ifndef SCM
#ifdef MP
#ifndef NISLQ
      klev=3+ntotal_
#else
      klev=3
#endif
      call mpnx2nk(gra,lonf22p_,lota,
     &             anf,lonf22_ ,lotas,latg2p_,levs_,levsp_,kau,kaus,
     &             klev)
      if( levlen(mype).lt.levsp_ ) then
        do lan=1,lats2
         do i=1,lonf22_
          anf(i,kaus+levsp_-1,lan)=0.0
          anf(i,kavs+levsp_-1,lan)=0.0
          anf(i,kats+levsp_-1,lan)=0.0
         enddo
        enddo
#ifndef NISLQ
        do lan=1,lats2
         do i=1,lonf22_
           do kc = 1,ntotal_
             anf(i,kars+levsp_*kc-1,lan)=0.0
           enddo
         enddo
        enddo
#endif
      endif
c
      lat1=jstr
      lat2=jend
      latdon=jstr-1
#define ANLS anf
#else
#define ANLS anl
#endif
c
#ifdef ORIGIN_THREAD
c$doacross share(latdon,lat1,lat2,anl,anf,latdef,lcapd,lonfd),
c$&        local(lat,lan,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(latdon,lat1,lat2,anl,anf,latdef,lcapd,lonfd)
cmic$1 private(lat,lan,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,lcapf,lonff)
#endif
c
      do 2200 lat=lat1,lat2
        lan=lat-latdon
#ifdef REDUCE_GRID
        lcapf=lcapd(latdef(lat))
        lonff=lonfd(latdef(lat))
#else
        lcapf=jcap1_
        lonff=lonf_
#endif
        call ftilonf(ANLS(1,1,lan),ANLS(1,1,lan),2*lotas,
     1               lcapf,lonff,latdef(lat),-1)
2200  continue
#undef ANLS
c
#ifdef MP
      call mpny2nl(anf,lonf22_ ,latg2p_,
     1             anl,lcap22p_,latg2_ ,lotas,kaus,lotas)
c
      lat1=1
      lat2=latg2_
      latdon=0
#define ZES zea
#else
#define ZES ze
#endif
c
#ifdef ORIGIN_THREAD
c$doacross share(latdon,lat1,lat2,anl,flp,flm,llens,rcs2,lcapdp,mype),
c$&        local(lat,lan,j,k,kaa,llensd)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(latdon,lat1,lat2)
cmic$1 shared(anl,flp,flm,llens,rcs2,lcapdp,mype)
cmic$1 private(lat,lan,j,k,kaa,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,kaa,llensd)
#endif
c
      do 2500 lat=lat1,lat2
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
c
        do k=1,LEVSS*2
          kaa=kaus-1+k
          do j=1,llensd*2
            anl(      j,kaa,lan)=anl(      j,kaa,lan)*rcs2(lat)
            anl(LCAPS+j,kaa,lan)=anl(LCAPS+j,kaa,lan)*rcs2(lat)
          enddo
        enddo
c
        call flpflm(flp(1,1,1,lan),flm(1,1,1,lan),anl(1,1,lan),
     1              llensd,lotas)
2500  continue
c
c no multi-threads should be given in following loop
c otherwise results will be non-reproduceble.
c
      do 3000 lat=lat1,lat2
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
c
#define DEFAULT
#ifdef FL2I
#undef DEFAULT
        call fl2i(flp(1,1,1,lan),flm(1,1,1,lan),ZES,qww(1,lat),
     1            llstr,llensd,lwvdef,lotas)
#endif
#ifdef DEFAULT
        call fl2ip(flp(1,1,1,lan),flm(1,1,1,lan),ZES,qww(1,lat),
     1             llstr,llensd,lwvdef,lotas)
#endif
c
        call uvsums(flp(1,1,kaus,lan),flm(1,1,kaus,lan),
     1              flp(1,1,kavs,lan),flm(1,1,kavs,lan),
     2              anltop(1,1,1),anltop(1,1,LEVSS+1),
     3              qvv(1,lat),wgt(lat),
     4              llstr,llensd,lwvdef,LEVSS)
c
3000  continue
#undef ZES
#endif                  /* not SCM */
c
#endif  /* endif for ifndef psplit */
#ifndef MP
      latdon=latdon+(lat2-lat1+1)
10000 continue
c
#endif
c
#ifdef PSPLIT
      return
#else
#ifndef SCM
c
c     input : ze=d(u)/d(t) di=d(v)/d(t)
c     output: uln=d(di)/d(t) vln=d(ze)/d(t)
c
#ifdef MP
#define UUS zea
#define VVS dia
#define XS ulna
#define WS vlna
#else
#define UUS ze
#define VVS di
#define XS uln
#define WS vln
#endif
      call uvtodz(UUS,VVS,XS,WS,anltop(1,1,1),
     1            anltop(1,1,LEVSS+1),llstr,llens,lwvdef)
#undef UUS
#undef VVS
#undef XS
#undef WS
#endif                     /* SCM */
c
#ifdef MP
#ifndef NISLQ
           klev=3+ntotal_   ! te,rq,uln,vln
#else
           klev=3           ! te,uln,vln
#endif
      call mpnk2nn(tea,lln22p_,levsp_, te,lnt22p_,levs_,klev)
#endif
c
      if(ladj) then
        do j=1,lnts2
          z(j)=0.0
        enddo
      endif
c
#ifdef ORIGIN_THREAD
c$doacross share(z,q,del,rq,rt,di,x,ze,w,te,y,uln,vln),
c$&        local(j,k)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(z,q,del,rq,rt,di,x,ze,w,te,y,uln,vln,lnts2)
cmic$1 private(j,k)
#endif
#ifdef OPENMP
!$omp parallel do private(j,k)
#endif
c
      do k=1,levs_
        do j=1,lnts2
          if(ladj) then
            di(j,k)=uln(j,k)-x(j,k)
            te(j,k)=te(j,k)-y(j,k)
          else
#ifdef SCM
            x(j,k)=di(j,k)
            y(j,k)=te(j,k)
          endif
          w(j,k)=ze(j,k)
#else
            x(j,k)=uln(j,k)
            y(j,k)=te(j,k)
          endif
          w(j,k)=vln(j,k)
#endif
        enddo
      enddo
c

#ifndef NISLQ
      do k=1,levh_
        do j=1,lnts2
          rt(j,k)=rq(j,k)
        enddo
      enddo
#endif

      if(ladj) then
        call impadj(x,y,q,di,te,z,uln,vln,lnts2,lnoffset)
      endif
c
#ifdef DG
#ifdef MP
      if(mype.eq.master) then
#endif
        print *,'(" gloopb t range ",2(4x,f6.1," @i,k,lat ",3i4))',
     &        tgmx,igmx,kgmx,jgmx,tgmn,igmn,kgmn,jgmn
#ifdef MP
      endif
#endif
#endif
c

c...................................................................
#ifdef CLM_CWF
#ifndef SCM
c
c ... find points within search radius & calculate wighting fuction ...
c
      if (stepone .and. lwconv) then
        print *, ' Find weight - lwconv ',lwconv
        nn = 67
        resol = real(jcap_)
        call fnam('../wgt_resol',12,resol,fno,ncho)
        open(unit=nn,file=fno(1:ncho),form='unformatted',
     +       status='old',err=900)
        write(*,*) ' Read Weighting File ',fno(1:ncho)
        read(nn) cwgt, inx, jnx, npnt
        close(nn)
        go to 910
900     write(*,*) ' Weighting File not Found ',fno(1:ncho)
        do j = 1, latg2_
#ifdef REDUCE_GRID
#ifdef MP
        lonsd2=lonfdp(j,mype)*2
#else
        lonsd2=lonfd(latdef(j))*2
#endif
#else
        lonsd2=lons2
#endif
        do i = 1, lonsd2
          xlam = xlon(i,j)
          ythe = xlat(i,j)
          xp = rerth_ * cos(xlam) * cos(ythe)
          yp = rerth_ * sin(xlam) * cos(ythe)
          zp = rerth_ * sin(ythe)
          call search(i,j,xlam,ythe,pi_,
     +                xp,yp,zp,lonf2_,latg2_,rerth_,xlon,xlat,
     +                cwgt,inx,jnx,npnt,npmx)
        enddo
        enddo
        open(unit=nn,file=fno(1:ncho),form='unformatted')
        write(nn) cwgt, inx, jnx, npnt
        close(nn)
910     lwconv = .false.
      endif
c
c ... find moisture convergence ...
      call rowsep(xlat)
      do 999 k = 1, levs_
        do j = 1, latg2_
#ifdef REDUCE_GRID
#ifdef MP
        lonsd2=lonfdp(j,mype)*2
#else
        lonsd2=lonfd(latdef(j))*2
#endif
#else
        lonsd2=lons2
#endif
        do i = 1, lonsd2
          wconv(i,j) = 0.0
        enddo
        enddo
c
c ... full grid for only qu & qv, but for conv ...
c
        call rowsep(qugs(1,k))
        call rowsep(qvgs(1,k))
        do j = 1, latg_
        do i = 1, lonf_
          ij = (j-1)*lonf_ + i
          qu(i,j) = qugs(ij,k)
          qv(i,j) = qvgs(ij,k)
        enddo
        enddo
        dellon = 2.e0*pi_ / lonf_
        do j = 2, latg_-1
        do i = 1, lonf_
          jmod = mod(j,2)
          if (jmod.eq.1) then
            ip =  i
          else
            ip = i + lonf_
          endif
          jp = int(j/2) + jmod
          jj = int(j/2) + jmod
          if (jmod.eq.1) then
            ii = i
            iip = i + lonf_
            iim = i + lonf_
            jjp = jj
            jjm = jj - 1
          else
            ii = i + lonf_
            iip = i
            iim = i
            jjp = jj + 1
            jjm = jj
          endif
          ythe = xlat(ii,jj)
          ythep = xlat(iip,jjp)
          ythem = xlat(iim,jjm)
          dellat = 0.5 * (ythem - ythep)
          dy = dellat * rerth_
          dx = rerth_ * cos(ythe) * dellon
          if (i.eq.lonf_) then
            uq2 = qu(1,j)
          else
            uq2 = qu(i+1,j)
          endif
          if (i.eq.1) then
            uq1 = qu(lonf_,j)
          else
            uq1 = qu(i-1,j)
          endif
          duq = (uq2-uq1)/(2.*dx)
          vq2 = qv(i,j-1)
          vq1 = qv(i,j+1)
          dvq = (vq2-vq1)/(2.*dy)
          conv(ip,jp) = -1.0*(duq + dvq)
        enddo
        enddo
c
        sum1 = 0.0
        sum2 = 0.0
        do i = 1, lonf_
         sum1 = sum1 + conv(i+lonf_,1)
         sum2 = sum2 + conv(i,int((latg_-1)/2)+1)
        enddo
        do i = 1, lonf_
          conv(i,1) = sum1/real(lonf_)
          conv(i+lonf_,int(latg_/2)) = sum2/real(lonf_)
        enddo
c
        call row1ns(conv)   ! LSMC for each grid point
c
c                        ... weighted mean ...
        do j = 1, latg2_
#ifdef REDUCE_GRID
#ifdef MP
        lonsd2=lonfdp(j,mype)*2
#else
        lonsd2=lonfd(latdef(j))*2
#endif
#else
        lonsd2=lons2
#endif
        do i = 1, lonsd2
          dsum = 0.0
          wsum = 0.0
          do n = 1, npnt(i,j)
            ii = inx(n,i,j)
            jj = jnx(n,i,j)
            dsum = dsum + conv(ii,jj)*cwgt(n,i,j)
            wsum = wsum + cwgt(n,i,j)
          enddo
          if (wsum.eq.0.0) then
            print *, 'WGT is zero at i=',i,' j=',j,' NPNT=',npnt(i,j)
          else
            wconv(i,j) = dsum/wsum
          endif
        enddo
        enddo
        do j = 1, latg2_
#ifdef REDUCE_GRID
#ifdef MP
        lonsd2=lonfdp(j,mype)*2
#else
        lonsd2=lonfd(latdef(j))*2
#endif
#else
        lonsd2=lons2
#endif
        do i = 1, lonsd2
          cgs(i,k,j) = wconv(i,j)
        enddo
        enddo
999   continue
      call row1ns(xlat)
c
#else                           /* SCM */
c
      lwconv = .false.
      do 999 k = 1, levs_
        do j = 1, latg2_
#ifdef REDUCE_GRID
#ifdef MP
        lonsd2=lonfdp(j,mype)*2
#else
        lonsd2=lonfd(latdef(j))*2
#endif
#else
        lonsd2=lons2
#endif
        do i = 1, lonsd2
          ij = (j-1)*latg2_ + i
          cgs(i,k,j) = -1.0*rt(ij,k)*wdiv(ij,k) + hadq(ij,k)
      print 1, k, j, i, ij, rt(ij,k), wdiv(ij,k), hadq(ij,k),cgs(i,k,j)
1     format('CWF-SCM',4i5,4e13.5)
        enddo
        enddo
999   continue
#endif                                  /* SCM */
#endif                                  /* CLM_CWF */
c
      return
#endif                                  /* PSPLIT */
      end
#ifdef CLM_CWF
c----------------------------------------------------------------------
      subroutine search(ip,jp,xlamp,ythep,ppi,
     +                  xp,yp,zp,ig,jg,earth,xlon,ylat,
     +                  cwgt,inx,jnx,npnt,npmx)
c
c     to find points and their weights inside the searching radius
c
c program history log:
c   04-04-10 young-hwa byun
c----------------------------------------------------------------------
#include <paramodel.h>
c
      dimension xlon(ig,jg), ylat(ig,jg)
      dimension cwgt(npmx,ig,jg),
     +          inx(npmx,ig,jg), jnx(npmx,ig,jg), npnt(ig,jg)
      dimension work(100000), iwrk(100000), jwrk(100000), rwrk(100000)
c
#ifdef REDUCE_GRID
#include <comreduce.h>
#include <comfgrid.h>
#endif
#ifdef MP
      lons2=lonlen(mype)*2
#else
      lons2=lonf2_
#endif
c
      er = 750000. ! Searching Radius = 750km
      np = 0
      do js = 1, jg
#ifdef REDUCE_GRID
#ifdef MP
        lonsd2=lonfdp(js,mype)*2
#else
        lonsd2=lonfd(latdef(js))*2
#endif
#else
        lonsd2=lons2
#endif
      do is = 1, lonsd2
        xlam = xlon(is,js)
        ythe = ylat(is,js)
        x = earth * cos(xlam) * cos(ythe)
        y = earth * sin(xlam) * cos(ythe)
        z = earth * sin(ythe)
        dx = xp-x
        dy = yp-y
        dz = zp-z
        r = sqrt(dx*dx+dy*dy+dz*dz)
        if (r.le.er) then
          np = np + 1
          iwrk(np) = is
          jwrk(np) = js
          t = r/er
          work(np) = exp(-t*t)
          rwrk(np) = r
        endif
      enddo
      enddo
c
      do i = 1, np
        k = i
        do 100 j = i+1,np
          if (rwrk(k).lt.rwrk(j)) go to 100
          k = j
100     continue
        ii = iwrk(i)
        iwrk(i) = iwrk(k)
        iwrk(k) = ii
        jj = jwrk(i)
        jwrk(i) = jwrk(k)
        jwrk(k) = jj
        ww = work(i)
        work(i) = work(k)
        work(k) = ww
        rr = rwrk(i)
        rwrk(i) = rwrk(k)
        rwrk(k) = rr
      enddo
c
      npnt(ip,jp) = np
      if (np.gt.npmx) npnt(ip,jp) = npmx
      do n = 1, npnt(ip,jp)
        cwgt(n,ip,jp) = work(n)
        inx(n,ip,jp) = iwrk(n)
        jnx(n,ip,jp) = jwrk(n)
      enddo
c
      return
      end
#endif                          /* CLM_CWF */
