#define G g
#define JCAP1 jcap1
#define JCAP1P jcap1p
#define JCAP jcap
#define JCAPP jcapp
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
#define LNT2P lnt2p
#define LNTP lntp
#define LONF22 lonf22
#define LONF22P lonf22p
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NCLDG ncldg
#define NCPUS ncpus
#define NPES npes
#define NTOTAL ntotal
#define NTRAC ntrac
#define RERTH rerth
#define TWOJ1 twoj1
      subroutine cof2grd(fgz,fq,fte,fdi,fvt,frq,
     &                   ggz,ggzx,ggzy,gq,gqx,gqy,
     1                   gte,guu,gvv,gdiv,gvot,grq)
c$$$  subprogram documentation block
c
c subprogram:    cof2grd  computes grid point values for post processing
c
c usage:    call cof2grd
c
c attributes:
c   language: fortran 90
c
c$$$
c
#include <paramodel.h>
#include <comio.h>
#include <comgpln.h>
#include <comfgrid.h>
#include <comfver.h>
c
#ifdef MP
#define JCAP1S jcap1p_
#define LNT2S lnt2p_
#define LNT22S lnt22p_
#define LONF2S lonf2p_
#define LATG2S latg2p_
#define LEVSS levsp_
#define LEVHS levhp_
#define LCAPS lcapp_
#define LCAP22S lcap22p_
#else
#define JCAP1S jcap1_
#define LNT2S lnt2_
#define LNT22S lnt22_
#define LONF2S lonf2_
#define LATG2S latg2_
#define LEVSS levs_
#define LEVHS levh_
#define LCAPS lonf_
#define LCAP22S lonf22_
#endif
c
c syn(1, 0*levs_+1, lan)  ze
c syn(1, 1*levs_+1, lan)  di
c syn(1, 2*levs_+1, lan)  te
c syn(1, 3*levs_+1, lan)  rq
c syn(1, 4*levs_+1, lan)  uln
c syn(1, 5*levs_+1, lan)  vln
c syn(1, 6*levs_+1, lan)  dpdphi
c syn(1, 6*levs_+2, lan)  dpdlam
c syn(1, 6*levs_+3, lan)  q
c syn(1, 6*levs_+4, lan)  dgzdphi
c syn(1, 6*levs_+5, lan)  dgzdlam
c syn(1, 6*levs_+8, lan)  gz
c
      parameter(lots    =5*levs_+levh_+6,lotst=2*levs_+2,
     &          ksz     =1,
     &          ksd     =1*levs_+1,
     &          kst     =2*levs_+1,
     &          ksr     =3*levs_+1,
     &          ksu     =3*levs_+levh_+1,kstb=3*levs_+levh_+1,
     &          ksv     =4*levs_+levh_+1,
     &          kspphi  =5*levs_+levh_+1,
     &          ksgzphi =5*levs_+levh_+2,
     &          ksplam  =5*levs_+levh_+3,
     &          ksp     =5*levs_+levh_+4,
     &          ksgzlam =5*levs_+levh_+5,
     &          ksgz    =5*levs_+levh_+6)
c
      parameter(lotss   =5*LEVSS+LEVHS+6,lotsts=2*LEVSS+2,
     &          kszs    =1,
     &          ksds    =1*LEVSS+1,
     &          ksts    =2*LEVSS+1,
     &          ksrs    =3*LEVSS+1,
     &          ksus    =3*LEVSS+LEVHS+1,kstbs=3*LEVSS+LEVHS+1,
     &          ksvs    =4*LEVSS+LEVHS+1,
     &          kspphis =5*LEVSS+LEVHS+1,
     &          ksgzphis=5*LEVSS+LEVHS+2,
     &          ksplams =5*LEVSS+LEVHS+3,
     &          ksps    =5*LEVSS+LEVHS+4,
     &          ksgzlams=5*LEVSS+LEVHS+5,
     &          ksgzs   =5*LEVSS+LEVHS+6)
c
#ifdef MP
#include <commpi.h>
c
      dimension ffa(lln22p_,lotss)
      dimension syf(lonf22_ ,lotss,latg2p_)
      dimension grs(lonf22p_,lots ,latg2p_)
#endif
      dimension ff (LNT22S,lots)
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
c
#ifdef MP
#define NCPUSS latg2_
#else
#define NCPUSS ncpus_
#endif
      dimension syn(LCAP22S,lotss,NCPUSS),
     &          syntop(2,JCAP1S,lotsts)
#undef NCPUSS
c
      dimension fgz(LNT22S),
     &          fq (LNT22S),
     &          fte(LNT22S,levs_),
     &          fdi(LNT22S,levs_),
     &          fvt(LNT22S,levs_),
     &          frq(LNT22S,levh_)
c
      dimension ggz (LONF2S,LATG2S),
     &          ggzx(LONF2S,LATG2S),
     &          ggzy(LONF2S,LATG2S),
     &          gq  (LONF2S,LATG2S),
     &          gqx (LONF2S,LATG2S),
     &          gqy (LONF2S,LATG2S),
     &          gte (LONF2S,LATG2S,levs_),
     &          guu (LONF2S,LATG2S,levs_),
     &          gvv (LONF2S,LATG2S,levs_),
     &          gdiv(LONF2S,LATG2S,levs_),
     &          gvot(LONF2S,LATG2S,levs_),
     &          grq (LONF2S,LATG2S,levh_)
c
      real, allocatable :: spec(:),specp(:)
      allocate (spec(lnt22_))
c
      ga2=g_/(rerth_*rerth_)
#ifdef MP
      call mpsp2f(fgz,lnt22p_,spec,lnt22_,1)
      if(mype.eq.master) then
      call spcshfli(spec,lnt22_,1,jcap_,lwvdef)
#else
      do n=1,lnt2_
        spec(n)=fgz(n)
      enddo
#endif
      spec(1)=z00
      spec(2)=0.
      do n=3,lnt2_
        spec(n)=spec(n)/snnp1(n)/ga2
      enddo
      do i=lnt2_+1,lnt22_
        spec(i)=0.
      enddo
#ifdef MP
      call spcshflo(spec,lnt22_,1,jcap_,lwvdef)
      endif
      allocate (specp(lnt22p_))
      call mpsf2p(spec,lnt22_,specp,lnt22p_,1)
#define SPEC specp
#else
#define SPEC spec
#endif
      do n=1,LNT22S
        ff(n,ksgz)=SPEC(n)
      enddo
      deallocate(spec)
#ifdef MP
      deallocate(specp)
#endif
#undef SPEC
c
      do n=1,LNT22S
        ff(n,ksp )=fq (n)
      enddo
      do k=1,levs_
        do n=1,LNT22S
          ff(n,kst+k-1)=fte(n,k)
          ff(n,ksd+k-1)=fdi(n,k)
          ff(n,ksz+k-1)=fvt(n,k)
        enddo
      enddo
      do k=1,levh_
        do n=1,LNT22S
          ff(n,ksr+k-1)=frq(n,k)
        enddo
      enddo
c
#ifdef MP
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
c
#ifdef MP
      call mpnn2nk(ff (1,ksz  ),lnt22p_,levs_,
     1             ffa(1,kszs ),lln22p_,levsp_,3+ntotal_)
      call mpnn2n (ff (1,ksp  ),lnt22p_,
     1             ffa(1,ksps ),lln22p_,1)
      call mpnn2n (ff (1,ksgz ),lnt22p_,
     1             ffa(1,ksgzs),lln22p_,1)
#define QS ffa(1,ksps)
#define DPDPHIS ffa(1,kspphis)
#define DPDLAMS ffa(1,ksplams)
#define GZS ffa(1,ksgzs)
#define DGZDPHIS ffa(1,ksgzphis)
#define DGZDLAMS ffa(1,ksgzlams)
#define DIS ffa(1,ksds)
#define ZES ffa(1,kszs)
#define ULNS ffa(1,ksus)
#define VLNS ffa(1,ksvs)
#else
#define QS ff(1,ksp)
#define DPDPHIS ff(1,kspphi)
#define DPDLAMS ff(1,ksplam)
#define GZS ff(1,ksgz)
#define DGZDPHIS ff(1,ksgzphi)
#define DGZDLAMS ff(1,ksgzlam)
#define DIS ff(1,ksd)
#define ZES ff(1,ksz)
#define ULNS ff(1,ksu)
#define VLNS ff(1,ksv)
#endif
      call dellnp(QS ,DPDPHIS ,
     &            syntop(1,1,2*LEVSS+1),DPDLAMS ,
     &            llstr,llens,lwvdef)
      call dellnp(GZS,DGZDPHIS,
     &            syntop(1,1,2*LEVSS+2),DGZDLAMS,
     &            llstr,llens,lwvdef)
      call dzuvle(DIS,ZES,ULNS,VLNS,
     &            syntop(1,1,1),
     &            syntop(1,1,LEVSS+1),
     &            llstr,llens,lwvdef)
#undef QS
#undef DPDPHIS
#undef DPDLAMS
#undef QLAPS
#undef DIS
#undef ZES
#undef ULNS
#undef VLNS
c
#ifndef MP
c compute latitude band limits
      last=mod(latg2_,ncpus_)
      nggs=(latg2_-last)/ncpus_
      if(last.ne.0)nggs=nggs+1
      inclat=ncpus_
      lat1=1-ncpus_
      lat2=0
      latdon=0
c
      do ngg=1,nggs
        if((ngg.eq.nggs).and.(last.ne.0)) inclat=last
        lat1=lat1+ncpus_
        lat2=lat2+inclat
#endif
c
#ifdef MP
        lat1=1
        lat2=latg2_
        latdon=0
#define ZES ffa(1,kszs)
#define TES ffa(1,ksts)
#else
#define ZES ff(1,ksz)
#define TES ff(1,kst)
#endif
c first lat loop
#ifdef ORIGIN_THREAD
c$doacross share(syntop,syn,qtt,qvv,lat1,lat2,latdon,
c$&        colrad,ze,zea,llstr,llens,lwvdef,lcapdp,mype,lcapd)
c$&        local(lat,lan,llensd)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(syntop,syn,qtt,qvv,lat1,lat2,latdon)
cmic$1 shared(colrad,ze,zea,llstr,llens,lwvdef,lcapdp,mype,lcapd)
cmic$1 private(lat,lan,llensd)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,llensd)
#endif
c
        do lat=lat1,lat2
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
          call sums2i(ZES,syn(1,1,lan),qtt(1,lat),
     &              llstr,llensd,lwvdef,lotss)
          call sumtop(syn(1,kstbs,lan),syntop,qvv(1,lat),
     &              llstr,llensd,lwvdef,lotsts)
        enddo 
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
c$&        local(lat,lan,k,j,lcapf,lonff)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$& shared(syf,syn)
cmic$& shared(lat1,lat2,latdon,latdef,lcapd,lonfd)
cmic$& private(lat,lan,k,j,lcapf,lonff)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,k,j,lcapf,lonff)
#endif
        do lat=lat1,lat2
#ifdef REDUCE_GRID
          lcapf=lcapd(latdef(lat))
          lonff=lonfd(latdef(lat))
#else
          lcapf=jcap1_
          lonff=lonf_
#endif
          lan=lat-latdon
          call ftilonf(SYNS(1,1,lan),SYNS(1,1,lan),lotss*2,
     1                 lcapf,lonff,latdef(lat),1)
        enddo
c
#undef SYNS
c
#ifdef MP
        call mpnk2nx(syf,lonf22_,lotss,
     &              grs,lonf22p_,lots,latg2p_,levsp_,levs_,1,1,
     &              5+ntotal_)
        call mpx2nx (syf,lonf22_,lotss,
     &              grs,lonf22p_,lots,latg2p_,kspphis,kspphi,6)
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
c
#ifdef ORIGIN_THREAD
c$doacross share(grs,syn,lat1,lat2,latdon,spdlat,
c$&              rbs2,del,rdel2,ci,p1,p2,h1,h2,tov,lons2,lonfdp,mype),
c$&              local(lat,lan,j,k,lonsd2)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(syn,grs,lat1,lat2,latdon,spdlat)
cmic$1 shared(rbs2,del,rdel2,ci,p1,p2,h1,h2,tov,lons2,lonfdp,mype)
cmic$1 private(lat,lan,j,k,lonsd2)
#endif
#ifdef OPENMP
!$omp parallel do private(lat,lan,j,k,lonsd2)
#endif
c
        do lat=lat1,lat2
          lan=lat-latdon
#ifdef MP
          rcs=sqrt(rbs2(lan))
#else
          rcs=sqrt(rbs2(lat))
#endif
#ifdef REDUCE_GRID
#ifdef MP
          lonsd2=lonfdp(lan,mype)*2
#else
          lonsd2=lonfd(latdef(lat))*2
#endif
#else
          lonsd2=lons2
#endif
c
#ifdef MP
          if( lonsd2.gt.0 ) then
#endif
#define LAN lan
#ifndef MP
#define LAN lat
#endif
            do i=1,lonsd2
              ggz (i,LAN)=SYNS(i,ksgz   ,lan)
              ggzx(i,LAN)=SYNS(i,ksgzlam,lan)*rcs
              ggzy(i,LAN)=SYNS(i,ksgzphi,lan)*rcs
              gq  (i,LAN)=SYNS(i,ksp    ,lan)
              gqx (i,LAN)=SYNS(i,ksplam ,lan)*rcs
              gqy (i,LAN)=SYNS(i,kspphi ,lan)*rcs
            enddo
            do i=lonsd2+1,LONF2S
              ggz (i,LAN)=ggz (lonsd2,lan)
              ggzx(i,LAN)=ggzx(lonsd2,lan)
              ggzy(i,LAN)=ggzy(lonsd2,lan)
              gq  (i,LAN)=gq  (lonsd2,lan)
              gqx (i,LAN)=gqx (lonsd2,lan)
              gqy (i,LAN)=gqy (lonsd2,lan)
            enddo
            do k=1,levs_
              do i=1,lonsd2
                gte (i,LAN,k)=SYNS(i,kst-1+k,lan)
                guu (i,LAN,k)=SYNS(i,ksu-1+k,lan)*rcs
                gvv (i,LAN,k)=SYNS(i,ksv-1+k,lan)*rcs
                gdiv(i,LAN,k)=SYNS(i,ksd-1+k,lan)
                gvot(i,LAN,k)=SYNS(i,ksz-1+k,lan)
              enddo
              do i=lonsd2+1,LONF2S
                gte (i,LAN,k)=gte (lonsd2,lan,k)
                guu (i,LAN,k)=guu (lonsd2,lan,k)
                gvv (i,LAN,k)=gvv (lonsd2,lan,k)
                gdiv(i,LAN,k)=gdiv(lonsd2,lan,k)
                gvot(i,LAN,k)=gvot(lonsd2,lan,k)
              enddo
            enddo
            do k=1,levh_
              do i=1,lonsd2
                grq (i,LAN,k)=SYNS(i,ksr-1+k,lan)
              enddo
              do i=lonsd2+1,LONF2S
                grq (i,LAN,k)=grq (lonsd2,lan,k)
              enddo
            enddo
c
#ifdef MP
          endif
#endif
        enddo
#ifndef MP
        latdon=latdon+(lat2-lat1+1)
      enddo
#endif
c
      return
      end
