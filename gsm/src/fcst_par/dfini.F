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
#define LNT2P lnt2p
#define LNTP lntp
#define LNUV lnuv
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define MTNVAR mtnvar
#define NCLDG ncldg
#define NTRAC ntrac
#define PI pi
#define TWOJ1 twoj1
       subroutine dfini(icall,hrini,chour,solsec,n1)
c.................................................................
c.....
c
#include <paramodel.h>
#include <comfibm.h>
#include <sfcfld.h>
c
#ifdef MP
#define LNT2S lnt2p_
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LNT2S lnt2_
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
      common/inigsm/   dthour , dshour, dchour, dsolsec, totsum
     2                        ,qs(LNT2S)
     2  ,tes(LNT2S,levs_),rqs(LNT2S,levh_)
     2  ,dis(LNT2S,levs_),zes(LNT2S,levs_)
c
      character*128 fno
      integer ncho
c
      if(numsum.ge.nummax) return
c
      if( icall.eq.0 ) then
#ifndef NOPRINT
      if (iope) then
      print *,' initial dfini '
      print *,' ini time is ',hrini,' hour.'
      endif
#endif
      do 10 k=1,levs_
      do 10 i=1,LNT2S
      dis(i,k) = 0.0
      zes(i,k) = 0.0
      tes(i,k) = 0.0
  10  continue
      do 11 k=1,levh_
      do 11 i=1,LNT2S
      rqs(i,k) = 0.0
  11  continue
      do 12 i=1,LNT2S
      qs(i) = 0.0
  12  continue
      totsum=0.0
      endif
c
      numsum=numsum+1
#ifndef NOPRINT
      if (iope) then
      print *,' ---- in dfini ---- numsum nummax ',numsum,nummax
      endif
#endif
      if( numsum.ne.0 ) then
        sc = pi_ / nummax
        sx= numsum*sc
        tx= numsum*pi_
        wx= tx/ ( nummax+1 )
        digfil= sin(wx)/wx * sin(sx)/tx
      else
        digfil = 1.0/nummax
      endif
      totsum = totsum + digfil
c
c
c------------------------do summation with window---
c
c first lat loop
#ifdef ORIGIN_THREAD
c$doacross share(di ,ze ,te ,rq , q,
c$&              dis,zes,tes,rqs, qs,
c$&              digfil),
c$&        local(j,k)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(di ,ze ,te ,rq , q )
cmic$1 shared(dis,zes,tes,rqs, qs)
cmic$1 shared(digfil)
cmic$1 private(j,k)
#endif
#ifdef OPENMP
!$omp parallel do private(j,k)
#endif
c      autoscope
c
c .......obtain full field values
      do 110 k=1,levs_
      do 110 j=1,LNT2S
      dis(j,k) = dis(j,k) + digfil*di(j,k)
      zes(j,k) = zes(j,k) + digfil*ze(j,k)
      tes(j,k) = tes(j,k) + digfil*te(j,k)
 110  continue
      do 111 k=1,levh_
      do 111 j=1,LNT2S
      rqs(j,k) = rqs(j,k) + digfil*rq(j,k)
 111  continue
      do 120 j=1,LNT2S
      qs(j) = qs(j) + digfil*q(j)
 120  continue
c................................................
c save
      if( numsum.eq.0 ) then
        dthour=thour
        dshour=shour
        dchour=chour
        dsolsec=solsec
#ifndef NOPRINT
        if (iope) then
        print *,' numsum=0, save thour= ',dthour
        endif
#endif
        call fnam('sfd',3,thour,fno,ncho)
        call fixio(n1,fno,sfcftyp,
     1             labs,idate(4),idate(2),idate(3),idate(1),thour,
     2             tsea,LONF2S,LATG2S,1)
      endif
c................................................
c restore
      if( numsum.eq.nummax ) then
#ifndef NOPRINT
        if (iope) then
        print *,' numsum=nummax reassign perturbation '
        print *,' with normalized factor=',totsum,' at hour=',dthour
c        print *,dshour,dchour,dsolsec
        endif
#endif
        hrini=0
        thour=dthour
        shour=dshour
        chour=dchour
        solsec=dsolsec
        call fnam('sfd',3,thour,fno,ncho)
        if (iope) then
        print*,'reloading ',fno(1:ncho)
        endif
        call fixio(n1,fno,sfcftyp,
     1             labs,idate(4),idate(2),idate(3),idate(1),thour,
     2             tsea,LONF2S,LATG2S,0)
#ifdef ORIGIN_THREAD
c$doacross share(di ,ze ,te ,rq , q,
c$&              dis,zes,tes,rqs, qs,
c$&              dim,zem,tem, rm, qm,
c$&              totsum),
c$&        local(j,k)
#endif
#ifdef CRAY_THREAD
cmic$ do all
cmic$1 shared(di ,ze ,te ,rq , q )
cmic$1 shared(dis,zes,tes,rqs, qs)
cmic$1 shared(dim,zem,tem, rm, qm)
cmic$1 shared(totsum)
cmic$1 private(j,k)
#endif
#ifdef OPENMP
!$omp parallel do private(j,k)
#endif
c      autoscope
        do 210 k=1,levs_
        do 210 j=1,LNT2S
        di (j,k) = dis(j,k) / totsum
        ze (j,k) = zes(j,k) / totsum
        te (j,k) = tes(j,k) / totsum
        dim(j,k) = dis(j,k) / totsum
        zem(j,k) = zes(j,k) / totsum
        tem(j,k) = tes(j,k) / totsum
 210    continue
        do 211 k=1,levh_
        do 211 j=1,LNT2S
        rq (j,k) = rqs(j,k) / totsum
        rm (j,k) = rqs(j,k) / totsum
 211    continue
        do 220 j=1,LNT2S
        qm(j) = qs(j) / totsum
        q (j) = qs(j) / totsum
 220    continue
      do 230 l=1,LATG2S
      do 230 j=1,LONF2S
      geshem(j,l)=0.5*geshem(j,l)
230   continue
      call mltflx(0.5,dusfc,dvsfc,dtsfc,dqsfc,dlwsfc,ulwsfc,
     1 bengsh,gflux,
     2 dugwd,dvgwd,psmean)
#ifdef DG3
      call mltdia(0.5)
#endif
#ifdef DGZ
      call znlmlt(0.5)
#endif
      do i=1,LONF2S
        do j=1,LATG2S
          do k=1,25
            fluxr(i,j,k)=0.5*fluxr(i,j,k)
          enddo
        enddo
      enddo
      do i=1,LONF2S
        do j=1,LATG2S
          cvavg(i,j)=0.5*cvavg(i,j)
        enddo
      enddo
      endif
c................................................
c
      return
      end
