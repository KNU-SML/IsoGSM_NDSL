#define G g
#define IGEN igen
#define JCAP jcap
#define LATG latg
#define LEVH levh
#define LEVP1 levp1
#define LEVS levs
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LONF lonf
#define NCLDG ncldg
#define NPES npes
#define NTRAC ntrac
#define RERTH rerth
#include <define.h>
#ifndef HYBRID      
      subroutine rdsig(n,fhour,idate,gz,q,te,di,ze,rq,sl,si,z00)
#else
      subroutine rdsig(n,fhour,idate,gz,q,te,di,ze,rq,ak5,bk5,z00)
#endif
c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    rdsig       reads sigma level spectral coefficients.
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-27
c
c abstract: reads a complete set of sigma level spectral coefficients
c   at a single time to be used to start the model forecast.
c   the subroutine compares si and sl (the models vertical
c   structure) computed in setsig with the si and sl of the
c   input coefficients in order to make sure the coefficients
c   were generated under the same vertical structure.
c
c program history log:
c   88-04-27  joseph sela
c   99-03-02  hann-ming henry juang  for mpi
c   00-05-02  song-you hong  for tracer and prognostic clouds
c
c usage:    call rdsig (n,fhour,idate,gz,q,te,di,ze,rq,sl,si,z00)
c   input argument list:
c     n        - fortran unit number for file containing the
c                complete set of sigma level coefficients.
c     sl       - sigma layers     computed in setsig.
c     si       - sigma interfaces computed in setsig.
c
c   output argument list:
c     fhour    - forecast hour of the set of coefficients read
c                from unit n.
c     idate    - idate(1)=initial hour (gmt) of forecast from
c                         which coefficients were made.
c                idate(2)=month (1-12).
c                idate(3)=day of the month.
c                idate(4)=year of the century.
c     gz       - laplacian of topography.
c                gz is multiplied by the constant snnp1 array.
c                gz is then multiplied by the gravity constant and
c                divided by the square of the radius of the earth.
c     q        - ln(psfc)          coefficients.
c     te       - temperature       coefficients.
c     di       - divergence        coefficients.
c     ze       - vorticity         coefficients.
c     rq       - specific humidity coefficients.
c     z00      - mean topography.
c                z00 is set equal to gz(1) after gz is read.
c
c   input files:
c     unit n   - complete set of sigma level spectral coefficients.
c
c   output files:
c     output   - print file.
c
c attributes:
c   language: fortran 200.
c   machine:  cyber 205.
c
c$$$
#include <paramodel.h>
#include <comio.h>
#ifdef MP
#include <commpi.h>
      real, allocatable ::
     1      spec(:),coef1(:,:),coef2(:,:)
#endif
#ifdef MP
#define LNT22S lnt22p_
#else
#define LNT22S lnt22_
#endif
      dimension
     1 gz(LNT22S),
     2  q(LNT22S),
     3 te(LNT22S,levs_),
     4 di(LNT22S,levs_),
     5 ze(LNT22S,levs_),
     6 rq(LNT22S,levh_)
      dimension idate(4),si(levp1_),sl(levs_)
      dimension          xi(levp1_),xl(levs_)
      dimension          ak5(levp1_),bk5(levp1_)
      dimension          xak(levp1_),xbk(levp1_)
c
      parameter(kdum=201-levs_-1-levs_,kdum2=21,kens=2)
      dimension dummy(kdum),dummy2(kdum2),ensemble(kens)
c
      dimension work(LNT22S,levs_)
c
c  addition by mk 4/5/90
c
      save ifp
      data ifp/0/
c
      if(ifp.eq.0) then
        call cmpind
        ifp=1
      endif
c
c  end addition
c
#ifdef MP
      allocate (spec(lnt22_))
      allocate (coef1(lnt22_,levs_))
      allocate (coef2(lnt22_,levs_))
#endif
c
c     spectral data file format
c     lab
c     hour,idate(4),si(levp1_),sl(levs_)
c     zln q te di ze
c
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef ASSIGN
       call assign('assign -R')
#endif
        rewind n
        read(n)lab
#ifndef NOPRINT
        print 3000,lab,n
3000    format(1x,'rdsig lab  ',4a10,' n=',i3)
#endif
#ifdef MP
      endif
#endif
c
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifndef HYBRID
      read(n,err=201)fhour,idate,(xi(k),k=1,levp1_),(xl(k),k=1,levs_)
#else
      read(n,err=201)fhour,idate,(xak(k),k=1,levp1_),(xbk(k),k=1,levp1_)
#endif
     &       ,dummy,waves,xlayers,trun,order,realform,gencode
     &       ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers
     &       ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid
     &       ,pdryini,dummy2,clouds
      if(tracers.eq.0.) then
        tracers=1.
        print *,'tracers reset to 1.'
      endif
      itrace=nint(tracers)
      icloud=nint(clouds)
#ifndef NOPRINT
      print *,'rdsig unit,fhour,idate=',n,fhour,idate
      print *,' number of tracers input = ',itrace
      print *,' number of cloud input   = ',icloud
#endif
      goto 202
201   continue
      rewind n
      read(n) lab
#ifdef MP
c
#endif
#ifndef HYBRID
      read(n,err=201)fhour,idate,(xi(k),k=1,levp1_),(xl(k),k=1,levs_)
#else
      read(n,err=201)fhour,idate,(xak(k),k=1,levp1_),(xbk(k),k=1,levp1_)
#endif
      do i=1,kdum
        dummy(i)=0.
      enddo
      waves=jcap_
      xlayers=levs_
      trun=1.
      order=2.
      realform=1.
      gencode=igen_
      rlond=lonf_
      rlatd=latg_
      rlonp=lonf_
      rlatp=latg_
      rlonr=lonf_
      rlatr=latg_
      tracers=1.
      clouds=0.
      pdryini=0.
      subcen=0.
      do i=1,kens
        ensemble(i)=0.
      enddo
      ppid=0.
      slid=0.
      vcid=0.
      vmid=0.
      vtid=0.
      do k=1,kdum2
        dummy2(k)=0.
      enddo
      itrace=nint(tracers)
      icloud=nint(clouds)
#ifndef NOPRINT
      print *,'rdsig old format unit,fhour,idate=',n,fhour,idate
      print *,' number of tracers input = ',itrace
      print *,' number of cloud input   = ',icloud
#endif
202   continue
c
#ifndef DRY_MODEL
      if(itrace.le.0) then
        print *,'moisture not predicted'
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
#else
      itrace=0
#endif
c
#ifdef MP
      endif
      call mpbcastr(fhour,1)
      call mpbcasti(idate,4)
#ifndef HYBRID
      call mpbcastr(xi,levp1_)
      call mpbcastr(xl,levs_)
#else
      call mpbcastr(xak,levp1_)
      call mpbcastr(xbk,levp1_)
#endif
      call mpbcasti(itrace,1)
      call mpbcasti(icloud,1)
#endif
c
      ga2=g_/(rerth_*rerth_)
#ifdef MP
      if( mype.eq.master ) then
        read(n)(spec(i),i=1,lnt2_)
        z00=spec(1)
#ifndef NOPRINT
        print *,' rdsig gz z00=',z00
#endif
        call spcshfli(spec,lnt22_,1,jcap_,lwvdef)
        do j=1,lnt2_
          spec(j)=spec(j)*snnp1(j)*ga2
        enddo
        call spcshflo(spec,lnt22_,1,jcap_,lwvdef)
      endif
      call mpbcastr(z00,1)
      call mpsf2p(spec,lnt22_,gz,lnt22p_,1)
#else
      read(n)(gz(i),i=1,lnt2_)
      z00=gz(1)
      do 300 j=1,lnt2_
        gz(j)=gz(j)*snnp1(j)*ga2
300   continue
#endif
c
#ifdef MP
      if( mype.eq.master ) then
        read(n)(spec(i),i=1,lnt2_)
#ifndef NOPRINT
        print *,' rdsig q '
#endif
      endif
      call mpsf2p(spec,lnt22_,q,lnt22p_,1)
#else
      read(n)(q(i),i=1,lnt2_)
#endif
c
#ifdef MP
      if( mype.eq.master ) then
        do k=1,levs_
          read(n)(coef1(i,k),i=1,lnt2_)
        enddo
#ifndef NOPRINT
        print *,' rdsig te '
#endif
      endif
      call mpsf2p(coef1,lnt22_,te,lnt22p_,levs_)
#else
      do 400 k=1,levs_
      read(n)(te(i,k),i=1,lnt2_)
400   continue
#endif
c
#ifdef MP
      if( mype.eq.master ) then
        do k=1,levs_
          read(n)(coef1(i,k),i=1,lnt2_)
          read(n)(coef2(i,k),i=1,lnt2_)
        enddo
#ifndef NOPRINT
        print *,' rdsig di ze '
#endif
      endif
      call mpsf2p(coef1,lnt22_,di,lnt22p_,levs_)
      call mpsf2p(coef2,lnt22_,ze,lnt22p_,levs_)
#else
      do 401 k=1,levs_
      read(n)(di(i,k),i=1,lnt2_)
      read(n)(ze(i,k),i=1,lnt2_)
401   continue
#endif
c
#ifdef MP
      if( mype.eq.master ) then
#endif
        do k=1,levs_
#ifdef MP
          read(n)(coef1(i,k),i=1,lnt2_)
        enddo
#ifndef NOPRINT
         print *,' rdsig rq '
#endif
      endif
      call mpsf2p(coef1,lnt22_,rq,lnt22p_,levs_)
#else
          read(n)(rq(i,k),i=1,lnt2_)
        enddo
#endif
       if(itrace.ge.2) then
         do nt = 2,itrace
           ktrace = (nt-1)*levs_+1
#ifdef MP
           if( mype.eq.master ) then
                do k=1,levs_
                  read(n)(coef1(i,k),i=1,lnt2_)
                enddo
           endif
           if(ntrac_.ge.nt) then
             call mpsf2p(coef1,lnt22_,rq(1,ktrace),lnt22p_,levs_)
           endif
#else
           do k=1,levs_
             read(n)(work(i,k),i=1,lnt2_)
           enddo
           if(ntrac_.ge.nt) then
             do k=1,levs_
               do i=1,lnt2_
                 rq(i,k+ktrace-1) = work(i,k)
               enddo
             enddo
           endif
#endif
         enddo
       endif
       if(ntrac_.gt.itrace) then
         do nt = itrace+1,ntrac_
           ktrace = (nt-1)*levs_+1
           do  k = ktrace,ktrace+levs_-1
             do  i = 1,LNT22S
               rq(i,k) = 0.0
             enddo
           enddo
         enddo
#ifdef CO2PROG
** co2 uniform intial **
           nt=3
           ktrace = (nt-1)*levs_+1
#ifdef MP
           if( mype.eq.master ) then
             do k=1,levs_
               do i=1,lnt2_
                 coef1(i,k)=0.
               enddo
               coef1(1,k)=376.e-6 * 44. / 29. * 2.**0.5 !! ppmv=moleco2/moleair
             enddo              
           endif
           call mpsf2p(coef1,lnt22_,rq(1,ktrace),lnt22p_,levs_)
#else
           do k=1,levs_
             do i=1,lnt2_
               rq(i,k+ktrace-1) = 0.
             enddo
             rq(1,k+ktrace-1)=376.e-6 * 44. / 29. * 2.**0.5 !! ppmv=moleco2/moleair
           enddo
#endif
#endif
       endif
c
       kcloud = ntrac_*levs_ + 1
       if(icloud.ge.1) then
         do nc = 1,icloud
           kc = (nc-1)*levs_+kcloud
#ifdef MP
       if( mype.eq.master ) then
         do k=1,levs_
           read(n)(coef1(i,k),i=1,lnt2_)
         enddo
       endif
         if(ncldg_.ge.nc) then
           call mpsf2p(coef1,lnt22_,rq(1,kc),lnt22p_,levs_)
         endif
#else
         do k=1,levs_
           read(n)(work(i,k),i=1,lnt2_)
         enddo
         if(ncldg_.ge.nc) then
           do k=1,levs_
             do i=1,lnt2_
               rq(i,k+kc-1) = work(i,k)
             enddo
           enddo
         endif
#endif
         enddo
       endif
       if(ncldg_.gt.icloud) then
         do nc = icloud+1,ncldg_
           ncloud = (nc-1)*levs_ + kcloud
           do  k = ncloud,ncloud+levs_-1
             do  i = 1,LNT22S
               rq(i,k) = 0.0
             enddo
           enddo
         enddo
        endif
c
#ifdef MP
      deallocate (spec)
      deallocate (coef1)
      deallocate (coef2)
#endif
c
      return
      end
