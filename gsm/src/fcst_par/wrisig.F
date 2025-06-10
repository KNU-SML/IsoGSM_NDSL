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
#define NROW nrow
#define NTOTAL ntotal
#define NTRAC ntrac
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp

#include <define.h>
#ifndef HYBRID
      subroutine wrisig(n,fhour,idate,q,te,di,ze,rq,sl,si,gz,z00,
#else
      subroutine wrisig(n,fhour,idate,q,te,di,ze,rq,ak5,bk5,gz,z00,
#endif
#ifdef SCM
     1                  vvel,
#endif
     1                  itpdt)

c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    wrisig      writes sigma level spectral coefficients.
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-29
c
c abstract: writes a complete set of forecast sigma level
c   spectral coefficients for all model variables.
c
c program history log:
c   88-04-29  joseph sela
c   88-11-02  mark rozwodoski  changed second record to selalabel.
c
c usage:    call wrisig (n,fhour,idate,z,q,te,di,ze,rq,sl,si,gz,z00)
c   input argument list:
c     n        - fortran unit number for file to be written to.
c     fhour    - forecast hour of the set of coefficients.
c     idate    - idate(1)=initial hour (gmt) of forecast.
c                idate(2)=month (1-12).
c                idate(3)=day of the month.
c                idate(4)=year of the century.
c     q        - ln(psfc)          coefficients.
c     te       - temperature       coefficients.
c     di       - divergence        coefficients.
c     ze       - vorticity         coefficients.
c     rq       - specific humidity coefficients.
c     sl       - sigma layers     computed in setsig.
c     si       - sigma interfaces computed in setsig.
c     gz       - laplacian of topography.
c     z00      - mean topography.
c     itpdt    - flag for time level t, t+dt or initialized (1, 2 or 3)
c                4 for sigit, 5 for sigitdt file names
c
c   output argument list:
c     z        - topography spectral coefficients.
c
c   output files:
c     unit n   - complete set of forecast sigma level spectral
c                coefficients.
c     output   - print file.
c
c attributes:
c   language: fortran 200.
c   machine:  cyber 205.
c
c$$$
#include <paramodel.h>
#include <comio.h>
#ifdef NISLQ
#include <comnislq.h>
#endif
#ifdef MP
#include <commpi.h>
      real, allocatable ::
     1      coef1(:,:),coef2(:,:)
#endif
#ifdef MP
#define LNT22S lnt22p_
#else
#define LNT22S lnt22_
#endif
#ifdef SCM
      dimension vvel(lnt22_,levs_)
#endif
      dimension gz(LNT22S),q(LNT22S),spec(lnt22_),
     1          te(LNT22S,levs_),di(LNT22S,levs_),ze(LNT22S,levs_),
     2          rq(LNT22S,levh_)
      dimension idate(4),si(levp1_),sl(levs_)
      dimension ak5(levp1_),bk5(levp1_)
      dimension ak5x(levp1_),bk5x(levp1_)
c
      dimension orog(lnt22_)
c
#ifndef HYBRID      
      parameter(kdum=201-levs_-1-levs_,kdum2=21,kens=2)
#else
      parameter(kdum=201-levs_-2-levs_,kdum2=21,kens=2)
#endif
      dimension dummy(kdum),dummy2(kdum2),ensemble(kens)
c
      dimension idatex(4)
c
      character*128 fni(5)
      dimension nchi(5)
      data fni/'sig','sigp','sigi','sigit','sigitdt'/
      data nchi/3,4,4,5,7/
      character*128 fno

#ifdef NISLQ
      if (itpdt .eq. 1 .or. itpdt .eq. 4 ) then
        call sph_fft_driver(1,slq_q1,rq,ntotal_)
      else if (itpdt.eq.2 .or. itpdt.eq.5) then
        call sph_fft_driver(1,slq_q2,rq,ntotal_)
      endif
#endif





c
c  get topography from history file
c
#ifdef MP
      if( mype.eq.master ) then
#endif
      close(n)
      rewind n
      open(unit=n,file='./sigit',form='unformatted',err=700)
      go to 701
  700 continue
      write(6,*) ' error in opening file sigit at wrisig'
#ifdef MP
#ifdef RSM
                call rmpabort
#else
                call mpabort
#endif
#else
                call abort
#endif
  701 continue
      rewind n
      read(n)
#ifndef NOPRINT
      print *,' read lab '
#endif
      read(n)
#ifndef NOPRINT
      print *,' read fhour idate  '
#endif
      read(n)( orog(i),i=1,lnt2_)
#ifndef NOPRINT
      print *,' read gz '
#endif
      close(n)
#ifdef MP
      endif
#endif
c
#ifdef MP
      if( mype.eq.master ) then
#endif
      if(itpdt.lt.4) then
        call fnam(fni(itpdt),nchi(itpdt),fhour,fno,ncho)
      else
        fno=fni(itpdt)
        ncho=nchi(itpdt)
      endif
c
#ifdef ASSIGN
      call assign('assign -R')
#endif
      close(n)
      open(unit=n,file=fno(1:ncho),form='unformatted',err=900)
      go to 901
  900 continue
      write(6,*) ' error in opening file at wrisig ',fno(1:ncho)
#ifdef MP
      call mpabort
#else
      call abort
#endif
  901 continue
#ifndef NOPRINT
      write(6,*) ' file ',fno(1:ncho),' opened. unit=',n
#endif
#ifdef MP
      endif
      allocate (coef1(lnt22_,levs_))
      allocate (coef2(lnt22_,levs_))
#endif
c
#ifdef MP
      if( mype.eq.master ) then
#endif
      rewind(n)
      write(n) lab
#ifndef NOPRINT
      print 3000,lab,n
3000  format(1x,' wrisig lab ',4a10,' n=',i3)
#endif
      do k=1,kdum
        dummy(k)=0.
      enddo
      waves=jcap_
      xlayers=levs_
      trun=1.
      order=2.
      realform=1.
      gencode=igen
      rlond=lonf_
      rlatd=latg_
      rlonp=lonf_
      rlatp=latg_
      rlonr=lonf_
      rlatr=latg_
      tracers=ntrac_
      clouds=ncldg_
      pdryini=1.
      subcen=icen2
      ensemble(1)=ienst
      ensemble(2)=iensi
      ppid=0.
      slid=0.
      vcid=0.
      vmid=0.
      vtid=0.
      do k=1,kdum2
        dummy2(k)=0.
      enddo
#ifdef RESET_FHOUR
      call vday(idate(4),idate(2),idate(3),idate(1),fhour,
     *          jy,jm,jd,jh,rjday)
      if (jy.eq.idate(4)+1.and.jm.eq.1.and.jd.eq.1.and.jh.eq.0 ) then
        fxhour=0.
        idatex(4)=jy
        idatex(2)=1
        idatex(3)=1
        idatex(1)=0
        print *,'fhour reset to zero.'
      else
        fxhour=fhour
        idatex(4)=idate(4)
        idatex(2)=idate(2)
        idatex(3)=idate(3)
        idatex(1)=idate(1)
      endif
#else
      fxhour=fhour
      idatex(4)=idate(4)
      idatex(2)=idate(2)
      idatex(3)=idate(3)
      idatex(1)=idate(1)
#endif
#ifdef HYBRID
      do k=1,levp1_
        ak5x(k)=ak5(levp1_+1-k)*1000. ! cb -> Pa
        bk5x(k)=bk5(levp1_+1-k)
      enddo
#endif
      if(ntrac_.eq.1.and.ncldg_.eq.0) then
#ifndef HYBRID
      write(n)fxhour,idatex,si,sl
#else
      write(n)fxhour,idatex,ak5x,bk5x
#endif
     &       ,dummy,waves,xlayers,trun,order,realform,gencode
     &       ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers
     &       ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid
     &       ,dummy2
      else
#ifndef HYBRID
      write(n)fxhour,idatex,si,sl
#else
      write(n)fxhour,idatex,ak5x,bk5x
#endif
     &       ,dummy,waves,xlayers,trun,order,realform,gencode
     &       ,rlond,rlatd,rlonp,rlatp,rlonr,rlatr,tracers
     &       ,subcen,ensemble,ppid,slid,vcid,vmid,vtid,runid,usrid
     &       ,pdryini,dummy2,clouds
      endif
#ifndef NOPRINT
      print *,' wrisig fhour idate ... '
#endif
c gz
      write(n)( orog(i),i=1,lnt2_)
#ifndef NOPRINT
      print *,' wrisig gz '
#endif
#ifdef MP
      endif
#endif
c  q
#ifdef MP
      call mpsp2f(q,lnt22p_,spec,lnt22_,1)
      if( mype.eq.master ) then
        write(n)(spec(i),i=1,lnt2_)
#ifndef NOPRINT
        print *,' wrisig q '
#endif
      endif
#else
      write(n)( q(i),i=1,lnt2_)
#endif
c  te
#ifdef MP
      call mpsp2f(te,lnt22p_,coef1,lnt22_,levs_)
      if( mype.eq.master ) then
        do k=1,levs_
          write(n)(coef1(i,k),i=1,lnt2_)
        enddo
#ifndef NOPRINT
        print *,' wrisig te '
#endif
      endif
#else
      do 15 k=1,levs_
      write(n)(te(i,k),i=1,lnt2_)
15    continue
#endif
c  di ze
#ifdef MP
      call mpsp2f(di,lnt22p_,coef1,lnt22_,levs_)
      call mpsp2f(ze,lnt22p_,coef2,lnt22_,levs_)
      if( mype.eq.master ) then
        do k=1,levs_
          write(n)(coef1(i,k),i=1,lnt2_)
          write(n)(coef2(i,k),i=1,lnt2_)
        enddo
#ifndef NOPRINT
        print *,' wrisig di ze '
#endif
      endif
#else
      do k=1,levs_
        write(n)(di(i,k),i=1,lnt2_)
        write(n)(ze(i,k),i=1,lnt2_)
      enddo
#endif
c  rq
      do nt = 1,ntotal_
      kt = (nt-1)*levs_+1
#ifdef MP
      call mpsp2f(rq(1,kt),lnt22p_,coef1,lnt22_,levs_)
      if( mype.eq.master ) then
        do k=1,levs_
          write(n)(coef1(i,k),i=1,lnt2_)
        enddo
#ifndef NOPRINT
        print *,' wrisig rq '
#endif
      endif
#else
      do k=1,levs_
        write(n)(rq(i,kt+k-1),i=1,lnt2_)
      enddo
#endif
      enddo
c
#ifdef SCM
c --------------------------------- vvel
      do k=1,levs_
        write(n)(vvel(i,k),i=1,lnt2_)
      enddo
#endif
c ---------------------------------
c
#ifndef NOPRINT
#ifdef MP
      if( mype.eq.master ) then
#endif
      close(n)
      print 3001,fhour,idate,n
      print *,'number of tracer = ',tracers
      print *,'number of clouds = ',clouds
3001  format(1x,'wrisig fhour=',f10.1,2x,4i6,2x,'n=',i2)
#ifdef MP
      endif
#endif
#endif
c
#ifdef MP
      deallocate (coef1)
      deallocate (coef2)
#endif
c
      return
      end
