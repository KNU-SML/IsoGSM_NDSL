#define IGEN igen
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
#define LNT2P lnt2p
#define LNTP lntp
#define LNUT2 lnut2
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
#define RERTH rerth
#define SLVARK slvark
#define TWOJ1 twoj1
#include <define.h>
      subroutine getcon(n1,
#ifdef RSM
     3 nrsmi1,nrsmi2,nrflip,
     & nrsmo1,nrsmo2,nrflop,nrsfli,nrsflx,nrinit,nrpken,
#endif
#ifdef LFM
     & nlfmsgi,nlfmsfi,nlfmsgo,nlfmsfo,klenp,weight,fwindow,
#endif
#ifdef COUPLE
     5 kpfix,ksfcx,ksig,ksfc,kpost,krestart,klfm,krsm,kcouple)
#else
     5 kpfix,ksfcx,ksig,ksfc,kpost,krestart,klfm,krsm)
#endif
#include <paramodel.h>
c....
c
c $Id: getcon.F,v 1.37 2011/09/27 08:30:14 kyoshimura Exp $
c
#include <comfibm.h>
#include <comcon.h>
#include <comgpd.h>
#ifdef MP
#include <commpi.h>
      dimension ilist(100),rlist(100)
      dimension rdex(lnt2_)
      dimension grid1(lonf2_,latg2_)
      dimension grid2(lonf2_,latg2_)
      dimension gridv(lonf2_,latg2_,mtnvar_)
      dimension tmpq(lnt2p_)
#endif
#ifdef REDUCE_GRID
#include <comreduce.h>
      data ndigit/1/
c      data ndigit/4/
#endif
      dimension tmpqtt(lnt2_),tmpqdd(lnt2_),tmpqvv(lnut2_)
c
#ifdef SCM
      common /comscmbdy/dtbdy, curtime, vvel(lnt22_,levs_)
#ifdef CLM_CWF
     +                , wdiv(lnt22_,levs_), hadq(lnt22_,levs_)
#endif
#endif
c
#ifdef LFM
      dimension weight(*)
      dimension weix(1000)
#endif
c
      namelist/namsmf/ con,num,labl,endhour,ldebug,filta,icen,igen,icen2
     &                ,ienst,iensi,runid,usrid,ncpus
     &                ,intpfix,intsfcx,intsig,intsfc,intpost,intrestart
     &                ,intlfm,intrsm
     &                ,ndigit,relx
#ifdef LFM
     &                ,critfs,filtwin
#endif
#ifdef COUPLE
     &                ,intcouple
#endif
#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
       limlow=1
       jcap=jcap_
       levs=levs_
       filta= 0.92
       dt80=939.14e0
       percut=27502.e0
       icen=7
       igen=igen_
       icen2=0
       ienst=0
       iensi=0
       runid=0
       usrid=0
       call gncpus(ncpus)
c..........................................
c
c  define unit numbers
c
c  input
      n1    = 11
c     nmtnv = 24
      nmtnv = 25
#ifdef LFM
      nlfmsgi=30
      nlfmsfi=31
#endif
#ifdef RSM
      nrsmi1= 30
      nrsmi2= 31
      nrflip= 32
#endif
#ifdef LFM
      nlfmsgo= 70
      nlfmsfo= 71
#endif
#ifdef RSM
      nrsmo1= 70
      nrsmo2= 71
      nrflop= 72
      nrsflx= 73
      nrsfli= 74
      nrpken= 77
      nrinit= 78
#endif
c
c  output interval control
c
      intpfix=24
      intsfcx=24
      intsig=24
      intsfc=24
      intpost=24
      intrestart=24
      intlfm=24
      intrsm=24
      intcouple=6
c....
c.... cmean,clstp control time averaging of convective clds in kuo
c....
      clstp=99.
c....
c...  averaging interval for conv cld approx 3 hrs (num of timesteps)
c....
#ifndef NOPRINT
      if(mype.eq.master) print 100, jcap, levs
100   format (1h0,'getcon ',i3,i3,'created april 92')
#endif
      filtb =(1.e0-filta) * 0.5e0
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef SCM
      open(unit=n1,file='basic.info',form='unformatted',status='old')
      read(n1) dtbdy, xlon, xlat, lsmask, ci, si, del, sl, cl, rpi
      close(n1)
      rad = pi_/180.
      do i = 1, lonf2_
      do j = 1, latg2_
        xlon(i,j) = xlon(i,j)*rad
        xlat(i,j) = xlat(i,j)*rad
        print *, xlon(i,j), xlat(i,j)
      enddo
      enddo
#else
#ifndef HYBRID
      call setsig(ci,si,del,sl,cl,rpi,n1)
#else
      call sethyb(ak5,bk5,ci,si,del,sl,cl,rpi,n1)
#endif
#endif
#ifdef MP
      endif
      call mpbcastr(ci,levs_)
      call mpbcastr(si,levs_)
      call mpbcastr(del,levs_)
      call mpbcastr(sl,levs_)
      call mpbcastr(cl,levs_)
      call mpbcastr(rpi,levs_)
#ifdef HYBRID
      call mpbcastr(ak5,levp1_)
      call mpbcastr(bk5,levp1_)
#endif
#endif
!! only when ak5(1)=ak5(2)=0.  
      sl1=sl(1)
      do 3 lev=1,levs_
      tov(lev)=300.e0
3     continue
#ifdef SCM
      rcs2(1) = 1.0
      wgt(1) = 1.0
#else
#ifndef HYBRID
      call ambmsv(levs_,si,sl,tov,am,bm,sv,gv,cm)
#else
      call ambmsv_hyb(levs_,ak5,bk5,am,bm,sv,gv,cm)
#endif      
      call glats(latg2_, colrad, wgt, wgtcs, rcs2)
#endif
#ifdef MP
      do j=1,latlen(mype)
        jj=latdef(latstr(mype)+j-1)
        colrab(j)=colrad(jj)
        wgb   (j)=wgt   (jj)
        wgbcs (j)=wgtcs (jj)
        rbs2  (j)=rcs2  (jj)
      enddo
#else
      do j=1,latg2_
        jj=latdef(j)
#ifndef SCM
        colrab(j)=colrad(jj)
        wgb   (j)=wgt   (jj)
        wgbcs (j)=wgtcs (jj)
#endif
        rbs2  (j)=rcs2  (jj)
      enddo
#endif
#ifndef SCM
      call epslon(eps, jcap_)
#endif

#ifdef NISLQ
!
! initialize nislq
!
#ifdef MP
      call nislq_init(nrow,myrow,colrad,rbs2)
#else
      call nislq_init(1,0,colrad,rbs2)
#endif
#endif /* NISLQ end */

c
c     rpi(k) = (sl(k+1)/sl(k))**rk  from setsig  k=1...levm1_
c
      do 9 k=1,levm1_
      rpirec(k) = 1.e0/rpi(k)
9     continue
      do 10 k=1,levs_
      rdel2(k)=0.5e0/del(k)
10    continue
      ind=0
      do 7 ll=1,jcap1_
      n=ll-2
      maxi=jcap1_+1-ll
      do 6 i=1,maxi
      ind=ind+1
      n=n+1
      ndex(ind*2-1) = n
      ndex(ind*2  ) = n
      fact=float(n*(n+1))
      snnp1(ind*2-1) = fact
      snnp1(ind*2  ) = fact
6     continue
7     continue
#ifdef MP
      call spcshfli(snnp1,lnt2_,1,jcap_,lwvdef)
      do n=1,lnt2_
        rdex(n)=ndex(n)
      enddo
      call spcshfli(rdex,lnt2_,1,jcap_,lwvdef)
      do n=1,lnt2_
        ndex(n)=rdex(n)
      enddo
#endif
c    initialize cv, cvt and cvb
      do 40 j=1,LATG2S
      do 40 i=1,LONF2S
        cv (i,j) = 0.e0
        cvt(i,j) = 0.e0
        cvb(i,j) = 100.e0
  40  continue
c
      do 1 i=1,28
1     num(i)=0
      num( 1) = 11
      num( 2) = 11
      num( 3) = 51
      num( 4) = 52
      num( 5) =  0
      num( 6) =  1
      num( 7) =  0
      num( 8) =  1
      num( 9) =  8
      num(10) = 15
      num(11) =  1
      num(12) = 23
cmk  turn off output of initialized sigma
      num(13) =  0
      num(14) = 55
      num(15) =  0
      num(16) = 11
      num(17) = 51
      num(18) =  4
      num(19) =  2
      num(20) =  6
      num(21) = 15
      num(22) = 10
      num(23) =  1
      num(24) =  0
      num(25) =  0
      num(26) =  0
      num(27) =  0
      num(28) =  0
ctemporarily set some cons and nums (setc block data not yet included)
      num(1)=0
      con(1)=0.
      con(3)=0.        ! gsm dfini initialization in hour or not (0.)
c hmhj
      con(4)=1.        ! dtswav in hour for gsm
c hmhj
      con(5)=3.        ! dtlwav in hour for gsm
      con(6)=0.
      con(7)=12.
      con(17)=120.     ! forecast ending hour
#ifdef RSM
      con(11)=400.     ! rsm deltim
      con(12)=21600.   ! rsm nesting period in second
      con(13)=6.0      ! rsm initialization step in hour or not (0.)
      con(14)=1.       ! rsm dtswav in hour
      con(15)=1.       ! rsm dtlwav in hour
      con(16)=0.       ! rsm start forecast period
      con(17)=36.      ! rsm ending forecast period
      con(18)=0.       ! rsm local diffusion (1) or not (0)
      con(19)=0.       ! rsm lateral boundary relaxation 1 blending 0
#endif
      num(31)=1
      num(32)=0
      num(30)=0
      endhour=0.
#ifdef MP
      if( mype.eq.master ) then
#endif
      nf=1
      close (nf)
      open (unit=nf,file='fcst.parm',form='formatted',status='unknown')
      read(nf,namsmf,end=199)
      goto 202
ctemporarily read from original input card if namelist is missing
199   continue
      print *,'eof on namelist read'
#ifdef MP
#ifdef RSM
      call rmpabort
#else
      call mpabort
#endif
#else
      call abort
#endif
202   continue
#ifdef MP
c--------
      ilist(1)=labl
      ilist(2)=ldebug
      ilist(3)=icen
      ilist(4)=igen
      ilist(5)=icen2
      ilist(6)=ienst
      ilist(7)=iensi
      ilist(8)=ncpus
      ilist(9)=intpfix
      ilist(10)=intsfcx 
      ilist(11)=intsig
      ilist(12)=intsfc
      ilist(13)=intpost
      ilist(14)=intrestart
      ilist(15)=intlfm
      ilist(16)=intrsm
      ilist(17)=ndigit
#ifdef COUPLE
      ilist(18)=intcouple
#endif
      rlist(1)=endhour
      rlist(2)=filta
      rlist(3)=runid
      rlist(4)=usrid
#ifdef LFM
      rlist(5)=critfs
      rlist(6)=filtwin
#endif
      endif
c
      call mpbcastr(con,1700)
      call mpbcasti(num,1700)
#ifdef COUPLE
      call mpbcasti(ilist,18)
#else
      call mpbcasti(ilist,17)
#endif
#ifdef LFM
      call mpbcastr(rlist,6)
#else
      call mpbcastr(rlist,4)
#endif
      labl   =ilist(1)
      ldebug =ilist(2)
      icen   =ilist(3)
      igen   =ilist(4)
      icen2  =ilist(5)
      ienst  =ilist(6)
      iensi  =ilist(7)
      ncpus  =ilist(8)
      intpfix=ilist(9)
      intsfcx=ilist(10)
      intsig =ilist(11)
      intsfc =ilist(12)
      intpost =ilist(13)
      intrestart =ilist(14)
      intlfm =ilist(15)
      intrsm =ilist(16)
      ndigit =ilist(17)
      endhour=rlist(1)
      filta  =rlist(2)
      runid  =rlist(3)
      usrid  =rlist(4)
#ifdef LFM
      critfs =rlist(5)
      filtwin=rlist(6)
#endif
c------
#endif
c
c  temporarily reset some cons and nums
c
#ifdef MP
      if( mype.eq.master ) then
#endif
      close(n1)
      open (unit=n1,file='sigit ',form='unformatted',err=999)
      go to 998
  999 continue
      print *,'error opening sigit in getcon'
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
      rewind n1
      read(n1)
      read(n1) fhour
#ifdef MP
      endif
      call mpbcastr(fhour,1)
#endif
      if(fhour.eq.0.) then
        if(num(5).eq.-1) num(5)=2
        if(num(5).eq.-2) num(5)=1
      else
        if(num(5).eq.-1) num(5)=0
        if(num(5).eq.-2) num(5)=0
        if(con(17)-fhour.lt.con(7)) then
          con(7)=con(17)-fhour
        endif
      endif
      if(con(1).le.0.) con(1)=dt80 *80./jcap
      if(num(7).le.0.and.con(7).ne.0.) then
        num(7)=3600.*con(7)/con(1)+0.99
        con(1)=nint(3600.*con(7)/num(7))
      else
        con(7)=num(7)*con(1)/3600.
      endif
      if(num(32).eq.0) num(32)=num(7)
      if(num(1).gt.0) con(6)=num(1)
c
c....
c.... dtswav is interval between short-wave heating calculations
c.... dtlwav is interval between long-wave heating calculations
c....
      dtswav=con(4)
      dtlwav=con(5)
      if(dtswav.gt.float(intpost)) then
        print *,'dtswav.gt.intpost. dtswav=',dtswav,' intpost=',intpost
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
      if(dtlwav.gt.float(intpost)) then
        print *,'dtlwav.gt.intpost. dtlwav=',dtlwav,' intpost=',intpost
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
c
c  check dtswav and dtlwav if they are reasonable
c
      if(mod(int(dtswav*3600.),int(con(1))).ne.0) then
        if(mype.eq.master) 
     1    print *,'dtswav must be a multiple of timestep'
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
      if(mod(24*60,int(dtswav*60.)).ne.0) then
        if(mype.eq.master) 
     1   print *,'24*60 must be must be a multiple of dtswav*60'
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
      if(mod(int(dtlwav*3600.),int(con(1))).ne.0) then
        if(mype.eq.master) 
     1   print *,'dtlwav must be a multiple of timestep'
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
      if(mod(24*60,int(dtlwav*60.)).ne.0) then
        if(mype.eq.master) 
     1   print *,'24*60 must be must be a multiple of dtlwav*60'
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
c
      cowave=0.
      dtwave=0.
c
c>yh  cvmint - maximum conv. cld accumulation time interval in hours
c....          currently hardwired as 3 hours, but may be as an
c....          input variable.
c     cvmint = 3.e0
      cvmint = con(4)
      dtcvav =  min (cvmint,  max (dtswav,dtlwav))
c
#ifdef DBG
      if(mype.eq.master) print 201,(num(i),i=1,28)
201   format(1h0,'num=',28(1x,i2))
      if(mype.eq.master) print *,'con'
      if(mype.eq.master) print *,(con(i),i=1,10)
#endif
c
      if (intpfix.eq.0.) then
        kpfix=0
      else
        kpfix=3600*intpfix/con(1)+0.5
      endif
      ksfcx=3600*intsfcx/con(1)+0.5
      ksig=3600*intsig/con(1)+0.5
      ksfc=3600*intsfc/con(1)+0.5
      kpost=3600*intpost/con(1)+0.5
      krestart=3600*intrestart/con(1)+0.5
      krsm=3600*intrsm/con(1)+0.5
      klfm=3600*intlfm/con(1)+0.5
      kcouple=3600*intcouple/con(1)+0.5
c
      dk=num(9)
      dk=dk*(10.e0)**num(10)
      tk=num(20)
      tk=tk*(10.e0)**num(21)
      if(num(20).eq.0)tk=dk
#ifndef NOPRINT
      if(mype.eq.master) print 105,con(1),filta,dk,tk
105   format(1h ,5x,f5.0,1x,f4.2,1x,e8.2,1x,e8.2)
#endif
#ifdef KEN
      npoint=num(1300)
      if(npoint.lt.0.or.npoint.gt.nptken) then
        print *,'ken points disabled - grid points exceed ',nptken
        npoint=0
      endif
      isave=0
      itnum=0
      if(npoint.ne.0) then
        isave=1
        itnum=1
        isshrt=num(1301)
        ilshrt=num(1302)
        ikfreq=num(1303)
        call kenpre(con,colrad,lonf_,latg2_,n1)
      endif
#endif
      ncpus1=ncpus+1
      ncldb1=ncpus*lonf2_/lonf2_+1
cc
cc    call  cmpind  to set common/comind/ for subs. transi,transo.
#ifndef SCM
      call  cmpind
#ifdef DBG
      if(mype.eq.master) print *,' done cmpind '
#endif
#endif
cc
      call gpvs
#ifdef DBG
      if(mype.eq.master) print *,' done gpvs '
#endif
      call gtdp
#ifdef DBG
      if(mype.eq.master) print *,' done gtdp '
#endif
      call gthe
#ifdef DBG
      if(mype.eq.master) print *,' done gthe '
#endif
      call gtma
#ifdef DBG
      if(mype.eq.master) print *,' done gtma '
#endif
#ifndef SCM
      call gpln2i
#ifdef DBG
      if(mype.eq.master) print *,' done gpln2i '
#endif
      call epsilo(epsi,jcap_)
#ifdef DBG
      if(mype.eq.master) print *,' done epsilo '
#endif
      call ggozrm(epsi)
#ifdef DBG
      if(mype.eq.master) print *,' done ggozrm '
#endif
#endif
#ifdef REDUCE_GRID
      qmaxall=0.0
#endif
#ifdef CLD3
      call cldpar
#endif
#ifdef CLD5
      call cldpar
#endif
#ifdef CLD6
      call cldpar
#endif
#ifdef CLDWSMMPS
      call cldpar
#endif
c
#ifndef SCM
      do lat=1,latg2_
c
        call pln2i(tmpqtt,tmpqvv,colrad,lat)
#ifdef REDUCE_GRID
        do n=1,lnt2_
          qmaxall=max(qmaxall,abs(tmpqtt(n)))
        enddo
#endif
        call gozrin(tmpqtt,tmpqvv,tmpqdd,epsi)
        wcsa=rcs2(lat)/rerth_
#ifdef MP
        call mpsf2p (tmpqtt,lnt2_,tmpq,lnt2p_,1)
        call mpnn2n (tmpq,lnt2p_,qtt(1,lat),lln2p_,1)
        call mpsf2p (tmpqdd,lnt2_,tmpq,lnt2p_,1)
        call mpnn2n (tmpq,lnt2p_,qdd(1,lat),lln2p_,1)
        do n=1,lln2p_
          qww(n,lat)=qtt(n,lat)*wgt(lat)
          qdd(n,lat)=qdd(n,lat)*wcsa
        enddo
        llstr=lwvstr(mype)
        llens=lwvlen(mype)
#else
        do n=1,lnt2_
          qtt(n,lat)=tmpqtt(n)
          qww(n,lat)=tmpqtt(n)*wgt(lat)
          qdd(n,lat)=tmpqdd(n)*wcsa
        enddo
        llstr=0
        llens=jcap1_
#endif
        len=2*jcap1_
        j=len
        do l=0,jcap_
          tmpqvv(2*l+1)=tmpqvv(j+1)
          tmpqvv(2*l+2)=tmpqvv(j+2)
          j=j+len
          len=len-2
        enddo
        do ll=1,llens
          j=lwvdef(llstr+ll)
          l=ll-1
          qvv(2*l+1,lat)=tmpqvv(2*j+1)
          qvv(2*l+2,lat)=tmpqvv(2*j+2)
        enddo
c
      enddo
#ifdef DBG
      if(mype.eq.master) print *,' done qtt qdd qww and qvv '
#endif
#endif
#ifdef REDUCE_GRID
      qttcut=qmaxall/(10.**ndigit)
      if( ndigit.eq.0 ) qttcut=0.0
#ifdef MP
      if( mype.eq.master ) then
#endif
      print *,' reduce grid is on with ',ndigit,' digit accuracy.'
#ifdef MP
      endif
#endif
      do lat=1,latg2_
        call pln2i(tmpqtt,tmpqvv,colrad,lat)
c  row1ns, rowsep depend on  greduceg
        call greduceg(tmpqtt,lnt2_,jcap_,qttcut,
     1                lcapd(lat),lonfd(lat))
#ifdef DBG
#ifdef MP
        if( mype.eq.master ) then
#endif
        print *,' --- lat=',lat,' needs local jcap and lonf as ',
     1                lcapd(lat)-1,lonfd(lat)
#ifdef MP
        endif
#endif
#endif
      enddo
c
#ifdef MP
c --- mpgf2p and mpgp2f depend on preduceg
      call preduceg
#ifdef DBG
      if( mype.eq.master ) print *,' done preduceg '
#endif
#endif
c ------- finish setting of reduce grid --------------
#endif
#ifdef MP
#define HPRIMES gridv
#else
#define HPRIMES hprime
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
        open (unit=nmtnv,file='mtnvar',form='unformatted',err=999)
        read(nmtnv) HPRIMES
#ifdef DBG
#ifndef NOPRINT
        call maxmin(HPRIMES,lonf2_,latg2_,1,1,"mtnvar")
#endif
#endif
        do k=1,mtnvar_
          call row1ns(HPRIMES(1,1,k))
        enddo
#ifdef MP
      endif
      call mpgf2p(gridv,lonf2_,latg2_,hprime,lonf2p_,latg2p_,mtnvar_)
#endif
#ifdef DBG
      if( mype.eq.master ) print *,' done hprime '
#endif
#undef HPRIMES
c
#ifdef MP
#define XLONS grid1
#define XLATS grid2
#else
#define XLONS xlon
#define XLATS xlat
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifndef SCM
        call lonlat(XLONS,XLATS,colrad,lonf_,latg_)
        call row1ns(XLONS)
        call row1ns(XLATS)
#ifdef MP
      endif
      call mpgf2p(grid1,lonf2_,latg2_,xlon,lonf2p_,latg2p_,1)
      call mpgf2p(grid2,lonf2_,latg2_,xlat,lonf2p_,latg2p_,1)
#endif
#ifdef DBG
      if(mype.eq.master) print *,' done xlon xlat '
#endif
#undef XLONS
#undef XLATS
      do 25 j=1,latg2_
        sinlat(j) = cos(colrad(j))
  25  continue
#endif
#ifdef MP
#define SINLABS grid1
#define COSLABS grid2
#else
#define SINLABS sinlab
#define COSLABS coslab
#endif
#ifdef MP
      if( mype.eq.master ) then
#endif
#ifdef SCM
        call sincoslat(SINLABS,COSLABS,xlon,xlat,lonf_,latg_)
#else
        call sincoslat(SINLABS,COSLABS,colrad,lonf_,latg_)
#endif
        call row1ns(SINLABS)
        call row1ns(COSLABS)
#ifdef MP
      endif
      call mpgf2p(grid1,lonf2_,latg2_,sinlab,lonf2p_,latg2p_,1)
      call mpgf2p(grid2,lonf2_,latg2_,coslab,lonf2p_,latg2p_,1)
#endif
#ifdef DBG
      if(mype.eq.master) print *,' done sinlab coslab '
#endif
#undef SINLABS
#undef COSLABS
c
#ifndef SCM
      call gftlonf
#ifdef DBG
      if(mype.eq.master) print *,' done gftlonf '
#endif
      call grddf
#ifdef DBG
      if(mype.eq.master) print *,' done grddf '
#endif
#endif
      call grdkt
#ifdef DBG
      if(mype.eq.master) print *,' done grdkt '
#endif
      call gradini
#ifdef DBG
      if(mype.eq.master) print *,' done gradini '
#endif
c
c  filtwin in the unit of hour
c
#ifdef LFM
      sechr=60.*60.
      klenp=nint(filtwin*sechr/con(1))+1
      nlenp=(klenp-1)/2+1
      critfl=9.0e10
      clancz=1.0
      tinc=con(1)
      critsc=critfs*sechr
#ifndef NOPRINT
      if( mype.eq.master )
     1write(6,*) 'klenp,nlenp,tinc,critfs,critsc,critfl,clancz=',
     1            klenp,nlenp,tinc,critfs,critsc,critfl,clancz
#endif
      call filtcof(nlenp,tinc,critsc,critfl,clancz,weix)
      fwindow=filtwin
#endif
c
c  expand weights
c
#ifdef LFM
      weight(nlenp)=weix(1)
      do k=2,nlenp
        weight(nlenp+k-1)=weix(k)
        weight(nlenp-k+1)=weix(k)
      enddo
#ifndef NOPRINT
      do k=1,klenp
        if(mype.eq.master) print *,k,' weight=',weight(k)
      enddo
#endif
#endif
c
      return
      end
