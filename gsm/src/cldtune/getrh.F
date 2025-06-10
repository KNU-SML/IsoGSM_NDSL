#define JCAP jcap
#define LATG latg
#define LEVS levs
#define LONF lonf
#include <define.h>
      subroutine getrh(fnsig,vvcldin,fhour,idate,rhmdl,ieee)
#include <paramodel.h>
      save
c
c-    *****************************************************************
c-      .h30k09(getrh28v) adapted to new cloud code..
c-        vvel limits adjusted, ptop is function of latitude...but
c-        rhmaxx for low clouds still encompasses clds abuv llyr
c-    *****************************************************************
      parameter (idim=lonf_,jdim=latg_,
     1           mwave=jcap_,jromb=0,kdim=levs_,mcld=3,
     2           kdimp=kdim+1,ijdim=idim*jdim,
     3           mwvd=(mwave+1)*(mwave+2)/2*(1-jromb)+
     4                (mwave+1)*(mwave+1)*jromb,
     5           mwvd2=mwvd*2  )
c
      character*80 fnsig,asgnstr
c
#ifdef ASSIGN
      dimension gout(idim,jdim)
#else
      real*4 sigin(mwvd2)
      real*4 sfhour,ssi(kdimp),ssl(kdim)
      integer*4 idates(4)
      real*4 gout(idim,jdim)
#endif
c
      dimension lab(4)
      dimension idate(4),fac(kdim)
      dimension psfc(mwvd2)
      dimension temp(mwvd2,kdim),q(mwvd2,kdim)
      dimension div(mwvd2,kdim)
      dimension vor(mwvd2,kdim)
cwinds
      dimension fer(mwvd),qimer(mwvd),sinclt(jdim)
      dimension gaul(jdim)
cwinds
      dimension tgrid(idim,jdim,kdim),qgrid(idim,jdim)
      dimension pstar(idim,jdim),rh(idim,jdim,kdim)
      dimension tt(idim,kdim),rhcldy(idim,mcld),rhrh(idim,kdim)
cwinds
      dimension urot(ijdim),udiv(ijdim)
      dimension ggrid(ijdim)
cwinds
cvvel
      dimension frx(mwvd),fry(mwvd)
      dimension del(kdim)
      dimension dphi(idim,jdim),dlam(idim,jdim)
      dimension vvel(idim,jdim,kdim),divg(idim,kdim)
      dimension ddiivv(idim,jdim,kdim),
     1          uw(idim,jdim,kdim),vw(idim,jdim,kdim)
#ifdef DBG
      dimension uwd(idim,jdim,kdim),vwd(idim,jdim,kdim),
     1          psi(idim,jdim,kdim),chi(idim,jdim,kdim)
#endif
      dimension uwind(idim,kdim),vwind(idim,kdim)
      dimension vv(idim,kdim)
cvvel
      dimension rhmdl (idim,jdim,mcld)
      dimension si(kdimp),sl(kdim)
      dimension yave(jdim),globl(9),rlat(jdim),wgt(jdim)
      dimension jpts(jdim),wgta(jdim)
      dimension colrad(jdim)
      dimension colr(jdim)
      dimension rlatd(jdim)
      character*1 ptitle(20)
      character*12 lbl,ltitle(mcld,3)
      dimension nundg(mcld),nok(mcld),nvvel(mcld)
      data (ptitle(k),k=1,20) /' ',' ',' ',' ',' ',' ',' ',' ',
     1         ' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '/
      data ltitle/
     1  'rtnef lowcld','rtnef midcld','rtnef hi cld',
     1  'mrf hi  cld ','mrf mid cld ','mrf low cld ',
     1  'mrf low -rh-','mrf mid -rh-','mrf hi -rh- '/
c     equivalence(psfc,temp)
      equivalence(lbl,ptitle(1))
c
c----    begin here .........
c          this code accesses global sigma data from sganl...
c          for 4 cycles..0,6,12,18 z..(numcy=4)
c     read (5,7) numcy
c   7 format(i3)
c
      numcy=1
c
      khalf = jdim/2
      xhalf = khalf
      call glatss(khalf,colrad)
      do 152 l=1,khalf
        colr(l) = colrad(l)
        rlatd(l) = 3.14159265/2. - colrad(l)
        rlatd(jdim+1-l) = rlatd(l)
  152 continue
      ii = idim
      jj = jdim
      kk = kdim
      iromb = jromb
      kmaxp=kdimp
      mwve=(mwave+1)*(mwave+2)/2*(1-iromb)+(mwave+1)*(mwave+1)*iromb
      mwve2=mwve*2
c
c....   program sphert converts sig coeff to equidistant 2 deg grid
c                  if idir=-101,,to gaussian grid if idir=-1
c     for dplam idir,mlt,fac = -1,-1,float(m)/a
c     for dpphi idir,mlt,fac = -2, 1, 1 /a
c
      idir=-1
cwinds
      idirx=-2
      er    =6371.0e3
      eriv  =1./er
      ersq=er*er
      pi=4.*atan(1.)
      rad=pi/180.
      call gaulat(gaul,jdim)
      do 167 j=1,jdim
      gaul(j)=90.-gaul(j)
  167 continue
      do 69 j=1,jdim
      sinclt(j)=cos(gaul(j)*rad)
   69 continue
      mend1=mwave+1
      nm =0
      do 33 mm=1,mend1
      m=mm-1
      nend1=mend1
      if(iromb.eq.1) nend1=mm+mend1-1
      do 43 nn=mm,nend1
      n=nn-1
      x1=-float(n*(n+1))
      x4=0.
      if(n.ne.0) x4=ersq/x1
      nm=nm+1
      fer   (nm)=x4*eriv
      fry   (nm)= eriv
      qimer (nm)=x4*eriv*float(m)
      frx (nm)=eriv*float(m)
   43 continue
   33 continue
c
cwinds
c
      do 555 k=1,kdim
      do 555 j=1,jdim
      do 555 i=1,idim
         vvel(i,j,k) = 0.
  555 continue
      do 55 k=1,kdim
       fac(k) = 1.
   55 continue
c
c  setup in order to compute rh as in mrf model
c
      rd=287.05
      rv=461.5
      eps=rd/rv
      epsm1=rd/rv-1.
      call gpvs
c
      irhm = 16
      rewind irhm
c
        iunit= 11
c
#ifdef ASSIGN
        ieee=1
        if(ieee.eq.1) then
          write(asgnstr,'(22hassign -nieee -ff77 u:,i2,)') iunit
        elseif(ieee.eq.0) then
          write(asgnstr,'(29hassign -fcos -cascii -nibm u:,i2,)') iunit
        else
          write(asgnstr,'(9hassign u:,i2,)') iunit
        endif
        call assign('assign -R')
        call assign(asgnstr)
#endif
c
        open(unit=iunit,file=fnsig,form='unformatted',err=788)
        go to 787
  788   continue
        print *,'error in opening file ',fnsig
        call abort
  787   continue
#ifdef ASSIGN
        write(6,*) 'asgnstr=',asgnstr
#endif
        write(6,*) 'file ',fnsig(1:80),' opened. unit=',iunit
c
        rewind iunit
        read (iunit,err=998,end=999) lab
        write(6,3001) lab
        if(ieee.eq.2) then
          read (iunit,err=998,end=999) fhour,(idate(i),i=1,4),
     1            (si(k),k=1,kmaxp),(sl(k),k=1,kk)
        else
#ifdef ASSIGN
          read (iunit,err=998,end=999) fhour,(idate(i),i=1,4),
     1             (si(k),k=1,kmaxp),(sl(k),k=1,kk)
#else
          read (iunit,err=998,end=999) sfhour,(idates(i),i=1,4),
     1             (ssi(k),k=1,kmaxp),(ssl(k),k=1,kk)
          fhour=sfhour
          do i=1,4
           idate(i)=idates(i)
          enddo
          do k=1,kk
           sl(k)=ssl(k)
          enddo
          do k=1,kmaxp
           si(k)=ssi(k)
          enddo
#endif
        endif
        write(6,3002) idate,fhour
        iy=idate(4)+1900
        im=idate(2)
        id=idate(3)
        ih=idate(1)
        ift=fhour
c
c ground hgts
c
#ifdef ASSIGN
        read (iunit,err=998,end=999) zsfcx
#else
        read (iunit,err=998,end=999) sigin
#endif
c
c natural log sfc pressure
c
        if(ieee.eq.2) then
          read (iunit,err=998,end=999) (psfc(nm),nm=1,mwve2)
        else
#ifdef ASSIGN
          read (iunit,err=998,end=999) (psfc(nm),nm=1,mwve2)
#else
          read (iunit,err=998,end=999) sigin
          do i=1,mwve2
            psfc(i)=sigin(i)
          enddo
#endif
        endif
c
c natural log sfc pressure
c
        call sphert(idir,tgrid,psfc,0,fac,ii,jj,mwave,iromb)
c
c get units of sfc pressure in cb
c
        do 67 j=1,jdim
         do 67 i=1,idim
          pstar(i,j)=exp(tgrid(i,j,1))
   67   continue
        do 29 k=1,kk
        if(ieee.eq.2) then
          read (iunit,err=998,end=999) (temp(nm,k),nm=1,mwve2)
        else
#ifdef ASSIGN
          read (iunit,err=998,end=999) (temp(nm,k),nm=1,mwve2)
#else
          read (iunit,err=998,end=999) sigin
          do i=1,mwve2
            temp(i,k)=sigin(i)
          enddo
#endif
        endif
   29   continue
c
c get div and vor
c
        do 30 k=1,kk
        if(ieee.eq.2) then
          read (iunit,err=998,end=999) (div(nm,k),nm=1,mwve2)
          read (iunit,err=998,end=999) (vor(nm,k),nm=1,mwve2)
        else
#ifdef ASSIGN
          read (iunit,err=998,end=999) (div(nm,k),nm=1,mwve2)
          read (iunit,err=998,end=999) (vor(nm,k),nm=1,mwve2)
#else
          read (iunit,err=998,end=999) sigin
        do i=1,mwve2
          div(i,k)=sigin(i)
        enddo
          read (iunit,err=998,end=999) sigin
          do i=1,mwve2
            vor(i,k)=sigin(i)
          enddo
#endif
        endif
   30   continue
        do 31 k=1,kk
        if(ieee.eq.2) then
          read (iunit,err=998,end=999) (q(nm,k),nm=1,mwve2)
        else
#ifdef ASSIGN
          read (iunit,err=998,end=999) (q(nm,k),nm=1,mwve2)
#else
          read (iunit,err=998,end=999) sigin
          do i=1,mwve2
            q(i,k)=sigin(i)
          enddo
#endif
        endif
   31   continue
c
c get dlnps/dlam, dlnps/dphi
c
        call sphert(idir,dlam,psfc,-1,frx,ii,jj,mwave,iromb)
        call sphert(idirx,dphi,psfc,1,fry,ii,jj,mwave,iromb)
c
        do j=1,idim
          do i=1,idim
            dphi(i,j)=-dphi(i,j)
          enddo
        enddo
c
        do 393 k=1,kk
         del(k) = si(k) - si(k+1)
#ifdef DBG
         print 405,si(k),sl(k),del(k)
#endif
  393   continue
  405 format(1h ,'si,sl,del=',3(e15.5,5x))
c
c wind level loop starts
c
        do 845 k=1,kk
c
c  u-wind
c
          call sphert(idirx,urot,vor(1,k),1,fer,ii,jj,mwave,iromb)
          call sphert(idir ,udiv,div(1,k),-1,qimer,ii,jj,mwave,iromb)
          call sphert(idir,ddiivv(1,1,k),div(1,k),0,fac,ii,jj,
     1                mwave,iromb)
#ifdef DBG
          call sphert(idir,psi(1,1,k),vor(1,k),1,fer,ii,jj,mwave,iromb)
          call sphert(idir,chi(1,1,k),div(1,k),1,fer,ii,jj,mwave,iromb)
#endif
          do 501 j=1,jdim
           ix2=(j-1)*idim
           xx=0.
           if(sinclt(j).ne.0.) xx=1./sinclt(j)
           do 501 i=1,idim
            urot(ix2+i)= urot(ix2+i)*xx
            udiv(ix2+i)= udiv(ix2+i)*xx
  501     continue
          do 9043 i=1,ijdim
           ggrid(i) = urot(i) + udiv(i)
 9043     continue
          it = 1
          do 8043 j=1,jdim
           do 8043 i=1,idim
            uw(i,j,k) = ggrid(it)
#ifdef DBG
            uwd(i,j,k) = udiv(it)
#endif
            it = it + 1
 8043     continue
#ifdef DBG
          xmax = ggrid(1)
          xmin = ggrid(1)
          do 9033 i=1,ijdim
           if(ggrid(i).gt.xmax) xmax =ggrid(i)
           if(ggrid(i).lt.xmin) xmin =ggrid(i)
#endif
 9033     continue
#ifdef DBG
          print 911,k,xmax,xmin
#endif
c
c v-wind
c
          call sphert(idir ,urot,vor(1,k),-1,qimer,ii,jj,mwave,iromb)
          call sphert(idirx,udiv,div(1,k),1,fer,ii,jj,mwave,iromb)
          do 601 j=1,jdim
           ix2=(j-1)*idim
           xx=0.
           if(sinclt(j).ne.0.) xx=1./sinclt(j)
           do 601 i=1,idim
            urot(ix2+i)= urot(ix2+i)*xx
            udiv(ix2+i)=-udiv(ix2+i)*xx
  601     continue
          do 9053 i=1,ijdim
           ggrid(i) = urot(i) + udiv(i)
 9053     continue
          it = 1
          do 8053 j=1,jdim
           do 8053 i=1,idim
            vw(i,j,k) = ggrid(it)
#ifdef DBG
            vwd(i,j,k) = udiv(it)
#endif
            it = it + 1
 8053     continue
#ifdef DBG
          xmax = ggrid(1)
          xmin = ggrid(1)
          do 9034 i=1,ijdim
           if(ggrid(i).gt.xmax) xmax =ggrid(i)
           if(ggrid(i).lt.xmin) xmin =ggrid(i)
#endif
 9034     continue
#ifdef DBG
          print 912,k,xmax,xmin
#endif
c
  845   continue
c
c
#ifdef DBG
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
             gout(i,j)=uw(i,j,k)
            enddo
          enddo
          write(80) gout
        enddo
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
             gout(i,j)=vw(i,j,k)
            enddo
          enddo
          write(80) gout
        enddo
#endif
c
#ifdef DBG
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
             gout(i,j)=uwd(i,j,k)
            enddo
          enddo
          write(80) gout
        enddo
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
             gout(i,j)=vwd(i,j,k)
            enddo
          enddo
          write(80) gout
        enddo
#endif
c
#ifdef DBG
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
             gout(i,j)=ddiivv(i,j,k)
            enddo
          enddo
          write(80) gout
        enddo
#endif
c
#ifdef DBG
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
             gout(i,j)=psi(i,j,k)
            enddo
          enddo
          write(80) gout
        enddo
#endif
c
#ifdef DBG
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
             gout(i,j)=chi(i,j,k)
            enddo
          enddo
          write(80) gout
        enddo
#endif
c
cwinds
cvvel
c
#ifdef DBG
          do j=1,jdim
            do i=1,idim
             gout(i,j)=dlam(i,j)
            enddo
          enddo
          write(80) gout
#endif
c
#ifdef DBG
          do j=1,jdim
            do i=1,idim
             gout(i,j)=dphi(i,j)
            enddo
          enddo
          write(80) gout
#endif
c
        do 945 j=1,jdim
         do 817 k=1,kdim
         do 817 ik=1,idim
           uwind(ik,k)=uw(ik,j,k)
           vwind(ik,k)=vw(ik,j,k)
           divg (ik,k)=ddiivv(ik,j,k)
  817    continue
         call omegas(dphi(1,j),dlam(1,j),uwind,vwind,
     1                divg,del,vv,pstar(1,j),sl)
         do 818 k=1,kdim
          do 818 i=1,idim
           vvel(i,j,k) = vv(i,k)
  818     continue
  945   continue
c
#ifdef DBG
        print *,'vvel max/min'
        do k=1,kdim
          xmax=vvel(1,1,k)
          xmin=vvel(1,1,k)
          do j=1,jdim
            do i=1,idim
              if(xmax.lt.vvel(i,j,k)) xmax=vvel(i,j,k)
              if(xmin.gt.vvel(i,j,k)) xmin=vvel(i,j,k)
            enddo
          enddo
          print *,'k=',k,' max=',xmax,' min=',xmin
        enddo
#endif
c
#ifdef DBG
        do k=1,kdim
          do j=1,jdim
            do i=1,idim
             gout(i,j)=vvel(i,j,k)
            enddo
          enddo
          write(80) gout
        enddo
#endif
c
cvvel
c
        do 32 k=1,kk
          call sphert(idir,tgrid(1,1,k),temp(1,k),0,fac,ii,jj,
     1                mwave,iromb)
          call sphert(idir,qgrid,q(1,k),0,fac,ii,jj,mwave,iromb)
#ifdef DBG
          xmax = tgrid(1,1,k)
          xmin = tgrid(1,1,k)
          do 924 j=1,jdim
          do 924 i=1,idim
           if(tgrid(i,j,k).gt.xmax) xmax =tgrid(i,j,k)
           if(tgrid(i,j,k).lt.xmin) xmin =tgrid(i,j,k)
#endif
  924     continue
#ifdef DBG
          print 913,k,xmax,xmin
          xmax = qgrid(1,1)
          xmin = qgrid(1,1)
#endif
          do 926 j=1,jdim
          do 926 i=1,idim
           if(qgrid(i,j).lt.0.) qgrid(i,j)=0.
#ifdef DBG
           if(qgrid(i,j).gt.xmax) xmax =qgrid(i,j)
           if(qgrid(i,j).lt.xmin) xmin =qgrid(i,j)
#endif
  926     continue
#ifdef DBG
          print 914,k,xmax,xmin
#endif
c
c convert tv to t,calculate-qs-,-rh-as in lrgscl
c
          do 27 j=1,jdim
           do 27 i=1,idim
           tgrid(i,j,k)=tgrid(i,j,k)/(1. + 0.608 *qgrid(i,j))
           es=fpvs(tgrid(i,j,k))
           qs=eps*es/(sl(k)*pstar(i,j)+epsm1*es)
           rh(i,j,k)=qgrid(i,j)/qs
   27     continue
c
#ifdef DBG
          xmax = rh(1,1,k)
          xmin = rh(1,1,k)
          do 928 j=1,jdim
          do 928 i=1,idim
           if(rh(i,j,k).gt.xmax) xmax = rh(i,j,k)
           if(rh(i,j,k).lt.xmin) xmin = rh(i,j,k)
#endif
  928     continue
#ifdef DBG
          print 915,k,xmax,xmin
#endif
          do 8 j=1,jdim
           do 8 i=1,idim
            if(rh(i,j,k).gt.1.) rh(i,j,k) = 1.
            if(rh(i,j,k).lt.0.) rh(i,j,k) = 0.
    8     continue
c
   32   continue
c
c
c
        do k=1,mcld
          nundg(k)=0
          nok  (k)=0
          nvvel(k)=0
        enddo
        do 34 j=1,jdim
         do 40 k=1,kdim
           do 40 i=1,idim
             tt(i,k)=tgrid(i,j,k)
             rhrh(i,k)=rh(i,j,k)
             vv(i,k)=vvel(i,j,k)
   40    continue
         call rhmaxx(pstar(1,j),rhrh,tt,vv,si,sl,vvcldin,
     1               rhcldy,rlatd(j),nundg,nok,nvvel)
         do 41 kc=1,mcld
          do 41 i=1,idim
           rhmdl(i,j,kc) = rhcldy(i,kc)
   41    continue
   34   continue
c
       print *,'global statistics'
         mmcld=1
#ifdef DBG
         mmcld=mcld
#endif
         do k=1,mmcld
       print *,'k=',k,' totl=',idim*jdim,
     1         ' ctotal=',nundg(k)+nok(k)+nvvel(k),
     2         ' undg=',nundg(k),' ok=',nok(k),' nvvel=',nvvel(k)
         enddo
c
 1000 continue
c
#ifdef DBG
      do 357 kc=1,mcld
       xmax = rhmdl(1,1,kc)
       xmin = rhmdl(1,1,kc)
       do 359 j=1,jdim
        do 359 i=1,idim
         if(rhmdl(i,j,kc).gt.xmax) xmax =rhmdl(i,j,kc)
         if(rhmdl(i,j,kc).lt.xmin) xmin =rhmdl(i,j,kc)
#endif
  359  continue
#ifdef DBG
       print 923,kc,xmax,xmin
#endif
  357 continue
c
c write rh max in cld domains as h,m,l
c
 1001 continue
c     write (irhm) fhour,idate
      do 2901 k=1,mcld
       kmc = mcld+1-k
c
c      write(irhm) ((rhmdl(i,j,kmc),i=1,idim),j=1,jdim)
c
       lbl = ltitle(kmc,3)
       call avlat(yave,rhmdl(1,1,kmc),idim,jdim,1,kmc,jpts)
       call globlm(yave,globl,rlat,wgt,wgta,
     1              1,jdim,colr,ptitle,fhour,jpts)
 2901 continue
c
      return
c
  998  continue
       write(6,*) 'sigma file read error'
       call abort
  999  continue
       write(6,*) 'sigma file end of file'
       call abort
  911 format(1h ,' -u-  k = ',i4,' max = ',e15.5,' min =',e15.5)
  912 format(1h ,' -v-  k = ',i4,' max = ',e15.5,' min =',e15.5)
  913 format(1h ,' -t-  k = ',i4,' max = ',e15.5,' min =',e15.5)
  914 format(1h ,' -q-  k = ',i4,' max = ',e15.5,' min =',e15.5)
  915 format(1h ,' -rh- k = ',i4,' max = ',e15.5,' min =',e15.5)
  923 format(1h ,' -rhmdl- k= ',i4,' max = ',e15.5,' min =',e15.5)
 3001 format(5x,'input label= ',a8,2x,2z16,2x,a8)
 3002 format(5x,'input data dates time and forecast hour',/,10x,
     1       4i15,f7.1)
 3006 format(5x,'... read error.. time and forecast hour',/,10x,
     1       4i15,f7.1)
 3007 format(5x,'... endoffile .. time and forecast hour',/,10x,
     1       4i15,f7.1)
      end
