#define IMN imn
#define JCAP jcap
#define JMN jmn
#define LATG latg
#define LONF lonf
#define GSMTNRES gsmtnres
#define MTNVAR mtnvar
      subroutine tersub(im,jm,nm,nr,nw)
#include <paramodel.h>
      parameter(imn=imn_,jmn=jmn_,numvar=mtnvar_)
c
      integer zavg(imn,jmn),zvar(imn,jmn),zmax(imn,jmn),zslm(imn,jmn)
c
      parameter(nnw=(jcap_+1)*(jcap_+2))
c
      real cosclt(latg_),wgtclt(latg_),rclt(latg_),xlat(latg_)
      real slm(lonf_,latg_),oro(lonf_,latg_),var(lonf_,latg_)
      real ors(nnw),orf(lonf_,latg_)
      real orsmth(nnw)
c
      real var4(lonf_,latg_),oa(lonf_,latg_,4),ol(lonf_,latg_,4)
      dimension ist(lonf_,latg_),ien(lonf_,latg_),
     1          jst(latg_),jen(latg_)
      dimension iwork(lonf_,latg_,4),
     1     work(lonf_,latg_),work1(lonf_,latg_),
     1     work2(lonf_,latg_),work3(lonf_,latg_),work4(lonf_,latg_),
     2     work5(lonf_,latg_),work6(lonf_,latg_),
     3     glat(jmn_)
      dimension flat(lonf_,latg_),flon(lonf_,latg_)
      dimension hprime(lonf_,latg_,10)
      equivalence (var (1,1),hprime(1,1,1))
      equivalence (var4(1,1),hprime(1,1,2))
      equivalence (oa(1,1,1),hprime(1,1,3))
      equivalence (oa(1,1,2),hprime(1,1,4))
      equivalence (oa(1,1,3),hprime(1,1,5))
      equivalence (oa(1,1,4),hprime(1,1,6))
      equivalence (ol(1,1,1),hprime(1,1,7))
      equivalence (ol(1,1,2),hprime(1,1,8))
      equivalence (ol(1,1,3),hprime(1,1,9))
      equivalence (ol(1,1,4),hprime(1,1,10))
c
      real*4 orfviw(lonf_,latg_)
c
      logical lsmth,lenh
c
      data lsmth/.true./,rk/0.7e-6/,lenh/.false./,enhc/1.2/
c
c     restricted data 
      integer numi(latg_)
      data  numi/latg_*lonf_/
c     data numi  /
c    .   30,  30,  30,  40,  48,  56,  60,  72,  72,  80,  90,  90,
c    .   96, 110, 110, 120, 120, 128, 144, 144, 144, 144, 154, 160,
c    .  160, 168, 168, 180, 180, 180, 180, 180, 180, 192, 192, 192,
c    .  192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192 /
c     do i = latg_/2+1,latg_
c       numi(i) = numi(latg_+1-i)
c     enddo
c
c  set constants and zero fields
      pi=4.*atan(1.)
      degrad = 180./pi
c
c- open(unit=11,form='formatted',err=900) ! average 
c- open(unit=12,form='formatted',err=900) ! variance 
c- open(unit=13,form='formatted',err=900) ! maximum 
c- open(unit=14,form='formatted',err=900) ! sea-land-lake-mask 
c
      read(11,11) zavg
      read(12,11) zvar
      read(13,11) zmax
      read(14,12) zslm
11    format(20i4)
12    format(80i1)
c
#ifdef SCM
c
#define GXLON _gxlon_
#define GYLAT _gylat_
#define GLSMSK _lsmsk_
c                           ... observation position
      xlon = GXLON
      ylat = GYLAT
      lsmask = GLSMSK
c
      delxn = 360./imn      ! mountain data resolution
      ipos = int(xlon/delxn) + 1
      jpos = int((ylat+90.0)/delxn) + 1
      print *, 'POS ',imn,jmn,ipos,jpos,delxn
      istr = ipos-4
      iend = ipos+4
      jstr = jpos-4
      jend = jpos+4
c
      xnsum = 0.0
      xland = 0.0
      xwatr = 0.0
      xl1 = 0.0
      xs1 = 0.0
      xw1 = 0.0
      xw2 = 0.0
      xv2 = 0.0
      xw4 = 0.0
      do im = istr, iend
      do jm = jstr, jend
        xland = xland + float(zslm(im,jm))
        xwatr = xwatr + float(1-zslm(im,jm))
        xnsum = xnsum + 1.
        height = float(zavg(im,jm))
        if(height.lt.-990.) height = 0.0
        xl1 = xl1 + height * float(zslm(im,jm))
        xs1 = xs1 + height * float(1-zslm(im,jm))
        xw1 = xw1 + height
        xw2 = xw2 + height ** 2
        xv2 = xv2 + float(zvar(im,jm)) ** 2
        print 1967, im,jm,xland,xwatr,xnsum,height,xl1,xs1,xw1,xw2,xv2
1967    format(2i5,9e13.5)
      enddo
      enddo
      if(xnsum.gt.1.) then
         slm1p = float(nint(xland/xnsum))
         if (float(lsmask).eq.slm1p) then
           if(slm1p.ne.0.) then
              oro1p = xl1 / xland
           else
              oro1p = xs1 / xwatr
           endif
         else
           print *, 'Land-Sea mask is not consistent with "define.h" ! '
           call abort
         endif
         var1p=sqrt(max((xv2+xw2)/xnsum-(xw1/xnsum)**2,0.))
         print 1966, xnsum,slm1p,oro1p,var1p
1966     format(4e13.5)
      endif
      do i = 1, lonf_
      do j = 1, latg_
        oro(i,j) = oro1p
        var(i,j) = var1p
      enddo
      enddo
      write(53) var
      write(62) oro
#else

      call gausslat(jm,cosclt,wgtclt)
      do j=1,jm/2
        rclt(j)=acos(cosclt(j))
      enddo
      do j = 1,jm/2
         phi = rclt(j) * degrad
         xlat(j) = 90. - phi
         xlat(jm-j+1) =  phi - 90.
      enddo
c
      dlon=360./float(lonf_)
      do j=1,latg_
        do i=1,lonf_
          flat(i,j)=xlat(j)
          flon(i,j)=dlon*float(i-1)
        enddo
      enddo
c
c     compute mountain data : oro slm var oc
c
      call makemt(zavg,zvar,zslm,oro,slm,var,var4,glat,
     1            ist,ien,jst,jen,im,jm,imn,jmn,xlat,numi)
c
c     compute mountain data : oa ol
c
      call makeoa(zavg,zmax,var,glat,oa,ol,iwork,
     1            work1,work2,work3,work4,
     2            work5,work6,
     3            ist,ien,jst,jen,im,jm,imn,jmn,xlat,numi)
c
c  enhancement
c
      if(lenh) then
        do j=1,jm
          do i=1,im
            oro(i,j)=oro(i,j)+var(i,j)*enhc
          enddo
        enddo
      endif
c
c  remove isolated points
c
      do j=2,jm-1
        jn=j-1
        js=j+1
        do i=1,im
          iw=mod(i+im-2,im)+1
          ie=mod(i,im)+1
          slma=slm(iw,jn)+slm(i,jn)+slm(ie,jn)+
     &         slm(iw,j )          +slm(ie,j )+
     &         slm(iw,js)+slm(i,js)+slm(ie,js)
          if(slm(i,j).eq.0..and.slma.eq.8.) then
            slm(i,j)=1.
            oroa=(oro(iw,jn)+oro(i,jn)+oro(ie,jn)+
     &            oro(iw,j )          +oro(ie,j )+
     &            oro(iw,js)+oro(i,js)+oro(ie,js))/8.
            vara=(var(iw,jn)+var(i,jn)+var(ie,jn)+
     &            var(iw,j )          +var(ie,j )+
     &            var(iw,js)+var(i,js)+var(ie,js))/8.
            print '("sea ",2f8.0," modified to land",2f8.0," at ",2i8)',
     &       oro(i,j),var(i,j),oroa,vara,i,j
            oro(i,j)=oroa
            var(i,j)=vara
          elseif(slm(i,j).eq.1..and.slma.eq.0.) then
            slm(i,j)=0.
            oroa=(oro(iw,jn)+oro(i,jn)+oro(ie,jn)+
     &            oro(iw,j )          +oro(ie,j )+
     &            oro(iw,js)+oro(i,js)+oro(ie,js))/8.
            vara=(var(iw,jn)+var(i,jn)+var(ie,jn)+
     &            var(iw,j )          +var(ie,j )+
     &            var(iw,js)+var(i,js)+var(ie,js))/8.
            print '("land",2f8.0," modified to sea ",2f8.0," at ",2i8)',
     &       oro(i,j),var(i,j),oroa,vara,i,j
            oro(i,j)=oroa
            var(i,j)=vara
          endif
        enddo
      enddo
c
c     zero over ocean
c
      do i = 1,im
        do j = 1,jm
          if(slm(i,j).eq.0.) then
            var4(i,j) = 0.
            oa(i,j,1) = 0.
            oa(i,j,2) = 0.
            oa(i,j,3) = 0.
            oa(i,j,4) = 0.
            ol(i,j,1) = 0.
            ol(i,j,2) = 0.
            ol(i,j,3) = 0.
            ol(i,j,4) = 0.
          endif
       enddo
      enddo
c  spectrally truncate orography
      call sphert(1,oro,ors,0,0.,im,jm,nm,nr)
c
c spectral smoothing 
c
      if(lsmth) then
        l=1
        do m=1,jcap_+1
          do n=m,jcap_+1
c
c  hoskins(1990)
c
c           nnx=(n-1)*(n-1)*n*n
c           fac=exp(-rk*float(nnx))
c
c  navarra(1994)
c
c           fac=1.-float(n-1)/float(jcap_+1)
c
c  empirical (kanamitsu, 1998)
c
c           beta=1.
c           alfa=0.6
c           fac=1.-alfa*(float(n-1)/float(jcap_+1))**beta
c
c  empirical (ebisuzaki, 1998)
c
            width=0.6
            nmax=sqrt(float((jcap_)*(jcap_+1)))
            nmin=(1.-width)*nmax
            nn=sqrt(float((n-1)*n))
            if(nn.lt.nmin) then
              fac=1.
            else
              fac=(float(nmax)-float(nn))/float(nmax-nmin)
            endif
c
c  empirical (navarra, 1994)
c
c           alfa=32.
c           beta=8.
c           fac=exp(-alfa*(n/jcap_)**(2*beta))
c
            orsmth(l  )=fac*ors(l  )
            orsmth(l+1)=fac*ors(l+1)
            l=l+2
          enddo
        enddo
      endif
c
      call sphert(-1,orf,ors,0,0.,im,jm,nm,nr)
c  output fields
      write(51) slm
      write(51) flat
      write(51) flon
      do j=1,jm
        do i=1,im
          orfviw(i,j)=slm(i,j)
        enddo
      enddo
      write(63) orfviw
      call nntprt(slm,im,jm,1.)
      write(52) orf
      call nntprt(orf,im,jm,0.001)
      if(numvar.eq.1) then
         write(53) var
      else
         write(53) hprime
      endif
      do ivar = 1,10
      do i = 1,im
      do j = 1,jm
          work(i,j) = hprime(i,j,ivar)
          orfviw(i,j)=work(i,j)
      enddo
      enddo
      call nntprt(work,im,jm,1.)
      write(63) orfviw
      enddo
      write(54) ors
      write(55) oro
      call nntprt(oro,im,jm,0.001)
c
      write(61) orsmth
c
      do j=1,jm
        do i=1,im
          orfviw(i,j)=oro(i,j)
        enddo
      enddo
      write(63) orfviw
      do j=1,jm
        do i=1,im
          orfviw(i,j)=orf(i,j)
        enddo
      enddo
      write(63) orfviw
c
      call sphert(-1,orf,orsmth,0,0.,im,jm,nm,nr)
      write(62) orf
      do j=1,jm
        do i=1,im
          orfviw(i,j)=orf(i,j)
        enddo
      enddo
      write(63) orfviw
c
#endif
      return
      end
