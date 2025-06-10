#define IGRD12 igrd12
#define IGRD igrd
#define JGRD12 jgrd12
#define JGRD jgrd
#define LATG2 latg2
#define LATG latg
#define LONF2 lonf2
#define LONF lonf
#define PI pi
      subroutine nfdr_ll2rsm(rdat,imin,jmin,inttyp,rout,
     1                  lmask,rslmsk,slmask,
     2                  proji,orienti,truthi,cotrui,
     3                  delxi,delyi,rlat1i,rlat2i,rlon1i,rlon2i)
cfpp$ noconcur r
c
c  interpolation from one projection to rsm model projection
c
c  inttyp ... type of interpolation.  0 = bilinear
c                                     1 = take the closest point value
c                                     2 = take the predominant type
c
#include <paramodel.h>
#ifdef RSM
#define ILOTS igrd12_
#define JLOTS jgrd12_
#else
#define ILOTS lonf2_
#define JLOTS latg2_
#endif
      parameter(idim=ILOTS/2,jdim=JLOTS*2)
      parameter(ijdim=idim*jdim)
      dimension rdat(imin*jmin)
c
      parameter(io=184,jo=116)
      dimension rout(io,jo)
c
      real, allocatable :: glon(:),glat(:),glao(:),tmp(:,:)
c
      dimension gaul(3000)
c
c fire danger model grid
c
      real projo,trutho,oriento,delxo,delyo
      parameter(projo=5.,trutho=45.,oriento=-100.,
     1          delxo=25000.,delyo=25000.)
c
      parameter(numtype=100)
      dimension ntype(numtype)
c
      dimension rlat(io,jo),rlon(io,jo)
c
      dimension rslmsk(imin*jmin)
      dimension slmask(io*jo)
c
      logical lmask
c
      real d00(io,jo),d10(io,jo),d11(io,jo),d01(io,jo)
      integer n00(io,jo),n10(io,jo),n11(io,jo),n01(io,jo)
      save n00,n10,n11,n01,d00,d10,d11,d01
c
      integer ij1(io),icl(io),jcl(io)
c
      integer ifp
      data ifp/0/
      save ifp    
c
#ifdef DBG
      print *,'inttyp=',inttyp
#endif
c
      iminjmin = imin*jmin
      ijo = io*jo
      if(imin.eq.1.or.jmin.eq.1) then
        print*,' wrong imin jmin in nfdr_ll2rsm',imin,jmin
        call abort
      endif
c
      if(delxi.eq.0..or.delyi.eq.0.) then
        print *,'delxi or delyi is zero in nfdr_ll2rsm'
        call abort
      endif
c
c  compute latitude and longitude of nfdr grid
c
      call nfdr_setgrd(rlat,rlon,delxo,delyo)
c
      if(inttyp.eq.2) then
        ipm=nint((delxo/2.)/(delxi*111.))
        jpm=nint((delyo/2.)/(delyi*111.))
#ifdef DBG
        print *,'ipm,jpm=',ipm,jpm
#endif
      endif
c
c     find i and j index and weighting factor
c
      rlat1o = rlat(1,1)
      rlat2o = rlat(io,jo)
      rlon1o = rlon(1,1)
      rlon2o = rlon(io,jo)
c
#ifdef DBG
      print *, ' prepare grid to grid interpolation.'
      print *, ' input grid : '
      print *, ' proji=',proji,' orienti=',orienti,
     1      ' truthi=',truthi,' cotrui=',cotrui,
     2      ' delxi=',delxi,' delyi=',delyi,
     3      ' rlat1i=',rlat1i,' rlon1i=',rlon1i,
     4      ' rlat2i=',rlat2i,' rlon2i=',rlon2i
      print *, ' output grid : '
      print *, ' projo=',projo,' oriento=',oriento,
     1      ' trutho=',trutho,' cotruo=',cotruo,
     2      ' delxo=',delxo,' delyo=',delyo,
     3      ' rlat1o=',rlat1o,' rlon1o=',rlon1o,
     4      ' rlat2o=',rlat2o,' rlon2o=',rlon2o
#endif
c
      if(abs(projo).eq.4.) then
        print *,'no regional gaussian grid allowed in nfdr_ll2rsm'
        call abort
      endif
c
      allocate (glon(imin),glat(jmin),glao(jo),tmp(io,jo))
c
      if(ifp.eq.0) then
        call i2oini(proji,orienti,truthi,cotrui,
     1              delxi,delyi,rlat1i,rlon1i,imin,jmin,
     2              projo,oriento,trutho,cotruo,
     3              x00,y00,x00o,y00o,
     4              rlat1o,rlon1o,jo,glao,glon,glat,gaul)
        do j=1,jo
          call i2o(io,j,proji,orienti,truthi,cotrui,
     1             delxi,delyi,rlat1i,rlon1i,imin,jmin,
     2             projo,oriento,trutho,cotruo,
     3             delxo,delyo,rlat1o,rlon1o,glao(j),
     4             x00,y00,x00o,y00o,
     5             n00(1,j),n10(1,j),n11(1,j),n01(1,j),
     6             d00(1,j),d10(1,j),d11(1,j),d01(1,j),
     7             glon,glat,gaul)
        enddo
        ifp=1
      endif
c
      sum1=0.
      sum2=0.
      wei1=0.
      wei2=0.
      j1 = 0
      j2 = imin*(jmin-1)
      do i=1,imin
        sum1=sum1+rdat(i+j1) * rslmsk(i+j1)
        sum2=sum2+rdat(i+j2) * rslmsk(i+j2)
        wei1=wei1+rslmsk(i+j1)
        wei2=wei2+rslmsk(i+j2)
      enddo
#ifdef DBG
      print *,' step wei12 ',wei1,wei2
#endif
      if(rlat1i.gt.0.) then
        if(wei1.gt.0.) then
          sumn=sum1/wei1
        else
          sumn=0.
        endif
        if(wei2.gt.0.) then
          sums=sum2/wei2
        else
          sums=0.
        endif
      else
        if(wei1.gt.0.) then
          sums=sum1/wei1
        else
          sums=0.
        endif
        if(wei2.gt.0.) then
          sumn=sum2/wei2
        else
          sumn=0.
        endif
      endif
c
      ifill=0
      do j=1,jo
        if(inttyp.eq.1.or.inttyp.eq.2) then
          do i=1,io
            weimax = 0.
            if(weimax.le.d00(i,j)) then
              weimax = d00(i,j)
              ij1(i) = n00(i,j)
            elseif(weimax.le.d01(i,j)) then
              weimax = d01(i,j)
              ij1(i) = n01(i,j)
            elseif(weimax.le.d10(i,j)) then
              weimax = d10(i,j)
              ij1(i) = n10(i,j)
            elseif(weimax.le.d11(i,j)) then
              weimax = d11(i,j)
              ij1(i) = n11(i,j)
            endif
            jcl(i)=(ij1(i)-1)/imin+1
            icl(i)=ij1(i)-(jcl(i)-1)*imin
          enddo
        endif
c
c  inttyp=1  take the closest point value
c
        if(inttyp.eq.1) then
          do i=1,io
            ij=(j-1)*io+i
            rout(i,j)=rdat(ij1(i))
          enddo
c
c  inttyp=2 take the predominant type
c
        elseif(inttyp.eq.2) then
          do i=1,io
            ij=(j-1)*io+i
            jcl(i)=(ij1(i)-1)/imin+1
            icl(i)=ij1(i)-(jcl(i)-1)*imin
            isrs=max(icl(i)-ipm,1)
            isrf=min(icl(i)+ipm,imin)
            jsrs=max(jcl(i)-jpm,1)
            jsrf=min(jcl(i)+jpm,jmin)
            do n=1,numtype
              ntype(n)=0
            enddo
            do jx=jsrs,jsrf
              do ix=isrs,isrf
                ixjx=(jx-1)*imin+ix
                n=rdat(ixjx)+1
                if(n.gt.numtype) then
                  print *,'type .gt.numtype.  change numtype in la2ga'
                  call abort
                endif
                if(n.lt.0) then
                  print *,'type cannot be less than zero'
                  call abort
                endif
                ntype(n)=ntype(n)+1
              enddo
            enddo
            maxcount=0
            maxtyp=-1
            do n=1,numtype
              if(ntype(n).gt.maxcount) then
                maxcount=ntype(n)
                maxtyp=n
              endif
            enddo
            rout(i,j)=maxtyp-1
          enddo
c
c  inttyp=0  bilinear interpolation without mask
c
        elseif(inttyp.eq.0.and..not.lmask) then
          do i=1,io
            if(n00(i,j).eq.-999) then
                rout(i,j)=-999.
                go to 200
            endif
            ij=(j-1)*io+i
            a00 = rdat(n00(i,j))
            a10 = rdat(n10(i,j))
            a11 = rdat(n11(i,j))
            a01 = rdat(n01(i,j))
            wi1j1 = d00(i,j)
            wi2j1 = d10(i,j)
            wi1j2 = d01(i,j)
            wi2j2 = d11(i,j)
            wsum  = wi1j1+wi2j1+wi1j2+wi2j2
            wsumiv = 1./wsum
            if(n00(i,j).ne.n01(i,j)) then
              rout(i,j)=(wi1j1*a00+wi2j1*a10+
     1                  wi1j2*a01+wi2j2*a11)*wsumiv
            else
              if(n00(i,j).le.imin.and.rlat1i.gt.0..or.n00(i,j).gt.j2.
     1           and.rlat1i.lt.0.) then
                  rout(i,j)=(wi1j1*sumn + wi2j1*sumn +
     1                   wi1j2*a01+wi2j2*a11)*wsumiv
              endif
              if(n00(i,j).le.imin.and.rlat1i.lt.0..or.n00(i,j).gt.j2.
     1           and.rlat1i.gt.0.) then
                rout(i,j)=(wi1j1*a00 + wi2j1*a10 +
     1                 wi1j2*sums + wi2j2*sums)*wsumiv
              endif
            endif
  200     continue
          enddo
c
c  inttyp=0  bilinear interpolation with mask
c
        elseif(inttyp.eq.0.and.lmask) then
          do i=1,io
            ij=(j-1)*io+i
            a00 = rdat(n00(i,j))
            a10 = rdat(n10(i,j))
            a11 = rdat(n11(i,j))
            a01 = rdat(n01(i,j))
            if(slmask(ij).eq.rslmsk(n00(i,j)).and.
     1         slmask(ij).eq.rslmsk(n01(i,j)).and.
     2         slmask(ij).eq.rslmsk(n10(i,j)).and.
     3         slmask(ij).eq.rslmsk(n11(i,j)))then
              wi1j1 = d00(i,j)
              wi2j1 = d10(i,j)
              wi1j2 = d01(i,j)
              wi2j2 = d11(i,j)
            elseif(slmask(ij).eq.1.) then
              wi1j1 = d00(i,j) * rslmsk(n00(i,j))
              wi2j1 = d10(i,j) * rslmsk(n10(i,j))
              wi1j2 = d01(i,j) * rslmsk(n01(i,j))
              wi2j2 = d11(i,j) * rslmsk(n11(i,j))
            elseif(slmask(ij).eq.0.) then
              wi1j1 = d00(i,j) * (1.-rslmsk(n00(i,j)))
              wi2j1 = d10(i,j) * (1.-rslmsk(n10(i,j)))
              wi1j2 = d01(i,j) * (1.-rslmsk(n01(i,j)))
              wi2j2 = d11(i,j) * (1.-rslmsk(n11(i,j)))
            endif
            wsum  = wi1j1+wi2j1+wi1j2+wi2j2
            if(wsum.ne.0.) then
              wsumiv = 1./wsum
              if(n00(i,j).ne.n01(i,j)) then
                rout(i,j)=(wi1j1*a00+wi2j1*a10+
     1                   wi1j2*a01+wi2j2*a11)*wsumiv
              else
                if(n00(i,j).le.imin.and.rlat1i.gt.0..or.n00(i,j).gt.j2.
     1             and.rlat1i.lt.0.) then
                    rout(i,j)=(wi1j1*sumn + wi2j1*sumn +
     1                     wi1j2*a01+wi2j2*a11)*wsumiv
                endif
                if(n00(i,j).le.imin.and.rlat1i.lt.0..or.n00(i,j).gt.j2.
     1             and.rlat1i.gt.0.) then
                  rout(i,j)=(wi1j1*a00 + wi2j1*a10 +
     1                     wi1j2*sums + wi2j2*sums)*wsumiv
                endif
              endif
            else
              ifill=ifill+1
#ifdef DBG
              if(ifill.le.2) then
                write(6,*) 'n00 n01 n10 n11 = ',
     1                      n00(i,j),n01(i,j),n10(i,j),n11(i,j)
                write(6,*) 'rslmsk=',rslmsk(n00(i,j)),rslmsk(n01(i,j)),
     1                      rslmsk(n10(i,j)),rslmsk(n11(i,j))
                write(6,*) 'ij = ',ij,' slmask(ij)=',slmask(ij)
              endif
#endif
              ix = n11(i,j)
              iy = n11(i,j)
  103         continue
              if(ix.gt.iminjmin) go to 101
              if((slmask(ij).eq.1..and.
     1            slmask(ij).eq.rslmsk(ix)).or.
     2           (slmask(ij).eq.0..and.
     3           slmask(ij).eq.rslmsk(ix))) then
                 rout(i,j) = rdat(ix)
                go to 100
              endif
              ix=ix+1
  101         continue
              if(iy.lt.1) go to 102
              if((slmask(ij).eq.1..and.
     1            slmask(ij).eq.rslmsk(iy)).or.
     2           (slmask(ij).eq.0..and.
     3            slmask(ij).eq.rslmsk(iy))) then
                 rout(i,j) = rdat(iy)
                go to 100
              endif
              iy=iy-1
  102         continue
              if(ix.le.iminjmin.or.iy.ge.1) go to 103
              write(6,*) ' error! no filling value found in nfdr_ll2rsm'
              print *,"debug print rslmsk"
              call nntprt(rslmsk,imin,jmin)
              print *,"debug print slmask"
              call nntprt(slmask,io,jo)
              call abort
            endif
  100       continue
c
          enddo
        endif
      enddo
c
      do j=1,jo
        jj=jo-j+1
        do i=1,io
          tmp(i,jj)=rout(i,j)
        enddo
      enddo
      do j=1,jo
        do i=1,io
          rout(i,j)=tmp(i,j)
        enddo
      enddo
     
c
      if(ifill.gt.1) then
        write(6,*) ' unable to interpolate.  filled with nearest',
     1             ' point value at ',ifill,' points'
      endif
c
#ifdef DBG
      call rmaxmin(rout,io*jo,1,1,1,'nfdr_ll2rsm:rout')
c
#endif
      deallocate (glon,glat,tmp)
      return
      end
