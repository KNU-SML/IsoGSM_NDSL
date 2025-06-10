      subroutine la2ga(regin ,imxin ,jmxin , dloin,dlain, rlon, rlat,
     1                 gauout,imxout,jmxout,lmask,rslmsk,slmask,igau,
     2                 inttyp)
#include <paramodel.h>
c
c  interpolation from lat/lon or gaussian grid to other gaussian grid
c
c  igau=0  lat/lon input
c  igau=1  gaussian input
c
c  inttyp=0  regular interpolation
c  inttyp=1  take the closest point value
c  inttyp=2  take the dominant type
c
      save
c
      parameter(numtype=20)
      dimension ntype(numtype)
c
      dimension regin (imxin ,jmxin )
      dimension gauout(imxout,jmxout)
c
      dimension rslmsk(imxin,jmxin)
      dimension slmask(imxout,jmxout)
c
      dimension gaulo (jmxout)
      dimension gauli (jmxin)
      dimension rinlat(jmxin),outlat(jmxout)
      dimension rinlon(imxin)
c
      dimension iindx1(imxout)
      dimension iindx2(imxout)
      dimension jindx1(jmxout)
      dimension jindx2(jmxout)
c
      dimension ddx(imxout)
      dimension ddy(jmxout)
c
      logical lmask
c
      print *,'inttyp=',inttyp
c
      if(imxin.eq.1.or.jmxin.eq.1) then
        do j=1,jmxout
           do i=1,imxout
            gauout(i,j)=0.
           enddo
        enddo
        return
      endif
c
      if(dloin.eq.0..or.dlain.eq.0.) then
        print *,'dloin or dlain is zero .... check data cards'
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
      ifpi=imxin
      jfpi=jmxin
      ifpo=imxout
      jfpo=jmxout
      rfp1=rlon
      rfp2=rlat
c
      if(igau.eq.0) then
        do j=1,jmxin
          if(rlat.gt.0.) then
            rinlat(j)=rlat-float(j-1)*dlain
          else
            rinlat(j)=rlat+float(j-1)*dlain
          endif
        enddo
      else
        call gaulat(gauli,jmxin)
        do j=1,jmxin
          rinlat(j)=90.-gauli(j)
        enddo
        rlat=90.
      endif
c
c    compute gaussian latitude for output grid
c
      call gaulat(gaulo,jmxout)
c
      do j=1,jmxout
        outlat(j)=90.-gaulo(j)
      enddo
c
      do i=1,imxin
        rinlon(i)=rlon+float(i-1)*dloin
      enddo
c
c  find i-index for interplation
c
      do 30 i=1,imxout
      alamd=float(i-1)*360./float(imxout)
      if(rlon.lt.0.) then
      if(alamd.gt.180.) alamd=alamd-360.
      endif
      do 35 ii=1,imxin
      if(alamd.gt.rinlon(ii)) go to 35
      ix=ii
      go to 32
   35 continue
      i1=nint(360./dloin)
      i2=1
      go to 34
   32 continue
      if(ix.ge.2) go to 33
      i1=nint(360./dloin)
      i2=1
      go to 34
   33 continue
      i2=ix
      i1=i2-1
   34 continue
      iindx1(i)=i1
      iindx2(i)=i2
      denom=rinlon(i2)-rinlon(i1)
      if(denom.lt.0.) denom=denom+360.
      rnume=alamd-rinlon(i1)
      if(rnume.lt.0.) rnume=rnume+360.
      ddx(i)=rnume/denom
   30 continue
c
c  find j-index for interpolation
c
      jq=1
      do 40 j=1,jmxout
      aphi=outlat(j)
      do 50 jj=1,jmxin
      jx=jj
      if(rlat.lt.0.) jx=jmxin-jj+1
      if(aphi.lt.rinlat(jx)) go to 50
      jq=jx
      go to 42
   50 continue
      if(rlat.gt.0.) then
        j1=jmxin
        j2=jmxin
      else
        j1=1
        j2=1
      endif
      go to 44
   42 continue
      if(rlat.gt.0.) then
         if(jq.ge.2) go to 43
         j1=1
         j2=1
      else
         if(jq.lt.jmxin) go to 43
         j1=jmxin
         j2=jmxin
      endif
      go to 44
   43 continue
      if(rlat.gt.0.) then
      j2=jq
      j1=jq-1
      else
      j1=jq
      j2=jq+1
      endif
   44 continue
      jindx1(j)=j1
      jindx2(j)=j2
      if(j2.ne.j1) then
         ddy(j)=(aphi-rinlat(j1))/(rinlat(j2)-rinlat(j1))
      else
      if(j1.eq.1.and.rlat.gt.0..or.j1.eq.jmxin.and.rlat.lt.0.) then
            if(abs(90.-rinlat(j1)).gt.0.001) then
               ddy(j)=(aphi-rinlat(j1))/(90.-rinlat(j1))
            else
               ddy(j)=0.0
            endif
      endif
      if(j1.eq.1.and.rlat.lt.0..or.j1.eq.jmxin.and.rlat.gt.0.) then
            if(abs(-90.-rinlat(j1)).gt.0.001) then
               ddy(j)=(aphi-rinlat(j1))/(-90.-rinlat(j1))
            else
               ddy(j)=0.0
            endif
         endif
      endif
   40 continue
c
c increments of index for dominant type search
c
      if(inttyp.eq.2) then
        ipm=imxin/imxout/2
        jpm=jmxin/jmxout/2
        print *,'ipm=',ipm,' jpm=',jpm
      endif
c
      sum1=0.
      sum2=0.
      wei1=0.
      wei2=0.
      do i=1,imxin
        sum1=sum1+regin(i,1) * rslmsk(i,1)
        sum2=sum2+regin(i,jmxin) * rslmsk(i,jmxin)
        wei1=wei1+rslmsk(i,1)
        wei2=wei2+rslmsk(i,jmxin)
      enddo
      if(rlat.gt.0.) then
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
c  quasi-bilinear interpolation or closest point or dominant type
c
      ifill=0
      do 70 j=1,jmxout
      y=ddy(j)
      j1=jindx1(j)
      j2=jindx2(j)
      jcl=j1
      if(y.gt.0.5) jcl=j2
      do 70 i=1,imxout
      x=ddx(i)
      i1=iindx1(i)
      i2=iindx2(i)
      icl=i1
      if(x.gt.0.5) icl=i2
c
      if(inttyp.le.1) then
        if(lmask) then
          if(slmask(i,j).eq.rslmsk(i1,j1).and.
     1       slmask(i,j).eq.rslmsk(i2,j1).and.
     2       slmask(i,j).eq.rslmsk(i1,j2).and.
     3       slmask(i,j).eq.rslmsk(i2,j2)) then
            wi1j1=(1.-x)*(1.-y)
            wi2j1=    x *(1.-y)
            wi1j2=(1.-x)*      y
            wi2j2=    x *      y
          elseif(slmask(i,j).eq.1.) then
            wi1j1=(1.-x)*(1.-y)  *rslmsk(i1,j1)
            wi2j1=    x *(1.-y)  *rslmsk(i2,j1)
            wi1j2=(1.-x)*      y *rslmsk(i1,j2)
            wi2j2=    x *      y *rslmsk(i2,j2)
          elseif(slmask(i,j).eq.0.) then
            wi1j1=(1.-x)*(1.-y)  *(1.-rslmsk(i1,j1))
            wi2j1=    x *(1.-y)  *(1.-rslmsk(i2,j1))
            wi1j2=(1.-x)*      y *(1.-rslmsk(i1,j2))
            wi2j2=    x *      y *(1.-rslmsk(i2,j2))
          endif
        else
          wi1j1=(1.-x)*(1.-y)
          wi2j1=    x *(1.-y)
          wi1j2=(1.-x)*      y
          wi2j2=    x *      y
        endif
c
        if(inttyp.eq.1) then
           if(icl.eq.i1.and.jcl.eq.j1) then
             wi1j1=1.
             wi2j1=0.
             wi1j2=0.
             wi2j2=0.
           elseif(icl.eq.i2.and.jcl.eq.j1) then
             wi1j1=0.
             wi2j1=1.
             wi1j2=0.
             wi2j2=0.
           elseif(icl.eq.i1.and.jcl.eq.j2) then
             wi1j1=0.
             wi2j1=0.
             wi1j2=1.
             wi2j2=0.
           elseif(icl.eq.i2.and.jcl.eq.j2) then
             wi1j1=0.
             wi2j1=0.
             wi1j2=0.
             wi2j2=1.
           endif
        endif
c
        wsum=wi1j1+wi2j1+wi1j2+wi2j2
        if(wsum.ne.0.) then
          wsumiv = 1./wsum
c
          if(j1.ne.j2) then
            gauout(i,j)=(wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+
     1                   wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))*wsumiv
          else
           if(j1.eq.1.and.rlat.gt.0..or.j1.eq.jmxin.and.rlat.lt.0.) then
            gauout(i,j)=(wi1j1*sumn        +wi2j1*sumn        +
     1                   wi1j2*regin(i1,j2)+wi2j2*regin(i2,j2))*wsumiv
           endif
           if(j1.eq.1.and.rlat.lt.0..or.j1.eq.jmxin.and.rlat.gt.0.) then
            gauout(i,j)=(wi1j1*regin(i1,j1)+wi2j1*regin(i2,j1)+
     1                   wi1j2*sums        +wi2j2*sums        )*wsumiv
           endif
          endif
        else
          if(.not.lmask) then
            write(6,*) ' la2ga called with lmask=.true. but ',
     1               ' rslmsk or slmask bad.'
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
          ifill=ifill+1
c
          do jx=j1,jmxin
            do ix=i1,imxin
              if((slmask(i,j).eq.1..and.
     1            slmask(i,j).eq.rslmsk(ix,jx)).or.
     2           (slmask(i,j).eq.0..and.
     3            slmask(i,j).eq.rslmsk(ix,jx))) then
                gauout(i,j)=regin(ix,jx)
                go to 71
              endif
            enddo
            do ix=i1,1,-1
              if((slmask(i,j).eq.1..and.
     1          slmask(i,j).eq.rslmsk(ix,jx)).or.
     2         (slmask(i,j).eq.0..and.
     3          slmask(i,j).eq.rslmsk(ix,jx))) then
                gauout(i,j)=regin(ix,jx)
                go to 71
              endif
            enddo
          enddo
          do jx=j1,1,-1
            do ix=i1,imxin
              if((slmask(i,j).eq.1..and.
     1            slmask(i,j).eq.rslmsk(ix,jx)).or.
     2           (slmask(i,j).eq.0..and.
     3            slmask(i,j).eq.rslmsk(ix,jx))) then
                gauout(i,j)=regin(ix,jx)
                go to 71
              endif
            enddo
            do ix=i1,1,-1
              if((slmask(i,j).eq.1..and.
     1            slmask(i,j).eq.rslmsk(ix,jx)).or.
     2           (slmask(i,j).eq.0..and.
     3            slmask(i,j).eq.rslmsk(ix,jx))) then
                gauout(i,j)=regin(ix,jx)
                go to 71
              endif
            enddo
          enddo
          write(6,*) ' error!!! no filling value found in la2ga'
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
c  dominant type
c
      elseif(inttyp.eq.2) then
        isrs=icl-ipm
        isrf=icl+ipm
        jsrs=max(jcl-jpm,1)
        jsrf=min(jcl+jpm,jmxin)
        do n=1,numtype
          ntype(n)=0
        enddo
        do jx=jsrs,jsrf
          do ixx=isrs,isrf
            ix=ixx
            if(ix.lt.1) ix=imxin-ixx
            if(ix.gt.imxin) ix=ixx-imxin
            n=regin(ix,jx)+1
            if(n.gt.numtype) then
              print *,'type .gt.numtype.  change numtype in la2ga'
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
            if(n.lt.0) then
              print *,'type cannot be less than zero'
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
        gauout(i,j)=maxtyp-1
      endif
c
   71 continue
   70 continue
c
      if(ifill.gt.1) then
        write(6,*) ' unable to interpolate.  filled with nearest',
     1             ' point value at ',ifill,' points'
      endif
c
c#ifdef DBG
      write(6,*) ' rmaxmin of model grid arry:'
      call rmaxmin(gauout,imxout*jmxout,1,1,1,'gauout')
c#endif
c
      return
      end
