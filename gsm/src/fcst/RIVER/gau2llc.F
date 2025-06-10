#define IO2 io2
#define JO2 jo2
      subroutine gau2llc(gauin,imxin,jmxin,
     1                  xlonw,xlatn,dxlon,dxlat,regout,imxout,jmxout)
#include <paramodel.h>
c      parameter(io2_=io_,jo2_=jo_)
c      parameter(io2_=360,jo2_=181) !! for 1-deg but not T126
c      parameter(io2_=720,jo2_=361) !! for 0.5-deg but not T248
c
c  interpolation from lat/lon grid to other lat/lon grid
c
      dimension gauin (imxin,jmxin)
      dimension regout(imxout,jmxout)
c
      dimension gaul(jmxin),regl(jmxout)
c
      real iindx1(io2_),iindx2(io2_),
     1               jindx1(jo2_-1),jindx2(jo2_-1),
     2                  ddx(io2_),   ddy(jo2_-1)
      save iindx1,iindx2,jindx1,jindx2,ddx,ddy
c
      save ifp
      data ifp/0/
c
      if(ifp.ne.0) go to 111
      ifp=1
c
c      write(6,*) 'imxin=',imxin,' jmxin=',jmxin
c      write(6,*) 'xlatn=',xlatn,' xlonw=',xlonw,
c     1           ' dxlat=',dxlat,' dxlon=',dxlon
c      write(6,*) 'imxout=',imxout,' jmxout=',jmxout
c
      call gaulat(gaul,jmxin)
      do j=1,jmxin
        gaul(j)=90.-gaul(j)
      enddo
c      print *,'gaul=',gaul
c
      do j=1,jmxout
        regl(j)=xlatn-float(j-1)*dxlat
      enddo
c
      do 40 j=1,jmxout
      aphi=regl(j)
      do 50 jj=1,jmxin
      if(aphi.lt.gaul(jj)) go to 50
      j2=jj
      go to 42
   50 continue
      j2=jmxin
   42 continue
      if(j2.gt.2) go to 43
      j1=1
      j2=2
      go to 44
   43 continue
      if(j2.le.jmxin) go to 45
      j1=jmxin-1
      j2=jmxin
      go to 44
   45 continue
      j1=j2-1
   44 continue
      jindx1(j)=j1
      jindx2(j)=j2
      ddy(j)=(aphi-gaul(j1))/(gaul(j2)-gaul(j1))
   40 continue
c
      dxin =360./float(imxin)
c
      do i=1,imxout
        alamd=xlonw+float(i-1)*dxlon
        if(alamd.lt.0.) alamd=360.+alamd
        i1=alamd/dxin+1.001
        if(i1.gt.imxin) i1=1
        iindx1(i)=i1
        i2=i1+1
        if(i2.gt.imxin) i2=1
        iindx2(i)=i2
        ddx(i)=(alamd-float(i1-1)*dxin)/dxin
c        print*,'ddx',alamd,ddx(i)
      enddo
c
c
#ifdef DBG
      write(6,*) 'iindx1'
      write(6,*) (iindx1(n),n=1,imxout)
      write(6,*) 'iindx2'
      write(6,*) (iindx2(n),n=1,imxout)
      write(6,*) 'jindx1'
      write(6,*) (jindx1(n),n=1,jmxout)
      write(6,*) 'jindx2'
      write(6,*) (jindx2(n),n=1,jmxout)
      write(6,*) 'ddy'
      write(6,*) (ddy(n),n=1,jmxout)
      write(6,*) 'ddx'
      write(6,*) (ddx(n),n=1,jmxout)
#endif
  111 continue
c
      do j=1,jmxout
        y=ddy(j)
        j1=jindx1(j)
        j2=jindx2(j)
        do i=1,imxout
          x=ddx(i)
          i1=iindx1(i)
          i2=iindx2(i)
          regout(i,j)=(1.-x)*(1.-y)*gauin(i1,j1)+
     1                (1.-y)*x*gauin(i2,j1)+
     2                (1.-x)*y*gauin(i1,j2)+x*y*gauin(i2,j2)
        enddo
      enddo
c
c      sum1=0.
c      sum2=0.
c      do i=1,imxin
c        sum1=sum1+gauin(i,1)
c        sum2=sum2+gauin(i,jmxin)
c      enddo
c      sum1=sum1/float(imxin)
c      sum2=sum2/float(imxin)
c
c      do i=1,imxout
c        if(abs(regl(1)).eq.90.) then
c          regout(i,     1)=sum1
c        endif
c        if(abs(regl(jmxout)).eq.90.) then
c          regout(i,jmxout)=sum2
c        endif
c      enddo
c
      return
      end
