      subroutine gl2gl(gauin,imxin,jmxin,gauout,imxout,jmxout)
c
c  interpolation from gaussian grid to other gaussian grid
c
      dimension gauin (imxin ,jmxin )
c
      dimension gauout(imxout,jmxout)
      dimension gaulin(500),gaulou(500)
      dimension iindx1(1000)
      dimension iindx2(1000)
      dimension jindx1(500)
      dimension jindx2(500)
      dimension ddx(1000)
      dimension ddy(500)
c
      save iindx1,iindx2,jindx1,jindx2,ddx,ddy
c
      data ifp/0/
c
      if(imxin.eq.imxout.and.jmxin.eq.jmxout) then
        do 5 j=1,jmxin
        do 5 i=1,imxin
         gauout(i,j)=gauin(i,j)
    5   continue
        return
      endif
c
      if(ifp.ne.0) go to 111
      ifp=1
c
      call gaulat(gaulin,jmxin)
      do 10 j=1,jmxin
      gaulin(j)=90.-gaulin(j)
   10 continue
c
      call gaulat(gaulou,jmxout)
      do 20 j=1,jmxout
      gaulou(j)=90.-gaulou(j)
   20 continue
c
      dxin =360./float(imxin )
      dxout=360./float(imxout)
c
      do 30 i=1,imxout
      alamd=float(i-1)*dxout
      i1=alamd/dxin+1.001
      iindx1(i)=i1
      i2=i1+1
      if(i2.gt.imxin) i2=1
      iindx2(i)=i2
c     print *,'i,i1,i2',i,i1,i2
      ddx(i)=(alamd-float(i1-1)*dxin)/dxin
   30 continue
c
      j2=1
      do 40 j=1,jmxout
      aphi=gaulou(j)
      do 50 jj=1,jmxin
      if(aphi.lt.gaulin(jj)) go to 50
      j2=jj
      go to 42
   50 continue
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
c     print *,'j,j1,j2',j,j1,j2
      jindx1(j)=j1
      jindx2(j)=j2
      ddy(j)=(aphi-gaulin(j1))/(gaulin(j2)-gaulin(j1))
   40 continue
c
c     print *,'iindx1'
c     print *,(iindx1(n),n=1,imxout)
c     print *,'iindx2'
c     print *,(iindx2(n),n=1,imxout)
c     print *,'jindx1'
c     print *,(jindx1(n),n=1,jmxout)
c     print *,'jindx2'
c     print *,(jindx2(n),n=1,jmxout)
c     print *,'ddy'
c     print *,(ddy(n),n=1,jmxout)
c     print *,'ddx'
c     print *,(ddx(n),n=1,jmxout)
c
  111 continue
c
      do 60 j=1,jmxout
      y=ddy(j)
      j1=jindx1(j)
      j2=jindx2(j)
      do 60 i=1,imxout
      x=ddx(i)
      i1=iindx1(i)
      i2=iindx2(i)
      gauout(i,j)=(1.-x)*(1.-y)*gauin(i1,j1)+(1.-y)*x*gauin(i2,j1)+
     1           (1.-x)*y*gauin(i1,j2)+x*y*gauin(i2,j2)
   60 continue
c
      return
      end
