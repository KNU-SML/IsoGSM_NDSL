      subroutine valts(q,vst,vinc,p,y,m,n,yn)
      save
c-----------------------------------------------------
c     subroutine calculates interpolated values of function using
c        weights calcuated by trispl
c     m...number of grid points desired
c     yn...interpolated value of y(p) at grid points
c     dyn..c     dyn...analytic value of derivative at grid point
c     vst...starting location of interval
c     vinc...increment of regular grid
c     n...number of points input from trispl
c     q...weight function from trispl
c     ovh,sh...weights from trispl
c---------------------------------------------------
      common/spl/ ovh(100),sh(100),theta
      common/spli/ iflag,jflag
      dimension q(n),p(n),y(n),yn(m),vinc(m)
      vs=vst
      j=1
      i=1
    1 if(vs .gt. p(1)) go to 9
      yn(i)=99999.9
      i=i+1
      vs =vs +vinc(i-1)
      go to 1
c
    9 j=j+1
      if (j.gt.n) go to 90
    8 if (vs .gt.p(j)) go to 9
      jm=j-1
      yn(i)=q(jm)*sinh(theta*(p(j)-vs))*sh(jm)+(y(jm)-q(jm))*(p(j)-vs)
     +*ovh(jm)+q(j)*sinh(theta*(vs-p(jm)))*sh(jm)+(y(j)-q(j))*(vs-
     + p(jm))*ovh(jm)
      i=i+1
      vs =vs +vinc(i-1)
      if(m.ge.i) go to 8
      return
   90 yn(i)=99999.9
      i=i+1
      if (m.ge.i) go to 90
      return
      end
