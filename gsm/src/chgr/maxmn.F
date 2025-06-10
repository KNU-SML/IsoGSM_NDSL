      subroutine maxmn(f,imax,jmax,kmax,title)
c
      dimension f(imax,jmax,kmax)
      character*8 title
c
      print 99, title
   99 format(2x,'title=',a8)
c
      do 10 k=1,kmax
c
      fmax=f(1,1,k)
      fmin=f(1,1,k)
c
      do 20 j=1,jmax
      do 20 i=1,imax
      if(fmax.lt.f(i,j,k)) then
      fmax=f(i,j,k)
      iimax=i
      jjmax=j
      endif
      if(fmin.gt.f(i,j,k)) then
      fmin=f(i,j,k)
      iimin=i
      jjmin=j
      endif
   20 continue
c
      print 100, k,fmax,iimax,jjmax,fmin,iimin,jjmin
  100 format(2x,'level=',i3,' max=',e12.4,' at i=',i5,' j=',i5,
     1                      ' min=',e12.4,' at i=',i5,' j=',i5)
c
   10 continue
c
      return
      end
