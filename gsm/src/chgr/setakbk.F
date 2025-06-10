      subroutine setakbk(km,si,sl,ak,bk,del)
c
c  for sigma pressure re: hybrids implemented may 1st 2007 @ 12Z
c
      dimension si(km+1),sl(km),ak(km+1),bk(km+1),del(km)
c
c     print *,' set ak and bk to related to centibar '
c
      do k=1,km
        ak(k) = si(k) / 1000.  ! Pa -> centibar 
        bk(k) = sl(k)
      enddo
      ak(km+1) = 0.0
      bk(km+1) = 0.0
      do k=1,km
        del(k) = (bk(k)+ak(k)/100)-(bk(k+1)+ak(k+1)/100) ! temporary depth
      enddo
      do k=1,km
       print *,' k ak bk del',k,ak(k),bk(k),del(k)
      enddo
      return
      end
