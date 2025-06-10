      subroutine sg2sg(ps1,s1,var1,ps2,s2,var2,im,jm,km1,km2,in,tensn,
     1                 iuv)
      save
c **********
c this subroutine transfer from one defined sigma coordinate (s1) to
c  another sigma coordinate (s2).
c input    ps1     primary ground pressure   alog(psfc)
c          ps2     secondary gound pressure alog(psfc(cb))
c          s1      primary sigma coordinate
c          var1    primary 3 dimensional variable
c          s2      secondary sigma coordinate
c          var2    secondary 3 dimensional variable
c          im      dimension in x
c          jm      dimension in y
c          km1     dimension in z for primary coordinate and variable
c          km2     dimension in z for secondary coordinate and variable
c          in      control index:   0 for initial transformation
c                                   1 for transformation as previous
c          tensn   factor of tension 0 for cubic spline
c                                   50 for linear interpolation
c          iuv     index for wind avoid extrapolation
c
c local dimensions: st1 st2 sv1 sv2 q vinc
c
c common block /spl/ ovh sh iflag theta
c
c routines related: trispl valts
c
c written by hann-ming henry juang
c date july 1988
c
      common /spl/ ovh(100), sh(100), theta
      common /spli/ iflag, jflag
c
      dimension s1(im,jm,km1),s2(im,jm,km2)
      dimension var1(im,jm,km1),var2(im,jm,km2)
      dimension st1(100),st2(100),st0(100),ps1(im,jm),ps2(im,jm)
      dimension sv1(100),sv2(100),q(100),vinc(100)
c
*
c
c theta is tension factor
      theta=tensn

      if( in .eq. 0 ) then
      in=1
      print *, ' initial values for transformation.'
      iflag = 0
      print *,'  iflag=',iflag,' tension factor=',theta
      do 1 k=1,km1
      do 1 j=1,jm
      do 1 i=1,im
      kk=km1-k+1
      if(s1(i,j,kk).gt.0.0.and.s1(i,j,kk).lt.1.0 ) then
        st0(k)=log(s1(i,j,kk))
      else
        print *, ' error in s1(i,j,k) vals @ k=',k,' sigma=',s1(i,j,kk)
        print *, ' error                        i=',i,' j=',j
        call abort
      endif
 1    continue
      do 2 k=1,km2
      do 2 j=1,jm
      do 2 i=1,im
      kk=km2-k+1
      if( s2(i,j,kk).gt.0.0.and.s2(i,j,kk).lt.1.0 ) then
       st2(k)=log(s2(i,j,kk))
      else
       print *, ' error in s2(i,j,k) vals @ k=',k,' sigma=',s2(i,j,kk)
       print *, ' error                        i=',i,' j=',j
       call abort
      endif
 2    continue
      vst=st2(1)
      do 3 k=1,km2-1
      vinc(k)=st2(k+1)-st2(k)
 3    continue
      endif

      do 10 i=1,im
      do 10 j=1,jm

      do 20 k=1,km1
      kk=km1-k+1
      st1(k)=log(s1(i,j,kk))+ps1(i,j)-ps2(i,j)
      sv1(k)=var1(i,j,km1-k+1)
 20   continue
      do 21 k=1,km2
      kk=km2-k+1
      st2(k)=log(s2(i,j,kk))
 21   continue
      vst=st2(1)
      do 22 k=1,km2-1
      vinc(k)=st2(k+1)-st2(k)
 22   continue

      iflag=0
      call trispl(km1,st1,q,sv1)
      call valts(q,vst,vinc,st1,sv1,km2,km1,sv2)

c      if ((i.eq.193).and.(j.eq.1)) then
c        print'(2i5,2f15.7)',i,k,ps1(i,j),ps2(i,j)
c        do k=1,km2
c          print'(2i5,5f15.7)',i,k,sv1(k),sv2(k),st1(k),st2(k),q(k)
c        enddo
c      endif

      kmm=km2/2 + 1
      do 30 k=kmm,1,-1
      if( sv2(k) .eq. 99999.9 ) then
c
c  set to the lowest primary sigma level value
c
      do kk=1,km1
        if(st1(kk).gt.st2(k)) then
          diff1=st1(kk)-st2(k)
          if(kk.gt.1) then
            diff2=st1(kk-1)-st2(k)
          else
            diff2=999.
          endif
          if(diff1.lt.diff2) then
            sv2(k)=sv1(kk)
          else
            sv2(k)=sv1(kk-1)
          endif
          go to 700
        endif
      enddo
  700 continue
      endif
 30   continue
      do 40 k=kmm,km2
      if( sv2(k) .eq. 99999.9 ) then
c
c  set to the lowest primary sigma level value
c
      do kk=1,km1
        if(st1(kk).gt.st2(k)) then
          diff1=st1(kk)-st2(k)
          if(kk.gt.1) then
            diff2=st1(kk-1)-st2(k)
          else
            diff2=999.
          endif
          if(diff1.lt.diff2) then
            sv2(k)=sv1(kk)
          else
            sv2(k)=sv1(kk-1)
          endif
          go to 710
        else
          sv2(k)=sv1(km1)
        endif
      enddo
  710 continue
      endif
 40   continue
      do 50 k=1,km2
      var2(i,j,k)=sv2(km2-k+1)
 50   continue
 10   continue
      return
      end
