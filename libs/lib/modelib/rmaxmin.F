      subroutine rmaxmin(a,len,k,k1,k2,ch)
c$$$  subprogram documentation block
c                .      .    .                                       .
c abstract:  do print maximum and minimum of a given array.
c
c program history log:
c
c usage:    call  maxmin(a,len,k,k1,k2,ch)
c   input argument list:
c     a 	- array for computing max and min (len,k)
c     len	- the first dimension of a
c     k 	- the second dimension of a
c     k1	- lower limit of second dimension to print
c     k2	- upper limit to print
c     ch 	- charcter string to print
c                 fpl and fml
c
c   output argument list:
c
c   input files: none
c
c   output files:
c     standard output
c
c   subprograms called:
c     intrinsic functions: amax1 amin1
c
c   remark: none
c
c attributes:
c   language: fortran 77.
c   machine:  cray ymp.
c
c$$$

      dimension a(len,k)
      character ch*(*)
c
      parameter(crit=1.e20)
c
      aamax=0.
      aamin=0.
      do j=k1,k2
        do m=1,len
          if(abs(a(m,j)).lt.crit.and.a(m,j).ne.0.) then
            aamax = a(m,j)
            aamin = a(m,j)
            go to 2
          endif
        enddo
    2   continue
        nmis=0
        nzero=0
        do m=1,len
          if(a(m,j).eq.0.) then
            nzero=nzero+1
            go to 1
          endif
          if(abs(a(m,j)).gt.crit) then
            nmis=nmis+1
            go to 1
          endif
          aamax = max( aamax, a(m,j) )
          aamin = min( aamin, a(m,j) )
    1     continue
        enddo
        print 100,ch,aamax,aamin,j,nzero,nmis
  100   format(a12,' has max=',e10.4,' min=',e10.4,
     1             ' k=',i4,' nzero=',i7,' nmis=',i7)
      enddo
      return
      end
