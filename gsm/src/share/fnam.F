      subroutine fnam(fnin,nchin,fhour,fnout,nchout)
c
      character*10 cfhour
      character*80 fnin,fnout
c
      do n=1,80
        fnout(n:n)=' '
      enddo
c
      ifh=fhour
      if(ifh.lt.10) then
        ndig=2
        write(cfhour,'(1h0,i1)') ifh
      elseif(ifh.lt.100) then
        ndig=2
        write(cfhour,'(i2)') ifh
      elseif(ifh.lt.1000) then
        ndig=3
        write(cfhour,'(i3)') ifh
      elseif(ifh.lt.10000) then
        ndig=4
        write(cfhour,'(i4)') ifh
      elseif(ifh.lt.100000) then
        ndig=5
        write(cfhour,'(i5)') ifh
      elseif(ifh.lt.1000000) then
        ndig=6
        write(cfhour,'(i6)') ifh
      elseif(ifh.lt.10000000) then
        ndig=7
        write(cfhour,'(i7)') ifh
      else
        print *,'ifh=',ifh, 'is too large'
        call abort
      endif
c     fnout=fnin(1:nchin)//'.ft'//cfhour(1:ndig)
      fnout(1:nchin)=fnin(1:nchin)
      fnout(nchin+1:nchin+3)='.ft'
      fnout(1+nchin+3:ndig+nchin+3)=cfhour(1:ndig)
      do n=ndig+nchin+3+1,80
        fnout(n:n)=' '
      enddo
      nchout=nchin+3+ndig
c
      return 
      end
