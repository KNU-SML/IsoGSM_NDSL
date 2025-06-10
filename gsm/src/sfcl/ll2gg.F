      subroutine ll2gg(io,proj,glat,glon,ii,jj,delx,dely,rlat1,rlon1,
     1                 rlat,rlon,xlon,ylat)
c
      dimension glat(jj),glon(ii)
      real rlat(io),rlon(io),xlon(io),ylat(io)
c
      integer, allocatable :: lon(:)
      logical, allocatable :: flag(:)
      real, allocatable :: slon(:)
c
      allocate(lon (io))
      allocate(flag(io))
      allocate(slon(io))
c
      if(abs(proj).eq.4.) then
        delx=360./float(ii)*1000.
      endif
c
      do i=1,io
        slon(i)=rlon(i)
        if(slon(i).lt.0.) slon(i)=slon(i)+360.
      enddo
c
      ibeg = 1
      glon0 = glon(1)
      do ix = 1,ii
         if(glon(ix).le.glon0 ) then
            glon0 = glon(ix)
            ibeg = ix
         endif
      enddo
      iend = ibeg + ii - 1
c
      do i=1,io
        flag(i)=.true.
      enddo
      do ix = ibeg,iend
        if(ix.le.ii) then
          iglon = ix
        else
          iglon = ix - ii
        endif
        gloni = glon(iglon)
        do i=1,io
          if(flag(i).and.slon(i).lt.gloni) then
             xlon(i) = iglon + (slon(i)-gloni)/delx*1000.
c            xlon(i) = iglon + (slon(i)-gloni)/delx
             lon(i) = xlon(i)
             flag(i) = .false.
          endif
        enddo
      enddo
      do ix = ibeg,iend
        if(iend.le.ii) then
          do i=1,io
            if((.not.flag(i)).and.lon(i).lt.1) then
              lon(i) = iend
            endif
          enddo
        else
          do i=1,io
            if((.not.flag(i)).and.lon(i).lt.1) then
              lon(i) = iend - ii
            endif
          enddo
        endif
        do i=1,io
          if((.not.flag(i)).and.lon(i).lt.1) then
            xlon(i) = iend + xlon(i)
          endif
        enddo
      enddo
      do ix = ibeg,iend
        if(iend.le.ii) then
          do i=1,io
            if(flag(i)) then
              lon(i) = iend
            endif
          enddo
        else
          do i=1,io
            if(flag(i)) then
              lon(i) = iend-ii
            endif
          enddo
        endif
      enddo
      do ix = ibeg,iend
        do i=1,io
          if(flag(i)) then
            xlon(i) = lon(i) + (slon(i)-glon(lon(i)))/delx*1000.
c           xlon(i) = lon(i) + (slon(i)-glon(lon(i)))/delx
          endif
        enddo
      enddo
c
      do i=1,io
        flag(i)=.true.
      enddo
      if(proj.gt.0.) then
        do j = jj,2,-1
          do i=1,io
            if(flag(i).and.glat(j).le.rlat(i)) then
              ylat(i) = j + (rlat(i)-glat(j))/(glat(j+1)-glat(j))
              flag(i) = .false.
            endif
            if(rlat(i).le.glat(1)) then
              ylat(i) = 1
              flag(i) = .false.
            endif
          enddo
        enddo
      else
        do j = 2,jj
          do i=1,io
            if(flag(i).and.glat(j).le.rlat(i)) then
              ylat(i) = j-1 + (glat(j-1)-rlat(i))/(glat(j-1)-glat(j))
              flag(i) = .false.
            endif
            if(rlat(i).ge.glat(1)) then
              ylat(i) = 1
              flag(i) = .false.
            endif
          enddo
        enddo
      endif
c
      deallocate (lon,flag,slon)
      return
      end
