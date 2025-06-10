      subroutine reducepick(a,lonfd,lonf2,latg2)
c
c pick the reduce-grid value from the nearest regular grid
c then fill the tailing points to be the same as the last
c reduce-grid point.
c 				by h.-m. h. juang (2000)
c
      dimension a(lonf2,latg2)
      dimension lonfd(latg2)
      real,allocatable:: tmp(:)
      allocate(tmp(lonf2))
      lonf=lonf2/2
      dg=360./float(lonf)
      do j=1,latg2
        dr=360./float(lonfd(j))
        do i=1,lonfd(j)
          ii=nint((i-1.)*dr/dg + 1.0)
          tmp(i)=a(ii,j)
          tmp(i+lonfd(j))=a(ii+lonf,j)
        enddo
        ii=lonfd(j)*2
        do i=1,ii
          a(i,j)=tmp(i)
        enddo
        do i=ii+1,lonf2
          a(i,j)=tmp(ii)
        enddo
      enddo
      deallocate(tmp)
      return
      end
