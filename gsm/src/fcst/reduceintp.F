      subroutine reduceintp(a,lonfd,lonf2,latg2)
c
c return the reduce-grid value to the nearest regular grid
c then, the restored value is propogated eastward to the regular
c grid which has un-specified value.
c 				by h.-m. h. juang (2000)
c
      dimension a(lonf2,latg2)
      dimension lonfd(latg2)
      real,allocatable:: tmp(:)
      allocate(tmp(lonf2))
      lonf=lonf2/2
      dg=360./float(lonf)
      do j=1,latg2
        do i=1,lonf2
          tmp(i)=-999.999
        enddo
        dr=360./float(lonfd(j))
        do i=1,lonfd(j)
          ii=nint((i-1.)*dr/dg + 1.0)
          tmp(ii)=a(i,j)
          tmp(ii+lonf)=a(i+lonfd(j),j)
        enddo
        do i=1,lonf
          if( tmp(i).ne.-999.999 ) rep=tmp(i)
          tmp(i)=rep
        enddo
        do i=1,lonf
          if( tmp(i+lonf).ne.-999.999 ) rep=tmp(i+lonf)
          tmp(i+lonf)=rep
        enddo
        do i=1,lonf2
          a(i,j)=tmp(i)
        enddo
      enddo
      deallocate(tmp)
      return
      end
