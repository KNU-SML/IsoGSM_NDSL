      subroutine spcshfli(a,lnt,ntotal,jcap,lwvdef)
c$$$  subprogram documentation block
c
c subprogram: 	spcshfli
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: shafle spectral coefficients for balancing from input
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call spcshfli(a,lnt,ntotal,jcap,lwvdef)
c
c    input argument lists:
c	a	- real (lnt,ntotal) total spectral field 
c	lnt	- integer total spectral grid
c	jcap	- integer total wavenumber of truncation
c	ntotal	- integer total set of fields
c	lwvdef	- integer (jcap+1) distribution index
c
c    output argument list:
c	a	- real (lnt,ntotal) total field
c 
c subprograms called:
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
      implicit none
      integer lnt,lntpp,ntotal,jcap,lwvdef
     1       ,offset,l,n,m,k,j,ll
      real    a
      dimension a(lnt,ntotal),lwvdef(jcap+1)
c
      offset(n,l)=(jcap+1)*(jcap+2)-(jcap-l+1)*(jcap-l+2)+2*(n-l)
c
      real,allocatable::tmp(:,:)
c
      allocate(tmp(lnt,ntotal))
c
      m=0
      do ll=1,jcap+1
        l=lwvdef(ll)
        do n=l,jcap
          m=m+1
          j=offset(n,l)
          do k=1,ntotal
            tmp(2*m-1,k)=a(j+1,k)
            tmp(2*m  ,k)=a(j+2,k)
          enddo
        enddo
      enddo
c
      do m=1,lnt
        do k=1,ntotal
          a(m,k)=tmp(m,k)
        enddo
      enddo
c
      deallocate(tmp)
c
      return
      end
