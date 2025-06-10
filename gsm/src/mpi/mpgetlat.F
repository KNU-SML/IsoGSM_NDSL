#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
#define LONF lonf
      subroutine mpgetlat(a,latgp,b,latg,ntotal)
c$$$  subprogram documentation block
c
c subprogram: 	mpgetlat
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: transpose (ntotal,jp) to (ntotal,jf)
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpgetlat(a,latgp,b,latg,ntotal)
c
c    input argument lists:
c	a	- real (ntotal,latgp) any array
c	latgp	- integer partial number of latitudes
c	ntotal	- integer group of latitude
c
c    output argument list:
c	b	- real (ntotal,latg) any array
c	latg	- integer total number of latitudes
c 
c subprograms called:
c   mpi_gatherv	- as gather message with variated length
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
c
      implicit none
#include <paramodel.h>
#include <commpi.h>
      integer latgp,latg,ntotal,
     1        jj,n,j,k,mk,nk,ip,ierr
      real    a,b
      dimension a(ntotal,latgp),b(ntotal,latg)
c
      real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::len(:),loc(:)
c
      allocate(tmpsnd(latgp*ntotal))
      allocate(tmprcv(latg *ntotal))
      allocate(len(0:npes_-1))
      allocate(loc(0:npes_-1))
c
      if( mycol.eq.0 ) then
c
        mk=0
        do k=1,ntotal
          do j=1,latlen(mype)
            mk=mk+1
            tmpsnd(mk)=a(k,j)
            mk=mk+1
            tmpsnd(mk)=a(k,latgp+1-j)
          enddo
        enddo
c
        if( mype.eq.0 ) then
          mk=0
          do n=0,nrow-1
            loc(n)=mk
            ip=mype+ncol*n
            len(n)=ntotal*latlen(ip)*2
            mk=mk+len(n)
          enddo
        else
          len(mype)=ntotal*latlen(mype)*2
        endif
c
        call mpi_gatherv(tmpsnd,len(mype),MPI_REAL,
     1     tmprcv,len(0),loc(0),MPI_REAL,0,comm_column,ierr)
c
        if( mype.eq.0 ) then
          mk=0
          do n=0,nrow-1
            ip=mype+ncol*n
            do k=1,ntotal
              do j=1,latlen(ip)
                jj=j+latstr(ip)-1
                mk=mk+1
                b(k,jj)=tmprcv(mk)
                jj=latg+1-jj
                mk=mk+1
                b(k,jj)=tmprcv(mk)
              enddo
            enddo
          enddo
        endif
c
      endif
c
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(len)
      deallocate(loc)
c
      return
      end
