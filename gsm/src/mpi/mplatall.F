#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
      subroutine mplatall(a,latgp,b,latg,ntotal)
c$$$  subprogram documentation block
c
c subprogram: 	mplatall
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: transpose (jp,ntotal) to (jf,ntotal)
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mplatall(a,latgp,b,latg,ntotal)
c
c    input argument lists:
c	a	- real (latgp,ntotal) partial field for each pe 
c	latgp	- integer partial grid in latitude divide 2
c	latg	- integer total grid in latitude divide 2
c	ntotal	- integer total set of fields
c
c    output argument list:
c	b	- real (latg,ntotal) full field 
c 
c subprograms called:
c   mpi_scatterv	- scatter message from master pe to all pe
c   mpi_bcast		- boardcast to all
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
      implicit none
#include <paramodel.h>
#include <commpi.h>
      integer latgp,latg,ntotal
     1       ,jj,ji,n,j,k,mk,nk,lens,ip,ierr
      real    a,b
      dimension a(latgp,ntotal),b(latg,ntotal)
c
      real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::len(:),loc(:)
c
      allocate(tmpsnd(latgp*ntotal))
      allocate(tmprcv(latg *ntotal))
      allocate(len(0:nrow-1))
      allocate(loc(0:nrow-1))
c
c
      mk=0
      do k=1,ntotal
        do j=1,latlen(mype)
          mk=mk+1
          tmpsnd(mk)=a(j,k)
        enddo
      enddo
c
      if( myrow.eq.0 ) then
        mk=0
        do n=0,nrow-1
          loc(n)=mk
          ip=mype+ncol*n
          len(n)=ntotal*latlen(ip)
          mk=mk+len(n)
        enddo
      endif
      lens=ntotal*latlen(mype)
c
c --- not working for allgatherv, dont know why ----
c     call mpi_allgatherv(tmpsnd,lens,MPI_REAL,
c    1     tmprcv,len(0),loc(0),MPI_REAL,comm_column,ierr)
c
c --- works for gatherv then bcast ----
      call mpi_gatherv(tmpsnd,lens,MPI_REAL,
     1     tmprcv,len(0),loc(0),MPI_REAL,
     2     0,comm_column,ierr)
      lens=ntotal*latg
      call mpi_bcast(tmprcv,lens,MPI_REAL,
     1     0,comm_column,ierr)
c
      mk=0
      do n=0,nrow-1
        ip=ncol*n
        do k=1,ntotal
          do j=1,latlen(ip)
            jj=j+latstr(ip)-1
            ji=latdef(jj)
            mk=mk+1
            b(ji,k)=tmprcv(mk)
          enddo
        enddo
      enddo
c
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(len)
      deallocate(loc)
c
      return
      end
