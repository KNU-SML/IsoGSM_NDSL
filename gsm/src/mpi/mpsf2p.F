#define NPES  npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
      subroutine mpsf2p(a,lnt,b,lntpp,ntotal)
c$$$  subprogram documentation block
c
c subprogram: 	mpsf2p
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: transpose (lnf,kf) to (lnpp,kf)
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpsf2p(a,lnt,b,lntpp,ntotal)
c
c    input argument lists:
c	a	- real (lnt,ntotal) total field
c	lnt	- integer total spectral grid
c	lntpp	- integer sub partial spectral grid
c	ntotal	- integer total set of fields
c
c    output argument list:
c	b	- real (lntpp,ntotal) sub partial field
c 
c subprograms called:
c   spcshfli		- shafle spectral for balance from input
c   spcshflo		- shafle spectral for output
c   mpi_scatterv  - send to all pe in the same comm
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
      implicit none
#include <paramodel.h>
#include <commpi.h>
      integer lnt,lntpp,ntotal
     1       ,n,m,k,mk,ierr
      real    a,b
      dimension a(lnt,ntotal),b(lntpp,ntotal)
c
      real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::len(:),loc(:)
c
      allocate(tmpsnd(lnt*ntotal))
      allocate(tmprcv(lntpp*ntotal))
      allocate(len(0:npes_-1))
      allocate(loc(0:npes_-1))
c
c
      if( mype.eq.master ) then
        call spcshfli(a,lnt,ntotal,jcap_,lwvdef)
        mk=0
        do n=0,npes-1
          loc(n)=mk
          do m=1,lntlen(n)*2
            do k=1,ntotal
              mk=mk+1
              tmpsnd(mk)=a(m+lntstr(n)*2,k)
            enddo
          enddo
          len(n)=ntotal*lntlen(n)*2
        enddo
      else
        len(mype)=ntotal*lntlen(mype)*2
      endif
c
      call mpi_scatterv(tmpsnd,len(0),loc(0),MPI_REAL,
     1     tmprcv,len(mype),MPI_REAL,0,MPI_COMM_GSM,ierr)
c
      mk=0
      do m=1,lntlen(mype)*2
        do k=1,ntotal
          mk=mk+1
          b(m,k)=tmprcv(mk)
        enddo
      enddo
c
      if( mype.eq.0 ) call spcshflo(a,lnt,ntotal,jcap_,lwvdef)
c
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(len)
      deallocate(loc)
c
      return
      end
