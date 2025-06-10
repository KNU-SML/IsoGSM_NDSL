#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
      subroutine mpsp2f(a,lntpp,b,lnt,ntotal)
c$$$  subprogram documentation block
c
c subprogram:   mpsp2f
c
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: transpose (lnpp,kf) to (lnf,kf)
c
c program history log:
c    99-06-27  henry juang      finish entire test for gsm
c
c usage:        call mpsf2p(a,lnt,b,lntpp,ntotal)
c
c    input argument lists:
c       a       - real (lntpp,ntotal) sub partial field
c       lnt     - integer total spectral grid
c       lntpp   - integer sub partial spectral grid
c       ntotal  - integer total set of fields
c
c    output argument list:
c       b       - real (lnt,ntotal) total field
c
c subprograms called:
c   spcshflo	 - shafle spectral for output
c   mpi_gatherv  - gather all pe in the same comm
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
      dimension a(lntpp,ntotal),b(lnt,ntotal)
c
      real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::len(:),loc(:)
c
      allocate(tmpsnd(lntpp*ntotal))
      allocate(tmprcv(lnt*ntotal))
      allocate(len(0:npes_-1))
      allocate(loc(0:npes_-1))
c
c
      if( mype.eq.master ) then
        mk=0
        do n=0,npes-1
          loc(n)=mk
          do m=1,lntlen(n)*2
            do k=1,ntotal
              mk=mk+1
            enddo
          enddo
          len(n)=ntotal*lntlen(n)*2
        enddo
      else
        len(mype)=ntotal*lntlen(mype)*2
      endif
c
      mk=0
      do m=1,lntlen(mype)*2
        do k=1,ntotal
          mk=mk+1
          tmpsnd(mk)=a(m,k)
        enddo
      enddo
c
      call mpi_gatherv(tmpsnd,len(mype),MPI_REAL,
     1     tmprcv,len(0),loc(0),MPI_REAL,
     2     0,MPI_COMM_GSM,ierr)
c
      if( mype.eq.0 ) then
        mk=0
        do n=0,npes-1
          do m=1,lntlen(n)*2
            do k=1,ntotal
              mk=mk+1
              b(m+lntstr(n)*2,k)=tmprcv(mk)
            enddo
          enddo
        enddo
      endif
c
      if( mype.eq.master ) call spcshflo(b,lnt,ntotal,jcap_,lwvdef)
c
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(len)
      deallocate(loc)
c
      return
      end
