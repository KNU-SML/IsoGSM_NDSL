#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
      subroutine mpnn2n(a,lnt2p,b,lln2p,nvar)
c$$$  subprogram documentation block
c
c subprogram:   mpnn2n
c
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract:  transpose (lnt2p,nvar) to (lln2p,nvar)
c
c program history log:
c    99-06-27  henry juang      finish entire test for gsm
c
c usage:        call mpnn2n(a,lnt2p,b,lln2p,nvar)
c
c    input argument lists:
c       a       - real (lnt2p,nvar) sub partial field
c       lnt2p    - integer sub partial spectral grid
c       lln2p    - integer partial spectral grid
c       nvar    - integer total set of fields
c
c    output argument list:
c       b       - real (lln2p,nvar) partial field
c
c subprograms called:
c   mpi_alltoallv  - send and receive from all pe in the same comm
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
      implicit none
#include <paramodel.h>
#include <commpi.h>
      integer lnt2p,lln2p,nvar,status,lntstr0
     1       ,ii,ipe,ipe0,n,m,mn,len,ierr
      real    a,b
      dimension a(lnt2p,nvar),b(lln2p,nvar)
c
      real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::lensnd(:),lenrcv(:)
      integer,allocatable::locsnd(:),locrcv(:)
c
c option for 1-d decomposition
      if( ncol.eq.1 ) then
        do n=1,nvar
          do m=1,lnt2p
            b(m,n)=a(m,n)
          enddo
        enddo
        return
      endif
c
c otherwise
      allocate(tmpsnd(lnt2p*nvar*ncol))
      allocate(tmprcv(lnt2p*nvar*ncol))
      allocate(lensnd(ncol))
      allocate(lenrcv(ncol))
      allocate(locsnd(ncol))
      allocate(locrcv(ncol))
c
      ipe0=int(mype/ncol)*ncol
c
      mn=0
      do ii=1,ncol
        locsnd(ii)=mn
        do n=1,nvar
          do m=1,lntlen(mype)*2
            mn=mn+1
            tmpsnd(mn)=a(m,n)
          enddo
        enddo
        lensnd(ii)=mn-locsnd(ii)
      enddo
c
      mn=0
      do ii=1,ncol
        locrcv(ii)=mn
        ipe=(ii-1)+ipe0
        lenrcv(ii)=lntlen(ipe)*2*nvar
        mn=mn+lenrcv(ii)
      enddo
c
      call mpi_alltoallv(tmpsnd,lensnd,locsnd,MPI_REAL,
     1                   tmprcv,lenrcv,locrcv,MPI_REAL,
     2                   comm_row,ierr)
c received
      mn=0
      do ii=1,ncol
        ipe=(ii-1)+ipe0
        lntstr0=lntstr(ipe)*2-lnpstr(ipe0)*2
        do n=1,nvar
          do m=1,lntlen(ipe)*2
            mn=mn+1
            b(m+lntstr0,n)=tmprcv(mn)
          enddo
        enddo
      enddo
c
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(lensnd)
      deallocate(lenrcv)
      deallocate(locsnd)
      deallocate(locrcv)
c
      return
      end
