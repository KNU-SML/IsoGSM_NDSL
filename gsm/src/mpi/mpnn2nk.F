#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
      subroutine mpnn2nk(a,lnt2p,levs,b,lln2p,levsp,nvar)
c$$$  subprogram documentation block
c
c subprogram:   mpnn2nk
c
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract:  transpose (lnt2p,levs,nvar) to (lln2p,levsp,nvar)
c
c program history log:
c    99-06-27  henry juang      finish entire test for gsm
c
c usage:        call mpnn2nk(a,lnt2p,levs,b,lln2p,levsp,nvar)
c
c    input argument lists:
c       a       - real (lnt2p,levs,nvar) partial field
c       lnt2p    - integer sub partial spectral grid
c       levsp   - integer partial vertical grid
c       lln2p   - integer partial spectral grid
c       levs    - integer full vertical grid
c       nvar    - integer total set of fields
c
c    output argument list:
c       b       - real (lln2p,levsp,nvar) partial field
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
      integer lnt2p,levs,lln2p,levsp,nvar,status,lntstr0
     1       ,ii,ipe,ipe0,n,k,m,mn,len,ierr
      real    a,b
      dimension a(lnt2p,levs,nvar),b(lln2p,levsp,nvar)
c
      real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::lensnd(:),lenrcv(:)
      integer,allocatable::locsnd(:),locrcv(:)
c
c option for 1-d decomposition
      if( ncol.eq.1 ) then
        do n=1,nvar
          do k=1,levs
            do m=1,lnt2p
              b(m,k,n)=a(m,k,n)
            enddo
          enddo
        enddo
        return
      endif
c
c otherwise
      allocate(tmpsnd(lnt2p*levsp*nvar*ncol))
      allocate(tmprcv(lnt2p*levsp*nvar*ncol))
      allocate(lensnd(ncol))
      allocate(lenrcv(ncol))
      allocate(locsnd(ncol))
      allocate(locrcv(ncol))
c
c
      ipe0=int(mype/ncol)*ncol
c cuting  k
      mn=0
      do ii=1,ncol
        locsnd(ii)=mn
        ipe=(ii-1)+ipe0
        do n=1,nvar
          do k=1,levlen(ipe)
            do m=1,lntlen(mype)*2
              mn=mn+1
              tmpsnd(mn)=a(m,k+levstr(ipe)-1,n)
            enddo
          enddo
        enddo
        lensnd(ii)=mn-locsnd(ii)
      enddo
c
      mn=0
      do ii=1,ncol
        locrcv(ii)=mn
        ipe=(ii-1)+ipe0
        lenrcv(ii)=levlen(mype)*lntlen(ipe)*2*nvar
        mn=mn+lenrcv(ii)
      enddo
c
      call mpi_alltoallv(tmpsnd,lensnd,locsnd,MPI_REAL,
     1                   tmprcv,lenrcv,locrcv,MPI_REAL,
     2                   comm_row,ierr)
      if(ierr.ne.0) then
        print *,' error after alltoallv ',ierr
      endif
c restore  n
      mn=0
      do ii=1,ncol
        ipe=(ii-1)+ipe0
        lntstr0=lntstr(ipe)*2-lnpstr(ipe0)*2
        do n=1,nvar
          do k=1,levlen(mype)
            do m=1,lntlen(ipe)*2
              mn=mn+1
              b(m+lntstr0,k,n)=tmprcv(mn)
            enddo
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
