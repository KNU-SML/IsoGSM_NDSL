#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
      subroutine mpbcastr(a,len)
c$$$  subprogram documentation block
c
c subprogram: 	mpbcastr
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: broadcast real array to all pes
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	mpbcastr(a,len)
c
c    input argument lists:
c	a	- real (len) array from master pe 
c	len	- integer length of array a
c
c    output argument list:
c	a	- real (len) array from master pe 
c 
c subprograms called:
c   mpi_bcast	- to broadcast to all pe in the comm
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
#include <paramodel.h>
#include <commpi.h>
      integer len,i
      real a
      dimension a(len)
c
      real(_mpi_real_),allocatable::tmpa(:)
c
      allocate(tmpa(len))
c
      if( mype.eq.master ) then
        do i=1,len
          tmpa(i) = a(i)
        enddo
      endif
      call mpi_bcast(tmpa,len,MPI_REAL,
     1               0,MPI_COMM_GSM,ierr)
      do i=1,len
        a(i) = tmpa(i)
      enddo
      deallocate(tmpa)
c
      return
      end
