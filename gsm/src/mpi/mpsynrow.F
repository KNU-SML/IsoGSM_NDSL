#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
      subroutine mpsynrow
c$$$  subprogram documentation block
c
c subprogram: 	mpsynrow
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: set mpi barrier for comm_row nodes
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpsynrow
c
c subprograms called:
c   mpi_barrier  - set barrier for comm
c
c attributes:
c    library: mpi
c    language: fortran 90
#include <paramodel.h>
#include <commpi.h>
      call mpi_barrier(comm_row,ierr)
      if( ierr.ne.0 ) then
        print *,'PE',mype,':***** Error stop in mpsynrow ******** '
        print *,'PE',mype,':error code from mpi_barrier = ',ierr
        call mpabort
      endif
      return
      end
