#define NPES npes
#define JCAP jcap
#define LATG latg
      subroutine mpabort
c$$$  subprogram documentation block
c
c subprogram: 	mpabort
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: abort the mpi by any pe.
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpfine
c
c    input argument lists:
c
c    output argument list:
c 
c subprograms called:
c   mpi_abort	- to abort mpi
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
#include <paramodel.h>
#include <commpi.h>
      call mpi_abort(MPI_COMM_GSM,1,info)
      if( info.ne.0 ) then
        print *,'PE',mype,': ********* Error stop in mpabort ****** '
        print *,'PE',mype,': error code from mpi_abort =',info
        call abort 
      endif
      return
      end
