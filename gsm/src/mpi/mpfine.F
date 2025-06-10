#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
      subroutine mpfine(endwtime)
c$$$  subprogram documentation block
c
c subprogram: 	mpfine
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: finalizing the mpi by each pe.
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
c   mpi_finalize	- to end of mpi
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
#include <paramodel.h>
#include <commpi.h>
      endwtime=mpi_wtime()
      call mpi_finalize(info)
      if( info.ne.0 ) then
        print *,'PE',mype,': ********* Error stop in mpfine ******* '
        print *,'PE',mype,': error code from mpi_finalize =',info
        call mpabort 
      endif
      return
      end
