#define NPES  npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
      subroutine mpbcastl(n,len)
c$$$  subprogram documentation block
c
c subprogram: 	mpbcasti
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract:  broadcast integer array to all pes
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	mpbcasti(n,len)
c
c    input argument lists:
c	n	- integer (len) array from master pe 
c	len	- integer length of array n
c
c    output argument list:
c	n	- integer (len) array from master pe 
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
      logical n
      integer len
      dimension n(len)
c
      call mpi_bcast(n,len,mpi_logical,0,MPI_COMM_GSM,ierr)
c
      return
      end
