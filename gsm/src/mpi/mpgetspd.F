#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
#define LEVS levs
      subroutine mpgetspd(spdmax)
c$$$  subprogram documentation block
c
c subprogram: 	mpgetspd
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract:  get spdmax from all pe
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpgetspd(spdmax)
c
c    input argument lists:
c	spdmax	- real (levs) array of maximal speed for each pe
c
c    output argument list:
c	spdmax	- real (levs) array of maximal speed from all pe for
c                 master pe.
c 
c subprograms called:
c   mpi_allreduce - to gather message from all pe  to master
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
      implicit none
#include <paramodel.h>
#include <commpi.h>
      integer k,ierr
      real    spdmax
      dimension spdmax(levs_)
      real(_mpi_real_),allocatable::spdsnt(:),spdrcv(:)
c
      allocate(spdsnt(levs_))
      allocate(spdrcv(levs_))
c
      do k=1,levs_
        spdsnt(k)=spdmax(k)
      enddo
      call mpi_allreduce(spdsnt,spdrcv,levs_,MPI_REAL,
     1                   mpi_max,MPI_COMM_GSM,ierr)
      do k=1,levs_
        spdmax(k)=spdrcv(k)
      enddo
c
      deallocate(spdrcv)
      deallocate(spdsnt)
c
      return
      end
