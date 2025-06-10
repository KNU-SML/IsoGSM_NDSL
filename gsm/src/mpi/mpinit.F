#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
#define NCOL ncol
#define NROW nrow
      subroutine mpinit(strwtime)
c$$$  subprogram documentation block
c
c subprogram: 	mpinit
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: initializing the mpi by each pe.
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpinit
c
c subprograms called:
c   mpi_init		- to initial mpi call
c   mpi_comm_size	- to get total number of pe for the comm
c   mpi_comm_rank	- to get my own pe number 
c   mpi_comm_split	- to split entire comm into sub comm
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
#include <paramodel.h>
#include <commpi.h>
      integer color
      call mpi_init(info)
      strwtime=mpi_wtime()
      call mpi_comm_size(MPI_COMM_WORLD,npes,ierr)
      call mpi_comm_rank(MPI_COMM_WORLD,mype,jerr)
      if(mype.eq.master) then
#ifdef COUPLE
        print *,' This is coupled run.  Total number of pe is ',npes
#else
        print *,' This is uncoupled.  Nmber of pe is ',npes
#endif
      endif
      color=0
      call mpi_comm_split(MPI_COMM_WORLD,color,mype,mpi_comm_gsm,ierr)
      call mpi_comm_size(mpi_comm_gsm,npes,ierr)
      call mpi_comm_rank(mpi_comm_gsm,mype,jerr)
      if( ierr.ne.0 .or. jerr.ne.0 .or.
     1    info.ne.0 .or. npes.ne.npes_ ) then
        print *,'PE',mype,':********* Error stop in rmpinit ********* '
        print *,'PE',mype,':error code from mpi_init = ',info
        print *,'PE',mype,':error code from mpi_comm_size = ',ierr
        print *,'PE',mype,':error code from mpi_comm_rank = ',jerr
        print *,'PE',mype,':npes preset in code = ',npes_
        print *,'PE',mype,':npes request from environment = ',npes
        print *,'PE',mype,':******* End of output for mpinit ******* '
        call mpabort
      else
        master=0
        ncol=ncol_
        nrow=nrow_
        msgtag=0
        myrow=mype/ncol
        mycol=mod(mype,ncol)
        call mpi_comm_split(mpi_comm_gsm,myrow,mycol,comm_row,ierr)
        call mpi_comm_split(mpi_comm_gsm,mycol,myrow,comm_column,jerr)
        if( ierr.ne.0 .or. jerr.ne.0 ) then
          print *,'PE',mype,':********* Error stop in mpinit ********* '
          print *,'PE',mype,':error code for doing comm_row = ',ierr
          print *,'PE',mype,':error code for doing comm_column = ',jerr
          print *,'PE',mype,':******* End of output for mpinit ******* '
          call mpabort
        endif
      endif
      return
      end
