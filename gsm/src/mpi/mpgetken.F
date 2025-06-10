#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
#define LEVS levs
#define LPNT lpnt
#define SLVARK slvark
#define MLVARK mlvark
#define LTSTP ltstp
      subroutine mpgetken
c$$$  subprogram documentation block
c
c subprogram: 	mpgetken
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract: transpose (nvrken,npointp,itnum) to (nvrken,npointf,itnum)
c           for ken point output
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpgetken
c
c    input argument lists:
c
c    output argument list:
c 
c subprograms called:
c   mpi_gather	- to gather all message from pe to master pe
c   mpi_gatherv	- same as mpi_gather but variated length
c
c attributes:
c    ligrary: mpi
c    language: fortran 90
c$$$
c
      implicit none
#include <paramodel.h>
#include <commpi.h>
#include <comgpd.h>
#ifdef KEN
      integer nptotal,n,i,j,k,mk,itnum
      real(_mpi_real_), allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::itmppnt(:),itmpgrd(:)
      integer,allocatable::loc(:),len(:)
c
      allocate(tmpsnd(nvrken*nptken*nstken))
      allocate(tmprcv(nvrken*nptken*nstken))
      allocate(itmppnt(nptken))
      allocate(itmpgrd(nptken))
      allocate(loc(0:ncpus_-1))
      allocate(len(0:ncpus_-1))
c
      call mpi_gather(npoint,1,mpi_integer,
     1     itmppnt,1,mpi_integer,0,MPI_COMM_GSM,ierr)
c
      if( mype.eq.master ) then
        mk=0
        do n=0,npes-1
          loc(n)=mk
          len(n)=itmppnt(n)
          mk=mk+len(n)
        enddo
        nptotal=mk
      endif
c
      call mpi_gatherv(igrd,npoint,mpi_integer,
     1     itmpgrd,len(0),loc(0),mpi_integer,0,MPI_COMM_GSM,ierr)
      if( mype.eq.master ) then
        mk=0
        do n=0,npes-1
          do i=1,len(n)
            mk=mk+1
            igrd(mk)=itmpgrd(mk)
          enddo
        enddo
      endif
c
      call mpi_gatherv(jgrd,npoint,mpi_integer,
     1     itmpgrd,len(0),loc(0),mpi_integer,0,MPI_COMM_GSM,ierr)
      if( mype.eq.master ) then
        mk=0
        do n=0,npes-1
          do i=1,len(n)
            mk=mk+1
            jgrd(mk)=itmpgrd(mk)
          enddo
        enddo
      endif
c
      if( mype.eq.master ) then
        mk=0
        do n=0,npes-1
          loc(n)=mk
          do k=1,itnum
            do j=1,itmppnt(n)
              do i=1,nvrken
                mk=mk+1
              enddo
            enddo
          enddo
          len(n)=nvrken*itnum*itmppnt(n)
        enddo
      endif
c
      mk=0
      do k=1,itnum
        do j=1,npoint
          do i=1,nvrken
            mk=mk+1
            tmpsnd(mk)=svdata(i,j,k)
          enddo
        enddo
      enddo
c
      call mpi_gatherv(tmpsnd,mk,MPI_REAL,
     1     tmprcv,len(0),loc(0),MPI_REAL,0,MPI_COMM_GSM,ierr)
c
      if( mype.eq.0 ) then
        mk=0
        do n=0,npes-1
          do k=1,itnum
            do j=1,itmppnt(n)
              do i=1,nvrken
                mk=mk+1
                svdata(i,j,k)=tmprcv(mk)
              enddo
            enddo
          enddo
        enddo
        npoint=nptotal
      endif
c
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(itmppnt)
      deallocate(itmpgrd)
      deallocate(loc)
      deallocate(len)
#endif
c
      return
      end
