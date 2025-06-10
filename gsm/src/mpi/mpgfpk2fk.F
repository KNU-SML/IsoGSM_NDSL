#define NPES npes
#define JCAP jcap
#define LATG latg
#define JCAP1 jcap1
#define LATG2 latg2
#define LATGP latgp
#define LATG2P latg2p
#define LONF lonf
      subroutine mpgfpk2fk(a,idim,nlen,nstr,nend,
     1                     b,ntotal)
c$$$  subprogram documentation block
c
c subprogram: 	mpgfpk2fk
c            
c prgmmr: Hideki Kanamaru and Masao Kanamitsu org:ECPC/SIO   date:04-08-27
c
c abstract:  gather pgb records to master processor
c
c program history log:
c    04-08-27  H.Kanamaru and M.Kanamitsu  the first version
c
c usage:	call mpgfpk2fk(a,idim,nlen,nstr,nend,b,ntotal)
c
c    input argument lists:
c	a	- 
c	idim	- 
c	nlen	- 
c	nstr	- 
c	nend	- 
c	ntotal	- 
c
c    output argument list:
c	b	- 
c	nlen	- 
c 
c subprograms called:
c   mpi_gatherv	- gather message from all pe to master
c   mpi_allgatherv - gather message from all pe to all pe
c
c attributes:
c    library: mpi
c    language: fortran 90
c$$$
      implicit none
#include <paramodel.h>
#ifndef RSM
#include <commpi.h>
#else
#include <comrmpi.h>
#endif
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
      integer idim,nstr,nend,ntotal,
     1        n,i,j,k,mk,kk,ierr,nstrend,nlocac,lenac
      character a(idim,ntotal),b(idim,ntotal)
      dimension nstr(0:npes_-1),nend(0:npes_-1),nstrend(0:npes_-1)
      integer nlen(ntotal)
c
      character, allocatable:: tmpsnd(:)
      character, allocatable:: tmprcv(:)
c
      integer,allocatable::len(:),loc(:),mlen(:),nloc(:)
c
      allocate(tmpsnd(idim*(nend(mype)-nstr(mype)+1)))
      allocate(tmprcv(idim*ntotal))
      allocate(len(npes_))
      allocate(loc(0:npes_-1))
      allocate(mlen(ntotal))
      allocate(nloc(0:npes_-1))
c
c make nstrend from nstr and nend, number of record fileds in each pe
c also make nloc, accumulated number of record files from pe=0 
c
      nlocac=0 
      do i=0,npes_-1
         nloc(i)=nlocac
         nstrend(i)=nend(i)-nstr(i)+1
         nlocac=nlocac+nstrend(i)
      enddo
c
c gather and bcst nlen into mlen, array count for each record
c
      call mpi_allgatherv(nlen,nstrend(mype),mpi_integer,
     1    mlen,nstrend(0),nloc(0),mpi_integer,MPI_COMM_GSM,ierr)
c
c rearrange a into tmpsnd
c
      mk=0
      do k=nstr(mype),nend(mype)
        do i=1,nlen(k-nstr(mype)+1)
          mk=mk+1
          tmpsnd(mk)=a(i,k-nstr(mype)+1)
        enddo
      enddo
c
c make len from mlen, array count for each pe
c rearrange len into loc, for use in gatherv of tmpsnd
c
      do i=1,npes
         len(i)=0
      enddo
      lenac=0
      do j=0,npes-1
         do i=nloc(j)+1,nloc(j)+nstrend(j)
            len(j+1)=len(j+1)+mlen(i)
         enddo
         loc(j)=lenac
         lenac=lenac+len(j+1)
      enddo
c
c gather to master, tmpsnd to tmprcv
c
      call mpi_gatherv(tmpsnd,len(mype+1),mpi_character,
     1     tmprcv,len,loc(0),mpi_character,
     2     0,MPI_COMM_GSM,ierr)
c
c rearrange tmprcv into b
c
      if( mype.eq.0 ) then
        mk=0
          do k=1,ntotal
            do i=1,mlen(k)
               mk=mk+1
               b(i,k)=tmprcv(mk)
            enddo
          enddo
      endif
c
c pass mlen to nlen
c
      do i=1,ntotal
         nlen(i)=mlen(i)
      enddo
c
      deallocate(tmpsnd)
      deallocate(tmprcv)
      deallocate(len)
      deallocate(loc)
      deallocate(mlen)
      deallocate(nloc)
c
      return
      end
