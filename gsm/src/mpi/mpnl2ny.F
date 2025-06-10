#define NPES npes
#define JCAP jcap
#define LATG latg
#define LATGP latgp
#define LATG2P latg2p
#define JCAP1 jcap1
#define LATG2 latg2
#define LONF lonf
#define JCAPP jcapp
      subroutine mpnl2ny(a,lonfp,latg,b,lonf,latgp,ntotal,kstr,klen)
c$$$  subprogram documentation block
c
c subprogram: 	mpnl2ny
c            
c prgmmr: hann-ming henry juang    org:w/np51   date:99-05-01
c
c abstract:   transpose (ip,kp,jf) to (if,kp,jp)
c
c program history log:
c    99-06-27  henry juang 	finish entire test for gsm
c
c usage:	call mpnl2ny(a,lonfp,latg,b,lonf,latgp,ntotal,kstr,klen)
c
c    input argument lists:
c	a	- real (lonfp,ntotal,latg) partial field in i k
c	lonfp	- integer partial longitude grid x 2
c	latg	- integer total latitude grid / 2
c	lonf	- integer total longitude grid x 2
c	latgp	- integer partial latitude grid / 2
c	kstr	- integer starting set of fields
c	klen	- integer partial vertical layer
c	ntotal	- integer total set of fields
c
c    output argument list:
c	b	- real (lonf,ntotal,latgp) partial field in j k
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
#ifdef REDUCE_GRID
#include <comreduce.h>
#endif
      integer lonf,latgp,ntotal,lonfp,latg,status,lcaps
     1       ,i,i1,i2,j1,j,k,jj,jpe,jpe0,n,mn,len,ierr
     2       ,kstr,klen,ij,ji,llens
      real    a,b
      dimension a(lonfp,ntotal,latg),b(lonf,ntotal,latgp)
c
      real(_mpi_real_),allocatable::tmpsnd(:),tmprcv(:)
      integer,allocatable::lensnd(:),lenrcv(:)
      integer,allocatable::locsnd(:),locrcv(:)
c
      allocate(tmpsnd(lonfp*latgp*klen*nrow))
      allocate(tmprcv(lonfp*latgp*klen*nrow))
      allocate(lensnd(nrow))
      allocate(lenrcv(nrow))
      allocate(locsnd(nrow))
      allocate(locrcv(nrow))
c
      jpe0=mod(mype,ncol)
      lcaps=(jcapp_+1)*2
c cut in y
      mn=0
#ifndef REDUCE_GRID
      llens=lwvlen(mype)*2
#endif
      do jj=1,nrow
        locsnd(jj)=mn
        jpe=jpe0+(jj-1)*ncol
        do n=kstr,kstr+klen-1
          do j=1,latlen(jpe)
            j1=j+latstr(jpe)-1
            ji=latdef(j1)
#ifdef REDUCE_GRID
            llens=lcapdp(ji,mype)*2
#endif
            do i=1,llens
              mn=mn+1
              tmpsnd(mn)=a(i      ,n,ji)
              mn=mn+1
              tmpsnd(mn)=a(i+lcaps,n,ji)
            enddo
          enddo
        enddo
        lensnd(jj)=mn-locsnd(jj)
      enddo
c
      mn=0
      do jj=1,nrow
        locrcv(jj)=mn
        jpe=jpe0+(jj-1)*ncol
#ifdef REDUCE_GRID
        llens=0
        do j=1,latlen(mype)
          j1=j+latstr(mype)-1
          ji=latdef(j1)
          llens=llens+lcapdp(ji,jpe)
        enddo
        lenrcv(jj)=llens*4*klen
#else
        lenrcv(jj)=latlen(mype)*lwvlen(jpe)*4*klen
#endif
        mn=mn+lenrcv(jj)
      enddo
c
      call mpi_alltoallv(tmpsnd,lensnd,locsnd,MPI_REAL,
     1                   tmprcv,lenrcv,locrcv,MPI_REAL,
     2                   comm_column,ierr)
c restore l
      mn=0
      do jj=1,nrow
        jpe=jpe0+(jj-1)*ncol
#ifndef REDUCE_GRID
        llens=lwvlen(jpe)
#endif
        do n=kstr,kstr+klen-1
          do j=1,latlen(mype)
            i1=lwvstr(jpe)
#ifdef REDUCE_GRID
            j1=latdef(latstr(mype)+j-1)
            llens=lcapdp(j1,jpe)
#endif
            do i=1,llens
              ij=lwvdef(i+i1)*2
              mn=mn+1
              b(ij+1      ,n,j)=tmprcv(mn)
              mn=mn+1
              b(ij+1+lonf_,n,j)=tmprcv(mn)
              mn=mn+1
              b(ij+2      ,n,j)=tmprcv(mn)
              mn=mn+1
              b(ij+2+lonf_,n,j)=tmprcv(mn)
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
