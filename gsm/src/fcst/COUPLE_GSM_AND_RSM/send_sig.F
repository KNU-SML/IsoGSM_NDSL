#define JCAP jcap
#define LATG latg
#define LEVH levh
#define LEVP1 levp1
#define LEVS levs
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LONF lonf
#define NCLDG ncldg
#define NPES npes
#define NTOTAL ntotal
#define NTRAC ntrac

#include <define.h>
      subroutine send_sig(q,te,di,ze,rq)
#include <paramodel.h>
#ifdef MP

c$$$  subprogram documentation block
c                .      .    .                                       .
c subprogram:    send_sig   assimble sigma coefficients and send.
c
c usage:    call send_sig(q,te,di,ze,rq)
c   input argument list:
c     q        - ln(psfc)          coefficients.
c     te       - temperature       coefficients.
c     di       - divergence        coefficients.
c     ze       - vorticity         coefficients.
c     rq       - specific humidity coefficients.
c
c$$$
#include <comio.h>
#include <commpi.h>
#define LNT22S lnt22p_
      dimension q(LNT22S),
     1          te(LNT22S,levs_),di(LNT22S,levs_),ze(LNT22S,levs_),
     2          rq(LNT22S,levh_)
c
      parameter(lenrec=lnt2_*(1+levs_*(3+ntotal_)))
      real*_mpi_real_ send_array(lnt2_,1+levs_*(3+ntotal_))
c
      integer len_send_recv
      parameter(len_send_recv=8192)
      integer numchunk
      integer, allocatable :: len(:)
c
      integer ihead(2)
c
      real coef(lnt22_,levs_)
c
      numchunk=(lenrec/len_send_recv)+1
      if(mod(lenrec,len_send_recv).eq.0) then
        numchunk=numchunk-1
      endif
      allocate (len(numchunk))
      do n=1,numchunk
        len(n)=len_send_recv
        if(len_send_recv*n.gt.lenrec) then
            len(n)=lenrec-len_send_recv*(numchunk-1)
        endif
      enddo
c
c  q
c
      call mpsp2f(q,lnt22p_,coef,lnt22_,1)
      if(mype.eq.master) then
        do i=1,lnt2_
          send_array(i,1)=coef(i,1)
        enddo
      endif
c
c tv
c
      call mpsp2f(te,lnt22p_,coef,lnt22_,levs_)
      if(mype.eq.master) then
        do k=1,levs_
          do i=1,lnt2_
            send_array(i,1+k)=coef(i,k)
          enddo
        enddo
      endif
c
c  di
c
      call mpsp2f(di,lnt22p_,coef,lnt22_,levs_)
      if(mype.eq.master) then
        do k=1,levs_
          do i=1,lnt2_
            send_array(i,1+levs_*1+k)=coef(i,k)
          enddo
        enddo
      endif
c
c  ze
c
      call mpsp2f(ze,lnt22p_,coef,lnt22_,levs_)
      if(mype.eq.master) then
        do k=1,levs_
          do i=1,lnt2_
            send_array(i,1+levs_*2+k)=coef(i,k)
          enddo
        enddo
       endif
c
c  rq
c
      do nt = 1,ntotal_
        kt = (nt-1)*levs_+1
        call mpsp2f(rq(1,kt),lnt22p_,coef,lnt22_,levs_)
        if(mype.eq.master) then
          do k=1,levs_
            do i=1,lnt2_
              send_array(i,1+levs_*3+kt-1+k)=coef(i,k)
            enddo
          enddo
        endif
      enddo
c
c  send full sigma array to prog 2
c
      call mpi_barrier(mpi_comm_world,ierr)
      if( ierr.ne.0 ) then
        print *,'PE',mype,':***** Error stop in send_sig ******** '
        print *,'PE',mype,':error code from mpi_barrier = ',ierr
        call mpabort
      endif
c
      if(mype.eq.master) then
        master_ocean=npes_
        call MPI_COMM_RANK(MPI_COMM_WORLD,mype_world,jerr)
#ifdef DBG
        print *,mype_world,' sending sig to ',master_ocean,
     1         ' in chunks of ',len_send_recv,' for ',numchunk,' times'
#endif
        ihead(1)=len_send_recv
        ihead(2)=numchunk
#ifdef DBG
        print *,'len_send_recv,numchunk=',ihead
#endif
        call MPI_SEND(ihead,2,MPI_INTEGER,
     1                master_ocean,msgtag,MPI_COMM_WORLD,ierr)
        call MPI_SEND(len,numchunk,MPI_INTEGER,
     1                master_ocean,msgtag,MPI_COMM_WORLD,ierr)
        do n=1,numchunk
          ind=(n-1)*len_send_recv+1
#ifdef DBG
          print *,'sending sig chunk number=',n,' length=',len(n)
#endif
          call MPI_SEND(send_array(ind,1),len(n),MPI_REAL,
     1                  master_ocean,msgtag,MPI_COMM_WORLD,ierr)
        enddo
c
c  receive something from prog 2
c
c       call MPI_RECV(dummy,20,MPI_REAL,
c    1                master_ocean,0,MPI_COMM_WORLD,status,ierr)
      endif
c
      call mpi_barrier(mpi_comm_world,ierr)
      if( ierr.ne.0 ) then
        print *,'PE',mype,':***** Error stop in send_sig ******** '
        print *,'PE',mype,':error code from mpi_barrier = ',ierr
        call mpabort
      endif
c
      deallocate (len)
c
#endif
      return
      end
