#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JCAP jcap
#define JGRD12 jgrd12
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define JGRD jgrd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define NPES npes
      subroutine fixrecio(nsfc,svar,kdim,fld,ioflag)
c
c  fix field record i/o for non-mpi and mpi
c
c  nsfc .... integer unit number of surface file
c  svar .... character identifying the field (for print only)
c  kdim .... integer number of soil layers (or classifications, e.g. albedo)
c  fld  .... real array to write or read (full or partial depending on the MP)
c  ioflag .. =0 read, =1 write
c
      implicit none
#include <paramodel.h>
#ifndef RSM
#define LONF2F lonf2_
#define LATG2F latg2_
#ifdef MP
#include <commpi.h>
#define MPGF2P mpgf2p
#define MPGP2F mpgp2f
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
#else
#define LONF2F igrd12_
#define LATG2F jgrd12_
#ifdef MP
#include <comrmpi.h>
#define MPGF2P rmpgf2p
#define MPGP2F rmpgp2f
#define LONF2S igrd12p_
#define LATG2S jgrd12p_
#else
#define LONF2S igrd12_
#define LATG2S jgrd12_
#endif
#endif
c
#ifdef MP
#ifndef RSM
#define ABORT mpabort
#else
#define ABORT rmpabort
#endif
#else
#define ABORT abort
#endif
c
      integer nsfc,kdim,ioflag
      real fld(LONF2S,LATG2S,kdim)
      character*8 svar  
c
      integer ij,k
c
#ifdef MP
      real, allocatable :: grid(:,:,:)
#endif
#ifdef MP
      allocate (grid(LONF2F,LATG2F,kdim))
#endif
c
      if(ioflag.eq.1) go to 1000
c
c read record
c
#ifdef MP
      if( mype.eq.master ) then
        read(nsfc,err=101) 
     1        ((grid(ij,1,k),ij=1,LONF2F*LATG2F),k=1,kdim)
        go to 100
  101   continue
        print *,'sfc file read error for var=',svar
        call ABORT
  100   continue
#ifdef DBG
        call rmaxmin(grid,LONF2F*LATG2F,kdim,1,1,svar)
#endif
#ifndef RSM
        do k = 1, kdim
          call row1ns(grid(1,1,k))
        enddo
#else
        call shalfi(grid,kdim)
#endif
      endif
#ifdef RSM
      call rmpsynall
#else
      call mpsynall
#endif
      call MPGF2P(grid,LONF2F,LATG2F,fld,LONF2S,LATG2S,kdim)
#else
      read(nsfc,err=111) 
     1        ((fld(ij,1,k),ij=1,LONF2F*LATG2F),k=1,kdim)
      go to 110
  111 continue
      print *,'sfc file read error for var=',svar
      call ABORT
  110 continue
#ifdef DBG
      do k=1,kdim
        call rmaxmin(fld(1,1,k),LONF2F*LATG2F,1,1,1,svar)
      enddo
#endif
#ifndef RSM
      do k = 1, kdim
        call row1ns(fld(1,1,k))
      enddo
#else
      call shalfi(fld,kdim)
#endif
#ifndef NOPRINT
      print *,' fixrdrec completed for ',svar
#endif
#endif
#ifdef MP
      deallocate (grid)
#endif
      return
c--------------------------------------------------------------------
c
c write field
c
 1000 continue
#ifdef MP
#ifdef RSM
      call rmpsynall
#else
      call mpsynall
#endif
      call MPGP2F(fld,LONF2S,LATG2S,grid,LONF2F,LATG2F,kdim)
      if( mype.eq.master ) then
#ifndef RSM
        do k = 1, kdim
          call rowsep(grid(1,1,k))
        enddo
#else
        call shalfo(grid,kdim)
#endif
        write(nsfc) ((grid(ij,1,k),ij=1,LONF2F*LATG2F),k=1,kdim)
      endif
#else
#ifndef RSM
      do k = 1, kdim
        call rowsep(fld(1,1,k))
      enddo
#else
      call shalfo(fld,kdim)
#endif
      write(nsfc) ((fld(ij,1,k),ij=1,LONF2F*LATG2F),k=1,kdim)
#ifndef RSM
      do k = 1, kdim
        call row1ns(fld(1,1,k))
      enddo
#else
      call shalfi(fld,kdim)
#endif
#endif
#ifdef DBG
#ifdef MP
      if( mype.eq.master ) then
#endif
      print *,' fixrecio completed for ',svar
      do k=1,kdim
        call rmaxmin(fld(1,1,k),LONF2S*LATG2S,1,1,1,svar)
      enddo
#ifdef MP
      endif
#endif
#endif
#ifdef MP
      deallocate (grid)
#endif
      return
      end
