#define JCAP jcap
#define LATG latg
#define LATG2 latg2
#define LATG2P latg2p
#define LATGP latgp
#define LONF lonf
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONFP lonfp
#define NPES npes
#include <define.h>
      subroutine insur2(nflin,alvsf,alnsf,alvwf,alnwf,facsf,facwf
#ifndef SW99UPDATE
     1                  ,paerf)
#else
     1                  )
#endif
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c**************************************************************                 
c  the code reads in surface albedo and aerosol data. the albedo                
c  data are derived from matthews vegetation index by using a                  
c  modified brieglebs scheme.  the aerosol distribution data is                
c  based on matthews vegetation index.    --- y.hou  mar 7, 1995                
c**************************************************************                 
#include <paramodel.h>
#ifdef MP
#include <commpi.h>
      real,allocatable::grid(:,:)
#endif
#ifdef MP
#define LONF2S lonf2p_
#define LATG2S latg2p_
#else
#define LONF2S lonf2_
#define LATG2S latg2_
#endif
      parameter (ix=lonf2_,  iy=latg2_)                                         
                       d i m e n s i o n                                        
     1  alvsf(LONF2S,LATG2S,4),alnsf(LONF2S,LATG2S,4)
     2, alvwf(LONF2S,LATG2S,4),alnwf(LONF2S,LATG2S,4)
     3, facsf(LONF2S,LATG2S),  facwf(LONF2S,LATG2S)
#ifndef SW99UPDATE
     4, paerf(LONF2S,LATG2S,5)
#endif
c                                                                               
#ifdef MP
      if(mype.eq.master) then
#endif
         open(unit=nflin,file='albaer.snl',form='unformatted'
     1    ,status='old')
#ifdef MP
      endif
#endif
c
#ifdef MP
      allocate(grid(ix,iy))
      if(mype.eq.master) rewind nflin
      do k=1,4                                                                  
        if(mype.eq.master) read(nflin) ((grid(i,j),i=1,ix),j=1,iy)
        call row1ns  (grid)                                             
        call mpgf2p(grid,lonf2_,latg2_,alvsf(1,1,k),lonf2p_,latg2p_,1)
      enddo                                                                     
      do k=1,4                                                                  
        if(mype.eq.master) read(nflin) ((grid(i,j),i=1,ix),j=1,iy)
        call row1ns  (grid)                                             
        call mpgf2p(grid,lonf2_,latg2_,alvwf(1,1,k),lonf2p_,latg2p_,1)
      enddo                                                                     
      do k=1,4                                                                  
        if(mype.eq.master) read(nflin) ((grid(i,j),i=1,ix),j=1,iy)
        call row1ns  (grid)                                             
        call mpgf2p(grid,lonf2_,latg2_,alnsf(1,1,k),lonf2p_,latg2p_,1)
      enddo                                                                     
      do k=1,4                                                                  
        if(mype.eq.master) read(nflin) ((grid(i,j),i=1,ix),j=1,iy)
        call row1ns  (grid)                                             
        call mpgf2p(grid,lonf2_,latg2_,alnwf(1,1,k),lonf2p_,latg2p_,1)
      enddo                                                                     
      if(mype.eq.master) read(nflin) grid
      call row1ns(grid)                                                        
      call mpgf2p(grid,lonf2_,latg2_,facsf,lonf2p_,latg2p_,1)
      if(mype.eq.master) read(nflin) grid
      call row1ns(grid)                                                        
      call mpgf2p(grid,lonf2_,latg2_,facwf,lonf2p_,latg2p_,1)
#ifndef SW99UPDATE
      do 10 k=1,5                                                               
        if(mype.eq.master) read(nflin) ((grid(i,j),i=1,ix),j=1,iy)
        call row1ns(grid)                                               
        call mpgf2p(grid,lonf2_,latg2_,paerf(1,1,k),lonf2p_,latg2p_,1)
  10  continue                                                                  
#endif
c ....................................................................
#else
      rewind nflin                                                              
      do k=1,4                                                                  
        read(nflin) ((alvsf(i,j,k),i=1,ix),j=1,iy)                              
        call row1ns  (alvsf(1,1,k))                                             
      enddo                                                                     
      do k=1,4                                                                  
        read(nflin) ((alvwf(i,j,k),i=1,ix),j=1,iy)                              
        call row1ns  (alvwf(1,1,k))                                             
      enddo                                                                     
      do k=1,4                                                                  
        read(nflin) ((alnsf(i,j,k),i=1,ix),j=1,iy)                              
        call row1ns  (alnsf(1,1,k))                                             
      enddo                                                                     
      do k=1,4                                                                  
        read(nflin) ((alnwf(i,j,k),i=1,ix),j=1,iy)                              
        call row1ns  (alnwf(1,1,k))                                             
      enddo                                                                     
      read(nflin) facsf                                                         
      call row1ns(facsf)                                                        
      read(nflin) facwf                                                         
      call row1ns(facwf)                                                        
#ifndef SW99UPDATE
      do 10 k=1,5                                                               
        read(nflin) ((paerf(i,j,k),i=1,ix),j=1,iy)                              
        call row1ns(paerf(1,1,k))                                               
  10  continue                                                                  
#endif
#endif
c                                                                               
      return                                                                    
      end                                                                       
