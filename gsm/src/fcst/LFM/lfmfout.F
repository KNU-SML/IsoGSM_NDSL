#define JCAP jcap
#define LALBD lalbd
#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVH levh
#define LEVP1 levp1
#define LEVS levs
#define LNT22 lnt22
#define LNT22P lnt22p
#define LNT2 lnt2
#define LNTP lntp
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp
#define LSOIL lsoil
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine lfmfout(ifstep,thour,nz)                                       
c                                                                               
#include <paramodel.h>
#include <comio.h>
#include <comlfm.h>
c                                                                               
      dimension si(levp1_),sl(levs_),idate(4)                                   
c                                                                               
      dimension gz(lnt2_)                                                       
c     dimension fslmsk(lonf2_,latg2_)                                           
c     dimension fcvb(lonf2_,latg2_),fcvt(lonf2_,latg2_)                         
      dimension work(lonf2_*latg2_)                                             
c                                                                               
      rewind nz                                                                 
      read(nz)                                                                  
      read(nz) dummy,idate,(si(k),k=1,levp1_),(sl(k),k=1,levs_ )                
c                                                                               
#ifndef NOPRINT
      write(6,*) 'lfmfout ifstep=',ifstep,' klenp=',klenp                       
#endif
      if(ifstep.eq.klenp) then                                                  
        vhour=nint(thour-filtwin*0.5)                                           
#ifndef NOPRINT
        write(6,*) 'vhour,idate of filtered output sig=',vhour,idate            
#endif
        write(nlfmsgo) lab                                                      
        write(nlfmsgo) vhour,idate,(si(k),k=1,levp1_),(sl(k),k=1,levs_)         
        read(nz)(gz(i),i=1,lnt2_ )                                              
        write(nlfmsgo)(gz(i),i=1,lnt2_ )                                        
        write(nlfmsgo)(fq(i),i=1,lnt2_ )                                        
        do k=1,levs_                                                            
          write(nlfmsgo) (fte(i,k),i=1,lnt2_ )                                  
        enddo                                                                   
        do k=1,levs_                                                            
          write(nlfmsgo) (fdi(i,k),i=1,lnt2_ )                                  
          write(nlfmsgo) (fze(i,k),i=1,lnt2_ )                                  
        enddo                                                                   
        do k=1,levs_                                                            
          write(nlfmsgo) (frq(i,k),i=1,lnt2_ )                                  
        enddo                                                                   
c                                                                               
        write(nlfmsfo) lab                                                      
        write(nlfmsfo) vhour,idate                                              
        call rowsep(ftsea)                                                      
        write(nlfmsfo) ftsea                                                    
        do k = 1,lsoil_                                                         
          do j = 1, latg2_                                                      
            do i = 1, lonf2_                                                    
              ij = (j-1) * lonf2_  + i                                          
              work(ij) = fsmc(i,j,k)                                            
            enddo                                                               
          enddo                                                                 
          call rowsep(work)                                                     
          do j = 1, latg2_                                                      
            do i = 1, lonf2_                                                    
              ij = (j-1) * lonf2_  + i                                          
              fsmc(i,j,k) = work(ij)                                            
            enddo                                                               
          enddo                                                                 
        enddo                                                                   
        write(nlfmsfo) fsmc                                                     
        call rowsep(fsheleg)                                                    
        write(nlfmsfo) fsheleg                                                  
        do k = 1,lsoil_                                                         
          do j = 1, latg2_                                                      
            do i = 1, lonf2_                                                    
              ij = (j-1) * lonf2_  + i                                          
              work(ij) = fstc(i,j,k)                                            
            enddo                                                               
          enddo                                                                 
          call rowsep(work)                                                     
          do j = 1, latg2_                                                      
            do i = 1, lonf2_                                                    
              ij = (j-1) * lonf2_  + i                                          
              fstc(i,j,k) = work(ij)                                            
            enddo                                                               
          enddo                                                                 
        enddo                                                                   
        write(nlfmsfo) fstc                                                     
        call rowsep(ftg3)                                                       
        write(nlfmsfo) ftg3                                                     
        call rowsep(fzorl)                                                      
        write(nlfmsfo) fzorl                                                    
        call rowsep(fcv)                                                        
        write(nlfmsfo) fcv                                                      
        do j=1,latg2_                                                           
          do i=1,lonf2_                                                         
            if(wcvb(i,j).ge.0.5) then                                           
              fcvb(i,j)=fcvb(i,j)/wcvb(i,j)                                     
            endif                                                               
          enddo                                                                 
        enddo                                                                   
        call rowsep(fcvb)                                                       
        write(nlfmsfo) fcvb                                                     
        do j=1,latg2_                                                           
          do i=1,lonf2_                                                         
            if(wcvt(i,j).ge.0.5) then                                           
              fcvt(i,j)=fcvt(i,j)/wcvt(i,j)                                     
            endif                                                               
          enddo                                                                 
        enddo                                                                   
        call rowsep(fcvt)                                                       
        write(nlfmsfo) fcvt                                                     
        call rowsep(falbedo)                                                    
        write(nlfmsfo) falbedo                                                  
        do j=1,latg2_                                                           
          do i=1,lonf2_                                                         
            if(islmsk(i,j,3).ge.klenp/2) then                                   
              fslmsk(i,j)=2.                                                    
            elseif(islmsk(i,j,1).gt.islmsk(i,j,2)) then                         
              fslmsk(i,j)=0.                                                    
            else                                                                
              fslmsk(i,j)=1.                                                    
            endif                                                               
          enddo                                                                 
        enddo                                                                   
        call rowsep(fslmsk)                                                     
        write(nlfmsfo) fslmsk                                                   
        call rowsep(fplantr)                                                    
        write(nlfmsfo) fplantr                                                  
        call rowsep(fcanopy)                                                    
        write(nlfmsfo) fcanopy                                                  
        call rowsep(ff10m)                                                      
        write(nlfmsfo) ff10m                                                    
      else                                                                      
        write(nlfmsgo) lab                                                      
        write(nlfmsgo) ifstep,idate,(si(k),k=1,levp1_),(sl(k),k=1,levs_)        
#ifndef NOPRINT
        write(6,*) 'ifstep,idate of filtered sig=',ifstep,idate                 
#endif
        write(nlfmsgo)(gz(i),i=1,lnt2_ )                                        
        write(nlfmsgo)(fq(i),i=1,lnt2_ )                                        
        do k=1,levs_                                                            
          write(nlfmsgo) (fte(i,k),i=1,lnt2_ )                                  
        enddo                                                                   
        do k=1,levs_                                                            
          write(nlfmsgo) (fdi(i,k),i=1,lnt2_ )                                  
          write(nlfmsgo) (fze(i,k),i=1,lnt2_ )                                  
        enddo                                                                   
        do k=1,levs_                                                            
          write(nlfmsgo) (frq(i,k),i=1,lnt2_ )                                  
        enddo                                                                   
c                                                                               
c  surface file                                                                 
c                                                                               
        write(nlfmsfo) lab                                                      
        write(nlfmsfo) ifstep,idate                                             
#ifndef NOPRINT
        write(6,*) 'ifstep,idate of filtered sfc=',ifstep,idate                 
#endif
        write(nlfmsfo) ftsea                                                    
        write(nlfmsfo) fsmc                                                     
        write(nlfmsfo) fsheleg                                                  
        write(nlfmsfo) fstc                                                     
        write(nlfmsfo) ftg3                                                     
        write(nlfmsfo) fzorl                                                    
        write(nlfmsfo) fcv                                                      
        write(nlfmsfo) fcvb,wcvb                                                
        write(nlfmsfo) fcvt,wcvt                                                
        write(nlfmsfo) falbedo                                                  
        write(nlfmsfo) islmsk                                                   
        write(nlfmsfo) fplantr                                                  
        write(nlfmsfo) fcanopy                                                  
        write(nlfmsfo) ff10m                                                    
      endif                                                                     
      return                                                                    
      end                                                                       
