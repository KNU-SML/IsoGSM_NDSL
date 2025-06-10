#define CP cp
#define RD rd
#define TTP ttp
      subroutine gthe                                                           
#ifdef CRAY_THREAD
cfpp$ expand(fthex)                                                             
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: gthe        compute equivalent potential temperature table        
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute equivalent potential temperature table                      
c   as a function of lcl temperature and pressure over 100 kpa                  
c   to the kappa power for function fthe.                                       
c   equivalent potential temperatures are calculated in subprogram fthex        
c   the current implementation computes a table with a first dimension          
c   of 241 for temperatures ranging from 183.16 to 303.16 kelvin                
c   and a second dimension of 151 for pressure over 100 kpa                     
c   to the kappa power ranging from 0.04**rocp to 1.10**rocp.                   
c                                                                               
c program history log:                                                          
c   91-05-07  iredell                                                           
c   94-12-30  iredell             expand table                                  
c                                                                               
c usage:  call gthe                                                             
c                                                                               
c subprograms called:                                                           
c   (fthex)  - inlinable function to compute equiv. pot. temperature            
c                                                                               
c common blocks:                                                                
c   comthe   - scaling parameters and table for function fthe.                  
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      parameter(cp=cp_,rd=rd_,ttp=ttp_)                                         
      parameter(rocp=rd/cp)                                                     
      parameter(nx=241,ny=151)                                                  
      dimension tbthe(nx,ny)                                                    
      common/comthe/ c1xthe,c2xthe,c1ythe,c2ythe,tbthe                          
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      xmin=ttp-90.                                                              
      xmax=ttp+30.                                                              
      xinc=(xmax-xmin)/(nx-1)                                                   
      c1xthe=1.-xmin/xinc                                                       
      c2xthe=1./xinc                                                            
      ymin=0.04**rocp                                                           
chmhj                                                                           
cmk   ymin=0.05**rocp                                                           
      ymax=1.10**rocp                                                           
      yinc=(ymax-ymin)/(ny-1)                                                   
      c1ythe=1.-ymin/yinc                                                       
      c2ythe=1./yinc                                                            
      do jy=1,ny                                                                
        y=ymin+(jy-1)*yinc                                                      
        pk=y                                                                    
        do jx=1,nx                                                              
          x=xmin+(jx-1)*xinc                                                    
          t=x                                                                   
          tbthe(jx,jy)=fthex(t,pk)                                              
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
