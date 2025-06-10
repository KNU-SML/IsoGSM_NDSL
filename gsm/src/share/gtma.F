#define CP cp
#define RD rd
      subroutine gtma                                                           
#ifdef CRAY_THREAD
cfpp$ expand(ftmaxg)                                                            
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: gtma         compute moist adiabat tables                         
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute temperature and specific humidity tables                    
c   as a function of equivalent potential temperature and                       
c   pressure over 100 kpa to the kappa power for function ftma.                 
c   exact parcel temperatures are calculated in subprogram ftmaxg.              
c   the current implementation computes a table with a first dimension          
c   of 151 for equivalent potential temperatures ranging from 200 to 500        
c   kelvin and a second dimension of 121 for pressure over 100 kpa              
c   to the kappa power ranging from 0.01**rocp to 1.10**rocp.                   
c                                                                               
c program history log:                                                          
c   91-05-07  iredell                                                           
c   94-12-30  iredell             expand table                                  
c                                                                               
c usage:  call gtma                                                             
c                                                                               
c subprograms called:                                                           
c   (ftmaxg) - inlinable function to compute parcel temperature                 
c                                                                               
c common blocks:                                                                
c   comma    - scaling parameters and table for function ftma.                  
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      parameter(cp=cp_,rd=rd_)                                                  
      parameter(rocp=rd/cp)                                                     
      parameter(nx=151,ny=121)                                                  
      dimension tbtma(nx,ny),tbqma(nx,ny)                                       
      common/comma/ c1xma,c2xma,c1yma,c2yma,tbtma,tbqma                         
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      xmin=200.                                                                 
      xmax=500.                                                                 
      xinc=(xmax-xmin)/(nx-1)                                                   
      c1xma=1.-xmin/xinc                                                        
      c2xma=1./xinc                                                             
      ymin=0.01**rocp                                                           
      ymax=1.10**rocp                                                           
      yinc=(ymax-ymin)/(ny-1)                                                   
      c1yma=1.-ymin/yinc                                                        
      c2yma=1./yinc                                                             
      do jy=1,ny                                                                
        y=ymin+(jy-1)*yinc                                                      
        pk=y                                                                    
        t=xmin*y                                                                
        do jx=1,nx                                                              
          x=xmin+(jx-1)*xinc                                                    
          the=x                                                                 
          t=ftmaxg(t,the,pk,q)                                                  
          tbtma(jx,jy)=t                                                        
          tbqma(jx,jy)=q                                                        
        enddo                                                                   
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
