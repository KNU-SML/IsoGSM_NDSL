      subroutine gpvs                                                           
#ifdef CRAY_THREAD
cfpp$ expand(fpvsx,fpvsx0)
#endif
#include <paramodel.h>
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: gpvs         compute saturation vapor pressure table              
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute saturation vapor pressure table as a function of            
c   temperature for the table lookup function fpvs.                             
c   exact saturation vapor pressures are calculated in subprogram fpvsx.        
c   the current implementation computes a table with a length                   
c   of 7501 for temperatures ranging from 180. to 330. kelvin.                  
c                                                                               
c program history log:                                                          
c   91-05-07  iredell                                                           
c   94-12-30  iredell             expand table                                  
C   96-02-19  hong                ice effect
c                                                                               
c usage:  call gpvs                                                             
c                                                                               
c subprograms called:                                                           
c   (fpvsx)  - inlinable function to compute saturation vapor pressure          
c                                                                               
c common blocks:                                                                
c   compvs   - scaling parameters and table for function fpvs.                  
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
      parameter(nx=7501)                                                        
      dimension tbpvs(nx),tbpvs0(nx)
      common/compvs0/ c1xpvs0,c2xpvs0,tbpvs0
      common/compvs/ c1xpvs,c2xpvs,tbpvs
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      xmin=180.0
      xmax=330.0
      xinc=(xmax-xmin)/(nx-1)
      c1xpvs=1.-xmin/xinc
      c2xpvs=1./xinc
      c1xpvs0=1.-xmin/xinc
      c2xpvs0=1./xinc
      do jx=1,nx
        x=xmin+(jx-1)*xinc
        t=x
        tbpvs(jx)=fpvsx(t)
        tbpvs0(jx)=fpvsx0(t)
      enddo
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
