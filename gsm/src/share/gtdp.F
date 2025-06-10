      subroutine gtdp                                                           
#ifdef CRAY_THREAD
cfpp$ expand(ftdpxg)                                                            
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: gtdp         compute dewpoint temperature table                   
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: compute dewpoint temperature table as a function of                 
c   vapor pressure for inlinable function ftdp.                                 
c   exact dewpoint temperatures are calculated in subprogram ftdpxg.            
c   the current implementation computes a table with a length                   
c   of 5001 for vapor pressures ranging from 0.001 to 10.001 kilopascals        
c   giving a dewpoint temperature range of 208.0 to 319.0 kelvin.               
c                                                                               
c program history log:                                                          
c   91-05-07  iredell                                                           
c   94-12-30  iredell             expand table                                  
c                                                                               
c usage:  call gtdp                                                             
c                                                                               
c subprograms called:                                                           
c   (ftdpxg) - inlinable function to compute dewpoint temperature               
c                                                                               
c common blocks:                                                                
c   comtdp   - scaling parameters and table for function ftdp.                  
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
      parameter(nx=5001)                                                        
      dimension tbtdp(nx)                                                       
      common/comtdp/ c1xtdp,c2xtdp,tbtdp                                        
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      xmin= 0.001                                                               
      xmax=10.001                                                               
      xinc=(xmax-xmin)/(nx-1)                                                   
      c1xtdp=1.-xmin/xinc                                                       
      c2xtdp=1./xinc                                                            
      t=208.0                                                                   
      do jx=1,nx                                                                
        x=xmin+(jx-1)*xinc                                                      
        pv=x                                                                    
        t=ftdpxg(t,pv)                                                          
        tbtdp(jx)=t                                                             
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
