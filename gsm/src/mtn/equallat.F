      subroutine equallat(k,a,w)                                                
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:  equallat   compute equally-spaced latitudes                      
c   prgmmr: iredell       org: w/nmc23       date: 92-04-16                     
c                                                                               
c abstract: computes cosines of colatitude and gaussian weights                 
c   on equally-spaced latitudes.                                                
c                                                                               
c program history log:                                                          
c   92-04-16  iredell                                                           
c                                                                               
c usage:    call equallat(k,a,w)                                                
c                                                                               
c   input argument list:                                                        
c     k        - number of latitudes.                                           
c                                                                               
c   output argument list:                                                       
c     a        - real (k) cosines of colatitude.                                
c     w        - real (k) gaussian weights.                                     
c                                                                               
c   remarks: fortran 9x extensions are used.                                    
c                                                                               
c attributes:                                                                   
c   cray ymp.                                                                   
c                                                                               
c$$$                                                                            
      save                                                                      
      dimension a(k),w(k)                                                       
      pi=4.*atan(1.)                                                            
      kh=k/2                                                                    
      dlt=pi/(k-1)                                                              
#ifdef CRAY_THREAD
cdir$ ivdep                                                                     
#endif
      do j=1,kh                                                                 
        a(j)=cos((j-1)*dlt)                                                     
        a(k+1-j)=-a(j)                                                          
      enddo                                                                     
      w(1)=1.-cos(dlt*0.5)                                                      
      w(k)=w(1)                                                                 
      sindlt=2.*sin(dlt*0.5)                                                    
#ifdef CRAY_THREAD
cdir$ ivdep                                                                     
#endif
      do j=2,kh                                                                 
        w(j)=sin((j-1)*dlt)*sindlt                                              
        w(k+1-j)=w(j)                                                           
      enddo                                                                     
      if(k.ne.kh*2) then                                                        
        a(kh+1)=0.                                                              
        w(kh+1)=sindlt*0.5                                                      
      endif                                                                     
      return                                                                    
      end                                                                       
