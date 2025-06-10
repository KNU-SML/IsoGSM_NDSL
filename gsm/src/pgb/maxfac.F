      function maxfac(n)                                                        
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    maxfac      return maximum prime factor                        
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: determines the maximum prime factor of a positive integer.          
c           useful for determining fitness for fft factorization.               
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:   ...=maxfac(n)                                                        
c                                                                               
c   input argument list:                                                        
c     n        - integer number to factor                                       
c                                                                               
c   output argument list:                                                       
c     maxfac   - maximum prime factor of n                                      
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      in=n                                                                      
      k=2                                                                       
      m=1                                                                       
      dowhile(in.gt.1)                                                          
        inx=in/k                                                                
        if(in.eq.inx*k) then                                                    
          in=inx                                                                
          m=k                                                                   
        else                                                                    
          k=k+1                                                                 
          if(k.gt.3) k=k+1                                                      
          if(k.gt.inx) k=in                                                     
        endif                                                                   
      enddo                                                                     
      maxfac=m                                                                  
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
