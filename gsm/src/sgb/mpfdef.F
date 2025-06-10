      subroutine mpfdef(iptv,mpf)                                               
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram: mpfdef         sets default pole vector flags                     
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: sets field identifier defaults for various parameters.              
c   a flag of 0 means scalar, 1 means vector, and 2 means flag.                 
c   these identifiers are used in interpolation.                                
c                                                                               
c program history log:                                                          
c   93-10-21  iredell                                                           
c                                                                               
c usage:    call mpfdef(iptv,mpf)                                               
c   input arguments:                                                            
c     iptv         paramter table version (only 1 or 2 is recognized)           
c   output arguments:                                                           
c     mpf          integer (255) field parameter identifiers                    
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      dimension mpf(255)                                                        
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      do i=1,255                                                                
      mpf(i)=0                                                                  
      enddo                                                                     
      if(iptv.eq.1.or.iptv.eq.2) then                                           
        do i=33,34                                                              
c       mpf(033:034)=1                                                          
        mpf(i)=1                                                                
        enddo                                                                   
        do i=49,50                                                              
c       mpf(049:050)=1                                                          
        mpf(i)=1                                                                
        enddo                                                                   
        do i=95,96                                                              
c       mpf(095:096)=1                                                          
        mpf(i)=1                                                                
        enddo                                                                   
        do i=124,125                                                            
c       mpf(124:125)=1                                                          
        mpf(i)=1                                                                
        enddo                                                                   
        do i=181,182                                                            
c       mpf(181:182)=1                                                          
        mpf(i)=1                                                                
        enddo                                                                   
        do i=183,184                                                            
c       mpf(183:184)=1                                                          
        mpf(i)=1                                                                
        enddo                                                                   
        do i=247,248                                                            
c       mpf(247:248)=1                                                          
        mpf(i)=1                                                                
        enddo                                                                   
        mpf(081)=2                                                              
        mpf(091)=2                                                              
        mpf(140)=2                                                              
        mpf(141)=2                                                              
        mpf(142)=2                                                              
        mpf(143)=2                                                              
        mpf(173)=2                                                              
        mpf(174)=2                                                              
        mpf(175)=2                                                              
        mpf(209)=2                                                              
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
