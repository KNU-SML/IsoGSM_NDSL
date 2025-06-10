      subroutine wryte(lu,lc,c)                                                 
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    wryte       write data out by bytes                            
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: efficiently write unformatted a characeter array.                   
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:    call wryte(lu,lc,c)                                                 
c                                                                               
c   input argument list:                                                        
c     lu       - integer unit to which to write                                 
c     lc       - integer number of characters or bytes to write                 
c     c        - characeter (lc) data to write                                  
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      character c(lc)                                                           
c
      write(lu) c                                                               
c
      return                                                                    
      end                                                                       
