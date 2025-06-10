      subroutine c2z(l,c,z)                                                     
c$$$  subprogram documentation block                                            
c                                                                               
c subprogram:    c2z         convert byte to hexadecimal character pair         
c   prgmmr: iredell          org: w/nmc23     date: 92-10-31                    
c                                                                               
c abstract: converts an array of bytes to its hexadecimal representation        
c   (2 characters per byte) for diagnostic purposes.                            
c                                                                               
c program history log:                                                          
c   91-10-31  mark iredell                                                      
c                                                                               
c usage:    call c2z(l,c,z)                                                     
c                                                                               
c   input argument list:                                                        
c     l        - integer number of bytes to represent                           
c     c        - character (l) byte data to convert                             
c                                                                               
c   output argument list:                                                       
c     z        - character (2*l) hexadecimal representation                     
c                                                                               
c attributes:                                                                   
c   language: cray fortran                                                      
c                                                                               
c$$$                                                                            
      character c(l)*1,z(l)*2                                                   
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      do i=1,l                                                                  
        write(z(i),'(z2)') ichar(c(i))                                          
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
