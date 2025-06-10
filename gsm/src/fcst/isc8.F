      function isc8(n,ca,ci)                                                    
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    isc8        find a character*8 string.                         
c   prgmmr: mark iredell     org: w/nmc23    date: 91-03-15                     
c                                                                               
c abstract: integer function to return the index of a character*8 name.         
c                                                                               
c program history log:                                                          
c   91-03-15  mark iredell                                                      
c                                                                               
c usage:    i=isc8 (n,ca,ci)                                                    
c   input argument list:                                                        
c     n        - number of strings to search                                    
c     ca       - character*8 array to search                                    
c     ci       - character*8 string for which to search                         
c                                                                               
c   output argument list:                                                       
c     isc8     - index of first string found or 0 if not found                  
c                                                                               
c remarks: this routine may be inlined.                                         
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray ymp.                                                         
c                                                                               
c$$$                                                                            
      character*8 ca(n),ci                                                      
      do 10 i=1,n                                                               
      if(ca(i).eq.ci) then                                                      
        isc8=i                                                                  
        return                                                                  
      endif                                                                     
10    continue                                                                  
      isc8=0                                                                    
      return                                                                    
      end                                                                       
