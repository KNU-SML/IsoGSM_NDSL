#include <define.h>
      function isrchfltx(n,a,is,val)                                            
      dimension a(n)                                                            
      isrchfltx=n                                                               
      do nn=is,n                                                                
      if( a(nn).lt.val ) then                                                   
        isrchfltx=nn                                                            
        return                                                                  
      endif                                                                     
      enddo                                                                     
      return                                                                    
      end                                                                       
