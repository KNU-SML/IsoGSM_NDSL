      function isrchflt(n,a,is,val)                                             
      dimension a(n)                                                            
      isrchflt=n                                                                
      do nn=is,n                                                                
      if( a(nn).lt.val ) then                                                   
        isrchflt=nn                                                             
        return                                                                  
      endif                                                                     
      enddo                                                                     
      return                                                                    
      end                                                                       
