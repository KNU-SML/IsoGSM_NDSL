      function isrcheq(n,a,is,val)                                              
      dimension a(n)                                                            
      isrcheq=n                                                                 
      do nn=is,n                                                                
      if( a(nn).eq.val ) then                                                   
        isrcheq=nn                                                              
        return                                                                  
      endif                                                                     
      enddo                                                                     
      return                                                                    
      end                                                                       
