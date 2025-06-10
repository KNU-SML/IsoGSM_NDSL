      function path(a,b,c,e)                                                    
c....                                                                           
c     double precision xa,ca                                                    
      common/coefs/xa(109),ca(109),eta(109),sexpv(109),core,uexp,sexp           
      pexp=1./sexp                                                              
      path=((a-b)**pexp*(a+b+c))/(e*(a+b+c)+(a-b)**(pexp-1.))                   
      return                                                                    
      end                                                                       
