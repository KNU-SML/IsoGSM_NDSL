      subroutine avelat(xave,x,idim,jdim,kdim,krec,jpts)                        
      save                                                                      
      dimension xave(kdim,jdim),x(idim,jdim),jpts(jdim)                         
      ktp = 0                                                                   
      if (krec.ge.10.and.krec.le.11) ktp = 1                                    
      if (krec.ge.13.and.krec.le.14) ktp = 1                                    
      if (krec.ge.16.and.krec.le.17) ktp = 1                                    
      do 1 k=1,kdim                                                             
        do 2 j=1,jdim                                                           
          xave(j,k) = 0.                                                        
    2   continue                                                                
        do 4 j=1,jdim                                                           
         numpts = 0                                                             
         do 3 i=1,idim                                                          
          if (ktp.eq.1) then                                                    
c....    cloud top or bottom                                                    
           if(x(i,j).gt.0.) then                                                
            xave(j,k) = xave(j,k) + x(i,j)                                      
            numpts = numpts + 1                                                 
           end if                                                               
          else                                                                  
           xave(j,k) = xave(j,k) + x(i,j)                                       
           numpts = numpts + 1                                                  
          end if                                                                
    3    continue                                                               
         jpts(j) = numpts                                                       
c...   averaging done in subroutine globlm.....                                 
    4   continue                                                                
    1 continue                                                                  
      return                                                                    
      end                                                                       
