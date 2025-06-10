      subroutine avlat(xave,x,idim,jdim,kdim,krec,jpts)                         
      save                                                                      
      dimension xave(kdim,jdim),x(idim,jdim),jpts(jdim)                         
      do 1 k=1,kdim                                                             
        do 2 j=1,jdim                                                           
          xave(j,k) = 0.                                                        
    2   continue                                                                
        do 4 j=1,jdim                                                           
         numpts = 0                                                             
         do 3 i=1,idim                                                          
           if (x(i,j).ge.0.) then                                               
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
