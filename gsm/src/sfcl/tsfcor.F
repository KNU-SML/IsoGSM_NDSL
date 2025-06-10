      subroutine tsfcor(tsfc,orog,slmask,umask,ijdim)                    
c
      real rlapse
      parameter(rlapse=0.65e-2)
c                                                                               
      dimension tsfc(ijdim),orog(ijdim),slmask(ijdim)                           
c                                                                               
      do ij=1,ijdim                                                             
        if(slmask(ij).eq.umask) then                                            
          tsfc(ij)=tsfc(ij)-orog(ij)*rlapse                                     
        endif                                                                   
      enddo                                                                     
      return                                                                    
      end                                                                       
