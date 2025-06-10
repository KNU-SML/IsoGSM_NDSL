      subroutine qcsli(slianl,slifcs,ijdim)                                     
      dimension slianl(ijdim),slifcs(ijdim)                                     
c     write(6,*) ' '                                                            
c     write(6,*) 'qcsli'                                                        
      kount=0                                                                   
      do ij=1,ijdim                                                             
        if(slianl(ij).eq.1..and.slifcs(ij).eq.0.) then                          
          kount=kount+1                                                         
          slifcs(ij)=1.                                                         
        endif                                                                   
        if(slianl(ij).eq.0..and.slifcs(ij).eq.1.) then                          
          kount=kount+1                                                         
          slifcs(ij)=0.                                                         
        endif                                                                   
        if(slianl(ij).eq.2..and.slifcs(ij).eq.1.) then                          
          kount=kount+1                                                         
          slifcs(ij)=0.                                                         
        endif                                                                   
        if(slianl(ij).eq.1..and.slifcs(ij).eq.2.) then                          
          kount=kount+1                                                         
          slifcs(ij)=1.                                                         
        endif                                                                   
      enddo                                                                     
c     if(kount.gt.0) then                                                       
c       per=float(kount)/float(ijdim)*100.                                      
c       write(6,*) ' inconsistency of slmask between forecast and',             
c    1             ' analysis corrected at ',kount, ' points ',per,             
c    2             'percent'                                                    
c     endif                                                                     
      return                                                                    
      end                                                                       
