      subroutine chgres(maxwv,wave,maxwvc,wavec)                                
      save                                                                      
c                                                                               
      dimension wave(*),wavec(*)                                                
c                                                                               
      mdim=(maxwv+1)*(maxwv+2)                                                  
c                                                                               
      nmc=0                                                                     
      do 10 m=1,maxwvc+1                                                        
      nend=maxwvc+1                                                             
      do 10 n=m,nend                                                            
      nmc=nmc+1                                                                 
      nm=m*(maxwv+1)-(maxwv+1-n)-(m-1)*m/2                                      
      do 10 ii=1,2                                                              
        if(m.gt.maxwv+1.or.n.gt.maxwv+1) then
          wavec(nmc*2-2+ii)=0.
        else
          if(nm*2-2+ii.le.0) then
            print *,'logic error'
            call abort
          endif
          wavec(nmc*2-2+ii)=wave(nm*2-2+ii)
        endif
   10 continue                                                                  
c                                                                               
      return                                                                    
      end                                                                       
