      subroutine setes                                                          
      parameter(nx=7501)                                                        
      dimension tbsvp(nx), tbdsvp(nx)                                           
      common/comsvp/ c1xsvp,c2xsvp,tbsvp, tbdsvp                                
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
      xmin=180.0                                                                
      xmax=330.0                                                                
      xinc=(xmax-xmin)/(nx-1)                                                   
      c1xsvp=1.-xmin/xinc                                                       
      c2xsvp=1./xinc                                                            
#ifdef CRAY_THREAD
cfpp$ expand (esatx)                                                            
#endif
      do jx=1,nx                                                                
        x=xmin+(jx-1)*xinc                                                      
        t=x                                                                     
        call esatx(t, tbsvp(jx), tbdsvp(jx))                                    
      enddo                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
