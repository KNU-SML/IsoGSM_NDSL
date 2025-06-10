#define LATG2 latg2
#define LATG2P latg2p
#define LATG latg
#define LATGP latgp
#define LEVS levs
#define LONF2 lonf2
#define LONF2P lonf2p
#define LONF lonf
#define LONFP lonfp

      subroutine znlavz(latb2,lonb2,secphy,secrad,                              
     &                  dumtn,slmsk,sheleg,dswsfc,uswsfc)                       
#ifdef CRAY_THREAD
cfpp$ expand(znlais)                                                            
#endif
c                                                                               
#include <paramodel.h>
      dimension dumtn(lonb2,latb2)                                              
      dimension slmsk(lonb2,latb2),sheleg(lonb2,latb2)                          
      dimension dswsfc(lonb2,latb2),uswsfc(lonb2,latb2)                         
c                                                                               
#include <comznl.h>
c                                                                               
#ifdef DYNAMIC_ALLOC
      dimension ist(lonb2)                                                      
      dimension work1(lonb2)                                                    
#else
#ifdef MP
#define LONF2S lonf2p_
#else
#define LONF2S lonf2_
#endif
      dimension ist(LONF2S)                                                     
      dimension work1(LONF2S)                                                   
#endif
c                                                                               
#ifdef CRAY_THREAD
cfpp$ privatearray                                                              
#endif
      do lat=1,latb2                                                            
c                                                                               
        do i=1,lonb2                                                            
          if(slmsk(i,lat).eq.0.) then                                           
            ist(i)=6                                                            
          elseif(slmsk(i,lat).eq.1.) then                                       
            ist(i)=2                                                            
            if(sheleg(i,lat).gt.1.e-3) ist(i)=3                                 
          else                                                                  
            ist(i)=4                                                            
            if(sheleg(i,lat).gt.1.e-3) ist(i)=5                                 
          endif                                                                 
        enddo                                                                   
c                                                                               
        if(secrad.gt.0.) then                                                   
          do i=1,lonb2                                                          
            work1(i)=uswsfc(i,lat)-dswsfc(i,lat)                                
          enddo                                                                 
          call znlais(work1,ist,1./secrad,nssfcsw,lat)                          
        endif                                                                   
c                                                                               
        if(secphy.gt.0.) then                                                   
          call znlais(dumtn(1,lat),ist,1./secphy,nsuamtn,lat)                   
        endif                                                                   
c                                                                               
      enddo                                                                     
c                                                                               
      return                                                                    
      end                                                                       
