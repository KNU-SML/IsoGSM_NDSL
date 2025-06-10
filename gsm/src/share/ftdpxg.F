#define CLIQ cliq
#define CP cp
#define CVAP cvap
#define HVAP hvap
#define PSAT psat
#define RD rd
#define RV rv
#define TTP ttp
      function ftdpxg(tg,pv)                                                    
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: ftdpxg       compute saturation vapor pressure                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: exactly compute dewpoint temperature from vapor pressure.           
c   a guess dewpoint temperature must be provided.                              
c   the water model assumes a perfect gas, constant specific heats              
c   for gas and liquid, and neglects the volume of the liquid.                  
c   the model does account for the variation of the latent heat                 
c   of condensation with temperature.  the ice option is not included.          
c   the clausius-clapeyron equation is integrated from the triple point         
c   to get the formula                                                          
c       pvs=psatk*(tr**xa)*exp(xb*(1.-tr))                                      
c   where tr is ttp/t and other values are physical constants                   
c   the formula is inverted by iterating newtonian approximations               
c   for each pvs until t is found to within 1.e-6 kelvin.                       
c   this function can be expanded inline in the calling routine.                
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             exact computation                             
c                                                                               
c usage:   tdp=ftdpxg(tg,pv)                                                    
c                                                                               
c   input argument list:                                                        
c     tg       - real guess dewpoint temperature in kelvin                      
c     pv       - real vapor pressure in kilopascals (cb)                        
c                                                                               
c   output argument list:                                                       
c     ftdpxg   - real dewpoint temperature in kelvin                            
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      parameter(cp=cp_,rd=rd_,rv=rv_,                                           
     &          ttp=ttp_,hvap=hvap_,psat=psat_,                                 
     &          cliq=cliq_,cvap=cvap_)                                          
      parameter(psatk=psat*1.e-3)                                               
      parameter(dldt=cvap-cliq,xa=-dldt/rv,xb=xa+hvap/(rv*ttp))                 
      parameter(terrm=1.e-6)                                                    
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      t=tg                                                                      
      tr=ttp/t                                                                  
      pvt=psatk*(tr**xa)*exp(xb*(1.-tr))                                        
      el=hvap+dldt*(t-ttp)                                                      
      dpvt=el*pvt/(rv*t**2)                                                     
      terr=(pvt-pv)/dpvt                                                        
      t=t-terr                                                                  
      dowhile(abs(terr).gt.terrm)                                               
        tr=ttp/t                                                                
        pvt=psatk*(tr**xa)*exp(xb*(1.-tr))                                      
        el=hvap+dldt*(t-ttp)                                                    
        dpvt=el*pvt/(rv*t**2)                                                   
        terr=(pvt-pv)/dpvt                                                      
        t=t-terr                                                                
      enddo                                                                     
      ftdpxg=t                                                                  
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
