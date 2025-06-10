#define CLIQ cliq
#define CP cp
#define CVAP cvap
#define HVAP hvap
#define PSAT psat
#define RD rd
#define RV rv
#define TTP ttp
      function ftmaxg(tg,the,pk,qma)                                            
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: ftmaxg       compute moist adiabat temperature                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: exactly compute temperature and humidity of a parcel                
c   lifted up a moist adiabat from equivalent potential temperature             
c   at the lcl and pressure over 100 kpa to the kappa power.                    
c   a guess parcel temperature must be provided.                                
c   equivalent potential temperature is constant for a saturated parcel         
c   rising adiabatically up a moist adiabat when the heat and mass              
c   of the condensed water are neglected. the formula for                       
c   equivalent potential temperature (derived in holton) is                     
c       the=t*(pd**(-rocp))*exp(el*eps*pv/(cp*t*pd))                            
c   where t is the temperature, pv is the saturated vapor pressure,             
c   pd is the dry pressure p-pv, el is the temperature dependent                
c   latent heat of condensation hvap+dldt*(t-ttp), and other values             
c   are physical constants defined in parameter statements in the code.         
c   the formula is inverted by iterating newtonian approximations               
c   for each the and p until t is found to within 1.e-4 kelvin.                 
c   the specific humidity is then computed from pv and pd.                      
c   this function can be expanded inline in the calling routine.                
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             exact computation                             
c                                                                               
c usage:   tma=ftmaxg(tg,the,pk,qma)                                            
c                                                                               
c   input argument list:                                                        
c     tg       - real guess parcel temperature in kelvin                        
c     the      - real equivalent potential temperature in kelvin                
c     pk       - real pressure over 100 kpa to the kappa power                  
c                                                                               
c   output argument list:                                                       
c     ftmaxg   - real parcel temperature in kelvin                              
c     qma      - real parcel specific humidity in kg/kg                         
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
      parameter(rocp=rd/cp,cpor=cp/rd,psatb=psatk*1.e-2,eps=rd/rv,              
     &          dldt=cvap-cliq,xa=-dldt/rv,xb=xa+hvap/(rv*ttp))                 
      parameter(terrm=1.e-4)                                                    
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      t=tg                                                                      
      p=pk**cpor                                                                
      tr=ttp/t                                                                  
      pv=psatb*(tr**xa)*exp(xb*(1.-tr))                                         
      pd=p-pv                                                                   
      el=hvap+dldt*(t-ttp)                                                      
      expo=el*eps*pv/(cp*t*pd)                                                  
      thet=t*pd**(-rocp)*exp(expo)                                              
      dthet=thet/t*(1.+expo*(dldt*t/el+el*p/(rv*t*pd)))                         
      terr=(thet-the)/dthet                                                     
      t=t-terr                                                                  
      dowhile(abs(terr).gt.terrm)                                               
        tr=ttp/t                                                                
        pv=psatb*(tr**xa)*exp(xb*(1.-tr))                                       
        pd=p-pv                                                                 
        el=hvap+dldt*(t-ttp)                                                    
        expo=el*eps*pv/(cp*t*pd)                                                
        thet=t*pd**(-rocp)*exp(expo)                                            
        dthet=thet/t*(1.+expo*(dldt*t/el+el*p/(rv*t*pd)))                       
        terr=(thet-the)/dthet                                                   
        t=t-terr                                                                
      enddo                                                                     
      ftmaxg=t                                                                  
      tr=ttp/t                                                                  
      pv=psatb*(tr**xa)*exp(xb*(1.-tr))                                         
      pd=p-pv                                                                   
      qma=eps*pv/(pd+eps*pv)                                                    
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
