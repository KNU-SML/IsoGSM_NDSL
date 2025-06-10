#define CLIQ cliq
#define CP cp
#define CVAP cvap
#define HVAP hvap
#define PSAT psat
#define RD rd
#define RV rv
#define TTP ttp
      function fthex(t,pk)                                                      
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: fthex        compute saturation vapor pressure                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: exactly compute equivalent potential temperature at the lcl         
c   from temperature and pressure over 100 kpa to the kappa power.              
c   equivalent potential temperature is constant for a saturated parcel         
c   rising adiabatically up a moist adiabat when the heat and mass              
c   of the condensed water are neglected. the formula for                       
c   equivalent potential temperature (derived in holton) is                     
c       the=t*(pd**(-rocp))*exp(el*eps*pv/(cp*t*pd))                            
c   where t is the temperature, pv is the saturated vapor pressure,             
c   pd is the dry pressure p-pv, el is the temperature dependent                
c   latent heat of condensation hvap+dldt*(t-ttp), and other values             
c   are physical constants defined in parameter statements in the code.         
c   zero is returned if the input values make saturation impossible.            
c   this function should be expanded inline in the calling routine.             
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             exact computation                             
c                                                                               
c usage:   the=fthex(t,pk)                                                      
c                                                                               
c   input argument list:                                                        
c     t        - real lcl temperature in kelvin                                 
c     pk       - real lcl pressure over 100 kpa to the kappa power              
c                                                                               
c   output argument list:                                                       
c     fthex    - real equivalent potential temperature in kelvin                
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
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      p=pk**cpor                                                                
      tr=ttp/t                                                                  
      pv=psatb*(tr**xa)*exp(xb*(1.-tr))                                         
      pd=p-pv                                                                   
      if(pd.gt.0.) then                                                         
        el=hvap+dldt*(t-ttp)                                                    
c       expo=el*eps*pv/(cp*t*pd)                                                
c       fthex=t*pd**(-rocp)*exp(expo)                                           
        expo=el*eps*pv/(cp*t*pd)                                                
        expo = min(expo,100.0)                                                  
        fthex=t*pd**(-rocp)*exp(expo)                                           
      else                                                                      
        fthex=0.                                                                
      endif                                                                     
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -         
      return                                                                    
      end                                                                       
