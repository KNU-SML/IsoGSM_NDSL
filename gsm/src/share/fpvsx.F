#define CLIQ cliq
#define CICE cice
#define CP cp
#define CVAP cvap
#define HVAP hvap
#define HSUB hsub
#define PSAT psat
#define RD rd
#define RV rv
#define TTP ttp
      function fpvsx(t)                                                         
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c$$$     subprogram documentation block                                         
c                                                                               
c subprogram: fpvsx        compute saturation vapor pressure                    
c   author: n phillips            w/nmc2x2   date: 30 dec 82                    
c                                                                               
c abstract: exactly compute saturation vapor pressure from temperature.         
c   the water model assumes a perfect gas, constant specific heats              
c   for gas and liquid, and neglects the volume of the liquid.                  
c   the model does account for the variation of the latent heat                 
c   of condensation with temperature.  the ice option is not included.          
c   the clausius-clapeyron equation is integrated from the triple point         
c   to get the formula                                                          
c       pvs=psatk*(tr**xa)*exp(xb*(1.-tr))                                      
c   where tr is ttp/t and other values are physical constants                   
c   this function should be expanded inline in the calling routine.             
c                                                                               
c program history log:                                                          
c   91-05-07  iredell             made into inlinable function                  
c   94-12-30  iredell             exact computation                             
c                                                                               
c usage:   pvs=fpvsx(t)                                                         
c                                                                               
c   input argument list:                                                        
c     t        - real temperature in kelvin                                     
c                                                                               
c   output argument list:                                                       
c     fpvsx    - real saturation vapor pressure in kilopascals (cb)             
c                                                                               
c attributes:                                                                   
c   language: fortran 77.                                                       
c   machine:  cray.                                                             
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      parameter(cp=cp_,rd=rd_,rv=rv_,                                           
     &          ttp=ttp_,hvap=hvap_,psat=psat_,                                 
     &          cliq=cliq_,cvap=cvap_,cice=cice_,hsub=hsub_)
      parameter(psatk=psat*1.e-3)
      parameter(dldt=cvap-cliq,xa=-dldt/rv,xb=xa+hvap/(rv*ttp))
      parameter(dldti=cvap-cice,xai=-dldti/rv,xbi=xai+hsub/(rv*ttp))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      tr=ttp/t
      if(t.ge.ttp) then
      fpvsx=psatk*(tr**xa)*exp(xb*(1.-tr))
      else
      fpvsx=psatk*(tr**xai)*exp(xbi*(1.-tr))
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      return
      end
