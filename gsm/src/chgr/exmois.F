#define IJCAP ijcap
#define ILATG ilatg
#define ILEVS ilevs
#define ILONF ilonf
#define JCAP jcap
#define LATG latg
#define LDIAG ldiag
#define LEVS levs
#define LGRID lgrid
#define LONF lonf
#define LPNT lpnt
#define LSOIL lsoil
#define LTSTP ltstp
#define LVAR lvar
#define PSAT psat
#define RD rd
#define RV rv
#define T0C t0c
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine exmois(pstar,qin,t,chgrp)                                      
#include <paramodel.h>
      save                                                                      
c                                                                               
#include <paramter.h>
c                                                                               
#include <parmchgr.h>
c                                                                               
      parameter(kdimpi=kdimi+1,kdimmi=kdimi-1)                                  
c                                                                               
      dimension qin(idim,kdimi),t(idim,kdimi),qout(kdimi)                      
      dimension qsat(kdimi),prs(kdimi),pstar(idim)                              
      dimension chgrp(idim,kdimi)                                               
c hmhj 2                                                                        
      common /inpver/ siin(kdimpi),slin(kdimi),delin(kdimi),                    
     1                ciin(kdimpi),clin(kdimi),rpiin(kdimmi)                    
     2               ,ak(kdimpi),bk(kdimpi)                                     ! 20070501(12)
c---------------------------------------------------------------                
c----    input - sfc pres (cb) - pstar                                          
c              - moisture      - qin - levh lyrs only                           
c              - temperature   - t   - model virtual temperature                
c----    output - levs moisture lyrs - chgrp                                    
c               - thermodynamic temperature back into t                         
c               - virtual temperature should be restored (hmhj)                 
c----    this code extrapolates moisture up into model dry layers               
c          using exponential decrease up to a minimum valu of q                 
c----    first convert mdl virtual temp to thermodynmic temp                    
c---------------------------------------------------------------                
      data itst / 0 /                                                           
c                                                                               
c                                                                               
      if (itst.gt.0) go to 1                                                    
      pmin = 5.                                                                 
      qmin = 3.e-6                                                              
      xlnpm = log(pmin)                                                         
      xlnqm = log(qmin)                                                         
      itst = 1                                                                  
    1 continue                                                                  
c
      rmaxsiin=siin(1)
      do k=1,kdimpi
        rmaxsiin=max(rmaxsiin,siin(k))
      enddo
c
      do 50 i=1,idim                                                            
c----    model temp is virtual (since hydrostatic from hgts initially),         
c          so we convert to thermodynamic temp for radiation calculation        
      do 2 k=1,kdimi                                                           
        if (qin(i,k).le.0.0) go to 2                                            
        t(i,k) = t(i,k) * (1.0 + qin(i,k))                                      
     1                   /(1.0 + rv_/rd_ * qin(i,k))                            
   2  continue                                                                  
c                                                                               
c----    get lyr pressure and then saturated moisture                           
      do 5 k=1,kdimi                                                            
      if(rmaxsiin.gt.1. or. rmaxsiin.eq.0. )then
        prs(k) = pstar(ix) * bk(k) + ak(k)
      else
        prs(k) = pstar(ix) * slin(k)                                              
      endif
    5 continue                                                                  
c----    compute saturation specific humidity(dimensionless) from               
c          temperature t (deg k) and pressure (cb)                              
c        method is tetens formula for saturation vapor pressure (mb) ,          
c          given by haurwitz p. 9 - dynamic meteorology , 1941.                 
c----    conversion to specific humidity follows from definition                
c          ibid. , p.10                                                         
c                                                                               
      do 7 k=1,kdimi                                                            
      expon = 7.50*(t(i,k)-t0c_)/((t(i,k)-t0c_)+237.30)                         
      es = psat_*1.e-2 * 10.0**expon                                            
      qsat(k)=rd_/rv_*es/(prs(k)*10.0-(1.0-rd_/rv_)*es)                         
    7 continue                                                                  
c----    limit moisture in lowest levh layers--rh le 1. but ge .15              
c          the latter avoids the negative q problem...                          
      do 10 k=1,kdimi                                                          
        qout(k) = qin(i,k)                                                      
        rh = qin(i,k) / qsat(k)                                                 
   10 continue                                                                  
      ldry = kdimi + 1                                                         
c hmhj 1                                                                        
      if( ldry .lt. kdimi ) then                                                
c                                                                               
c----    extrapolate moisture to min valu(qmin) at pressure pmin(cb)            
c         use exponential decrease from layer levh----                          
c         i.e.  q=q(levh)*(p/p(levh)) ** xlamb                                  
c        where xlamb is computed to fit qmin at pmin and q,p at levh.           
      xlnpq = log(prs(kdimi))                                                  
      xlnqq = log(qout(kdimi))                                                 
      xlamb = (xlnqm-xlnqq) / (xlnpm-xlnpq)                                     
      do 15 k=ldry,kdimi                                                        
        qout(k) = qout(kdimi)*(prs(k)/prs(kdimi)) ** xlamb                    
        if (prs(k).le.pmin) qout(k) = qmin                                      
        if (qout(k).le.qmin) qout(k) = qmin                                     
        rh = qout(k) / qsat(k)                                                  
        if (rh.gt.1.0) qout(k) = qsat(k)                                        
   15 continue                                                                  
c hmhj 1                                                                        
      endif                                                                     
c----    store extrapolated moisture                                            
      do 20 k=1,kdimi                                                           
        chgrp(i,k) = qout(k)                                                    
   20 continue                                                                  
c ---- temp. is restored to virtual temp.                                       
      do 25 k=1,kdimi                                                          
        if (qin(i,k).le.0.0) go to 25                                           
        t(i,k) = t(i,k) * (1.0 + rv_/rd_ * qin(i,k))                            
     1                   /(1.0 + qin(i,k))                                      
   25 continue                                                                  
   50 continue                                                                  
      return                                                                    
      end                                                                       
