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
#define PI pi
#define PSAT psat
#define RD rd
#define RV rv
#define T0C t0c
#define NTOTAL ntotal
#define NCLDG ncldg
#define NTRAC ntrac
      subroutine exmoisn(pstar,qin,t,chgrp,lat,xlat)                            
#include <paramodel.h>
      save                                                                      
c                                                                               
c  this code extrapolates moisture up into model dry layers                     
c  using exponential decrease up to a specified q at 80 mb,                     
c  then use lims profile above(see middle atmos program handbook,               
c  vol 22,sep 86,ed. j.m.russell ---> j alpert has copy)....                    
c                                                                               
c  input - sfc pres (cb) - pstar                                                
c        - moisture      - qin - levh lyrs only                                 
c        - temperature   - t   - model virtual temperature                      
c        - xlat=lat in radians--                                                
c  output - levs moisture lyrs - chgrp (use lims data abuv 8 cb)                
c         - thermodynamic temperature back into t {revised to                   
c         - restore temp. to virtual temp. ( hmhj )}                            
c                                                                               
#include <paramter.h>
c                                                                               
#include <parmchgr.h>
c                                                                               
      parameter(kdimpi=kdimi+1,kdimmi=kdimi-1)                                  
c                                                                               
      parameter(idimt=idim*2)                                                   
c                                                                               
      real xlamb,expon                                                          
      dimension qin(idimt,kdimi),t(idimt,kdimi),qout(kdimi)                    
      dimension qsat(kdimi),prs(kdimi),pstar(idimt)                             
      dimension chgrp(idimt,kdimi)                                              
c hmhj 2                                                                        
      common /inpver/ siin(kdimpi),slin(kdimi),delin(kdimi),                    
     1                ciin(kdimpi),clin(kdimi),rpiin(kdimmi)                    
     2               ,ak(kdimpi),bk(kdimpi)                                     ! 20070501(12)
c                                                                               
c  qq stores lims (satellite-1986-see j.alpert) observed h2o above              
c  80 mb...(k,j)-k=1,8 for 80,70,..,10 mb                                       
c  j=1,19 for every 10 deg lat (starting at n.pole)                             
c                                                                               
      dimension  qq(8,19)                                                       
      data qq /32*5.  ,  4.8,4.55,4.5,4.6,4.7,4.75,4.75,4.7,                    
     1 4.,3.85,3.9,4.,4.25,4.5,4.5,4.6, 3.,2.75,3.,3.3,3.7,4.2,4.3,4.5,         
     2           2.1,2.,2.,2.5,3.,3.9,4.,4. ,                                   
     3 2.,2.,2.,2.3,2.8,3.5,3.7,3.75 , 2.,2.,2.,2.3,2.7,3.25,3.45,3.75,         
     4 2.,2.,2.,2.3,2.8,3.5,3.7,3.75 , 2.1,2.,2.,2.5,3.,3.9,4.,4. ,             
     5           3.,2.75,3.,3.3,3.7,4.2,4.3,4.5,                                
     6      4.,3.85,3.9,4.,4.25,4.5,4.5,4.6 ,                                   
     7           4.8,4.55,4.5,4.6,4.7,4.75,4.75,4.7 , 32*5.  /                  
c                                                                               
      data itst / 0 /                                                           
c                                                                               
c                                                                               
c  first convert mdl virtual temp to thermodynmic temp                          
c                                                                               
      if (itst.gt.0) go to 1                                                    
c                                                                               
      pmin = 8.0                                                                
      pmm  = 1.0                                                                
      xlnpm =  log(pmin)                                                        
c                                                                               
c  get correct units of lims moisture                                           
c                                                                               
      do 41 j=1,19                                                              
      do 41 k=1,8                                                               
      qq(k,j) = qq(k,j) * 1.e-6                                                 
   41 continue                                                                  
      itst = 1                                                                  
    1 continue                                                                  
c                                                                               
      rmaxsiin=siin(1)
      do k=1,kdimpi
        rmaxsiin=max(rmaxsiin,siin(k))
      enddo
c
c  loop for hemisphere                                                          
c                                                                               
      do 1000 ihem=1,2                                                          
c                                                                               
c  get latitude in degrees..                                                    
c                                                                               
      dlat = xlat * 180.0/pi_                                                   
      if (ihem.eq.2) dlat =-dlat                                                
c                                                                               
c  lon loop                                                                     
c                                                                               
      do 50 i=1,idim                                                            
c                                                                               
      ix=(ihem-1)*idim+i                                                        
c                                                                               
c  model temp is virtual (since hydrostatic from hgts initially),               
c                                                                               
      do 2 k=1,kdimi                                                           
      if (qin(ix,k).le.0.0) go to 2                                             
      t(ix,k)=t(ix,k)*(1.0+qin(ix,k))/(1.0+rv_/rd_*qin(ix,k))                   
   2  continue                                                                  
c                                                                               
c  get lyr pressure and then saturated moisture                                 
c                                                                               
      do 5 k=1,kdimi                                                            
c
      if(rmaxsiin.gt.1. or. rmaxsiin.eq.0. )then
        prs(k) = pstar(ix) * bk(k) + ak(k)
      else
        prs(k) = pstar(ix) * slin(k)                                              
      endif

    5 continue                                                                  
c                                                                               
c  compute saturation specific humidity(dimensionless) from                     
c  temperature t (deg k) and pressure (cb)                                      
c  conversion to specific humidity follows from definition                      
c                                                                               
      do 7 k=1,kdimi                                                            
      expon = 7.50*(t(ix,k)-t0c_)/((t(ix,k)-t0c_)+237.30)                       
      es = psat_*1.e-2 * 10.0**expon                                            
      qsat(k)=rd_/rv_*es/(prs(k)*10.0-(1.0-rd_/rv_)*es)                         
    7 continue                                                                  
c                                                                               
c  limit moisture in lowest levh layers--rh le 1. but ge .15                    
c          the latter avoids the negative q problem...                          
c                                                                               
      do 10 k=1,kdimi                                                          
      qout(k) = qin(ix,k)                                                       
      rh = qin(ix,k) / qsat(k)                                                  
   10 continue                                                                  
c                                                                               
c  obtain 8.cb (k=1) and 1.cb (k=8) valu by horiz interpolation                 
c                                                                               
      jdx = (90.-dlat)/10.+1.                                                   
      dx = (90.-dlat)/10.+1.-jdx                                                
      q8  = qq(1,jdx)*(1.-dx)+qq(1,jdx+1)*dx                                    
      q1  = qq(8,jdx)*(1.-dx)+qq(8,jdx+1)*dx                                    
      ldry = kdimi + 1                                                         
c hmhj 1                                                                        
      if( ldry .le. kdimi ) then                                                
c                                                                               
      xlnqm =  log(q8)                                                          
c                                                                               
c  extrapolate moisture to min valu(q8) at pressure pmin(80 mb)                 
c   use exponential decrease from layer levh----                                
c   i.e.  q=q(levh)*(p/p(levh)) ** xlamb                                        
c   where xlamb is computed to fit q8 at pmin and q,p at levh.                  
c                                                                               
      qqout = qin(ix,kdimi)                                                    
      rh = qqout / qsat(kdimi)                                                 
      if (rh.le..15) qqout = .15 * qsat(kdimi)                                 
      if (rh.gt. 1.0) qqout = qsat(kdimi)                                      
      xlnpq =  log(prs(kdimi))                                                 
      xlnqq =  log(qqout)                                                       
      xlamb = (xlnqm-xlnqq) / (xlnpm-xlnpq)                                     
      do 15 k=ldry,kdimi                                                        
      if (prs(k).ge.pmin) go to 21                                              
      if (prs(k).ge.pmm) go to 22                                               
c                                                                               
c  above pmm(10 mb) use constant value from table                               
c                                                                               
      qout(k) = q1                                                              
      go to 13                                                                  
   21 continue                                                                  
      qout(k) = qout(kdimi)*(prs(k)/prs(kdimi)) ** xlamb                      
      go to 13                                                                  
   22 continue                                                                  
c                                                                               
c  above 8 cb so complete linear interp from table..                            
c                                                                               
      kdy = 9. - prs(k)                                                         
      dy = 9. - prs(k) - kdy                                                    
      qout(k) = qq(kdy,jdx  )*(1.-dy)*(1.-dx)+qq(kdy+1,jdx+1)*dy*dx             
     1         +qq(kdy,jdx+1)*(1.-dy)*dx+qq(kdy+1,jdx)*dy*(1.-dx)               
   13 rh = qout(k) / qsat(k)                                                    
      if(rh.gt.1.) qout(k)=qsat(k)                                              
   15 continue                                                                  
c hmhj 1                                                                        
      endif                                                                     
c                                                                               
c  store extrapolated moisture                                                  
c                                                                               
      do 20 k=1,kdimi                                                           
      chgrp(ix,k) = qout(k)                                                     
   20 continue                                                                  
c ---- temp. is restored to virtual temp.                                       
      do 25 k=1,kdimi                                                          
        if (qin(ix,k).le.0.) go to 25                                           
        t(ix,k) = t(ix,k) * (1. + rv_/rd_ * qin(ix,k))                          
     1                   /(1. + qin(ix,k))                                      
   25 continue                                                                  
c                                                                               
   50 continue                                                                  
c                                                                               
 1000 continue                                                                  
c                                                                               
      return                                                                    
      end                                                                       
