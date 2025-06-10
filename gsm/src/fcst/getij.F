#define LATG latg
#define LATGP latgp
#define LONF lonf
#define LONFP lonfp
      subroutine getij (xlat,xlon,slmsk,blat,dx,                                
     1                  ils,idm,jdm,ki,kj)                                      
#include <paramodel.h>
#ifdef KEN
#ifdef MP
#define LATGS latgp_
#define LONFS lonfp_
#else
#define LATGS latg_
#define LONFS lonf_
#endif
      dimension blat(jdm),dist(4),keni(4),kenj(4)                               
      dimension ipsort(4),ipskp(4)                                              
      dimension slmsk(LONFS,LATGS)                                              
      jdm2 = jdm / 2                                                            
      jdmp1 = jdm2 + 1                                                          
      if (abs(xlat).gt.blat(1)) go to 70                                        
c----    get upper left gaussian point (ia,ja) on gridbox                       
c          surrounding the input lat/lon point........                          
        ia = xlon/dx + 1                                                        
        ib = ia + 1                                                             
        if (ia.ge.idm) ib = 1                                                   
        xi = xlon/dx + 1. - ia                                                  
        do 10 jak=2,jdmp1                                                       
         jb = jak - 1                                                           
         if(abs(xlat).gt.blat(jak)) go to 15                                    
  10    continue                                                                
  15    continue                                                                
        ja = jb                                                                 
c----   normalized distance from upper lat to gaussian lat                      
        xj = (blat(ja) - abs(xlat)) / (blat(ja)-blat(ja+1))                     
c       xout(i,lat)   = (1-xi)* xj   *xin(ia,ja+1) +                            
cc   1                 xi *  xj  *xin(ib,ja+1) +                                
ccc  2              (1-xi)*(1-xj)*xin(ia  ,ja  ) +                              
cccc 3               xi   *(1-xj)*xin(ib,ja  )                                  
c----    southern hemisphere                                                    
        if (xlat.lt.0.) then                                                    
         ja = jdm - ja                                                          
         xj = 1. - xj                                                           
        end if                                                                  
cccc    xout(i,jout+1-lat)=(1-xi)* xj   *xin(ia,ja+1) +                         
ccc  1                 xi *  xj  *xin(ib,ja+1) +                                
cc   2              (1-xi)*(1-xj)*xin(ia  ,ja  ) +                              
c    3               xi   *(1-xj)*xin(ib,ja  )                                  
c...     upper left point                                                       
      dist(1)= sqrt(xi**2+xj**2)                                                
      keni(1)= ia                                                               
      kenj(1)= ja                                                               
c...     upper right                                                            
      dist(2)= sqrt((1-xi)**2+xj**2)                                            
      keni(2)= ib                                                               
      kenj(2)= ja                                                               
c...     lower right                                                            
      dist(3)= sqrt((1-xi)**2+(1-xj)**2)                                        
      keni(3)= ib                                                               
      kenj(3)= ja+1                                                             
c...     lower left                                                             
      dist(4)= sqrt(xi**2+(1-xj)**2)                                            
      keni(4)= ia                                                               
      kenj(4)= ja+1                                                             
c---     now sort the distances (by index, shortest first)                      
      npt = 4                                                                   
      npt1 = npt + 1                                                            
      do 20 kd=1,npt                                                            
        ipskp(kd) = 0                                                           
   20 continue                                                                  
      do 40 kd=1,npt                                                            
        dd=100.                                                                 
c---     find shortest dist of the remaining unsorted data...                   
        do 25 kk=1,npt                                                          
         if(ipskp (kk).gt.0) go to 25                                           
          if(dist(kk).lt.dd) then                                               
            dd = dist(kk)                                                       
            jx = kk                                                             
          end if                                                                
   25   continue                                                                
c---   store sorted index                                                       
        ipsort(kd) = jx                                                         
        ipskp (jx) = 1                                                          
   40 continue                                                                  
#ifndef NOPRINT
      print 102,(dist(kk),kk=1,npt),(ipsort(kk),kk=1,npt)                       
  102 format(1h ,' distances=',4f8.4,' sorted indices=',4i4)                    
#endif
c---   end of distance sort                                                     
      xils = ils                                                                
      if (ils.lt.0) then                                                        
       ki = keni(ipsort(1))                                                     
       kj = kenj(ipsort(1))                                                     
       return                                                                   
      end if                                                                    
      if (ils.eq.0) then                                                        
c....     find nearest sea point                                                
       do 45 kd=1,npt                                                           
        ii = keni(ipsort(kd))                                                   
        jj = kenj(ipsort(kd))                                                   
        if (slmsk(ii,jj).le.xils) go to 46                                      
   45  continue                                                                 
c....     no sea points so default to nearest point                             
#ifndef NOPRINT
       print 49,xlat,xlon                                                       
   49  format(1h ,' asked for sea point but can t find one for',                
     1            ' lat lon=',2f9.2,'..so default to nearest')                  
#endif
       ki = keni(ipsort(1))                                                     
       kj = kenj(ipsort(1))                                                     
       return                                                                   
   46  continue                                                                 
       ki = ii                                                                  
       kj = jj                                                                  
       return                                                                   
      end if                                                                    
      if (ils.eq.1) then                                                        
c....     find nearest land/ice point                                           
       do 55 kd=1,npt                                                           
        ii = keni(ipsort(kd))                                                   
        jj = kenj(ipsort(kd))                                                   
        if (slmsk(ii,jj).ge.xils) go to 56                                      
   55  continue                                                                 
c....     no land points so default to nearest point                            
#ifndef NOPRINT
       print 59,xlat,xlon                                                       
   59  format(1h ,' asked for land sea point but can t find one for',           
     1            ' lat lon=',2f9.2,'..so default to nearest')                  
#endif
       ki = keni(ipsort(1))                                                     
       kj = kenj(ipsort(1))                                                     
       return                                                                   
   56  continue                                                                 
       ki = ii                                                                  
       kj = jj                                                                  
       return                                                                   
      end if                                                                    
c...   outside limit of gaussian polar rows so just take nearest                
c        point without regard to land and sea                                   
   70 continue                                                                  
      ia = xlon/dx + 1                                                          
      ib = ia + 1                                                               
      if (ia.ge.idm) ib = 1                                                     
      xi = xlon/dx + 1. - ia                                                    
      ja = 1                                                                    
      if (xlat.lt.0.) ja = jdm                                                  
      if (xi.gt.0.5) then                                                       
       ki = ib                                                                  
       kj = ja                                                                  
      else                                                                      
       ki = ia                                                                  
       kj = ja                                                                  
      end if                                                                    
#endif
      return                                                                    
      end                                                                       
