      subroutine ggavet(cttin,iin,jtwidl,jin,cttout,iout,jpout,jout,            
     2                  tt,sum,nn,ltwidl,latrd1,latinb)                         
c--   *****************************************************************         
c     * put cloud top temperature onto fcst model grid......          *         
c     *    only average those points which have cld (ie temp nonzero) *         
c-    *  j = 1 is just belo n.pole, i = 1 is greenwich (then go east).*         
c     * iin,jin are i,j dimensions of input grid--iout,jout for output*         
c     * jin2,jout2=jin/2,jout/2                                       *         
c     *                                     --k.campana - august 91   *         
c--   *****************************************************************         
      dimension cttin(iin,jtwidl)                                               
      dimension cttout(iout,jpout)                                              
      dimension tt(iout,4),sum(iout)                                            
      dimension nn(iout)                                                        
      iii = iin                                                                 
      jbb = jtwidl                                                              
      jjj = jin                                                                 
      iiiout = iout                                                             
      lbb = ltwidl                                                              
      lr1 = latrd1                                                              
      do 50 latout=1,jpout                                                      
       lat=latout+latinb-1                                                      
       if(lat.eq.1) then
         inslat=-1  
         wgtlat=0.
       else
         inslat=lat-1
         wgtlat=1.
       endif
ccc     print 100,lat,xlat                                                      
c===>    if output lat is poleward of input lat=1 ,then simpl average           
c          (small region and cld amt wouldn t extrapolate well)                 
       call gintp(iii,jbb,jjj,iiiout,                                           
     1            inslat,wgtlat,                          
     2            cttin,cttout(1,latout),tt,sum,nn,lbb,lr1)                     
   50 continue                                                                  
ck100 format(1h ,' row =',i5,'  lat =',e15.5)                                   
      return                                                                    
      end                                                                       
