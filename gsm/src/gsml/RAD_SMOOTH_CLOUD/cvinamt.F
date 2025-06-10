      subroutine cvinamt(cvin,iin,jtwidl,jin,cvout,iout,jpout,jout, 
     3                  xx,wgt,sum,nn,ltwidl,latrd1,latinb)                                   
c--   *****************************************************************         
c     *  code bilinearly interpolates cld amt between gaussian grids--*         
c     *  clone of cvintf without the cloud top/base interpolation    *         
c-    *  j = 1 is just belo n.pole, i = 1 is greenwich (then go east).*         
c     * iin,jin are i,j dimensions of input grid--iout,jout for output*         
c     * jin2,jout2=jin/2,jout/2                                       *         
c     *                           campana+katz+campana(again) nov94   *         
c--   *****************************************************************         
      dimension cvin(iin,jtwidl)                                                
      dimension cvout(iout,jpout)                                               
      dimension xx(iout,4),wgt(iout,4),sum(iout,4)        
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
         wgtlat=1.0
       endif
ccc     print 100,lat,xlat                                                      
c===>    if output lat is poleward of input lat=1 ,then simpl average           
c          (small region and cld amt wouldn t extrapolate well)                 
       call cinpamt(iii,jbb,jjj,iiiout,inslat,wgtlat,                   
     2            cvin,cvout(1,latout),xx,wgt,sum,nn,lbb,lr1)                                  
   50 continue                                                                  
ck100 format(1h ,' row =',i5,'  lat =',e15.5)                                   
      return                                                                    
      end                                                                       
