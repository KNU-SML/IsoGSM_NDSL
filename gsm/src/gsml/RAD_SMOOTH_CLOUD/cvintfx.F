      subroutine cvintfx(cvin,cvtin,cvbin,iin,jtwidl,jin,
     1                  cvout,cvtout,cvbout,iout,jpout,jout,                    
     3                  xx,wgt,tt,bb,sum,nn,                                    
     4                  ltwidl,latrd1,latinb)                                   
c--   *****************************************************************         
c     *  code bilinearly interpolates cld amt between gaussian grids--*         
c     *  clone of ggintp for interpolation of convective cld amt (cv).*         
c     *    special interp procedure for tops(cvt) and bots(cvb)...    *         
c-    *  j = 1 is just belo n.pole, i = 1 is greenwich (then go east).*         
c     * iin,jin are i,j dimensions of input grid--iout,jout for output*         
c     * jin2,jout2=jin/2,jout/2                                       *         
c     *                                     --k.campana - june 1988   *         
c--   *****************************************************************         
      dimension cvin(iin,jtwidl),cvtin(iin,jtwidl),cvbin(iin,jtwidl)            
      dimension cvout(iout,jpout)                                               
      dimension cvtout(iout,jpout),cvbout(iout,jpout)                           
      dimension xx(iout,4),wgt(iout,4),tt(iout,4),bb(iout,4),sum(iout,4)        
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
c===>    if output lat is poleward of input lat=1 ,then simpl average           
c          (small region and cld amt wouldn t extrapolate well)                 
       call cintpx(iii,jbb,jjj,iiiout,inslat,wgtlat,                   
     2            cvin,cvtin,cvbin,cvout(1,latout),                             
     3            cvtout(1,latout),cvbout(1,latout),                            
     3            xx,wgt,tt,bb,sum,nn,lbb,lr1)                                  
   50 continue                                                                  
ck100 format(1h ,' row =',i5,'  lat =',e15.5)                                   
      return                                                                    
      end                                                                       
