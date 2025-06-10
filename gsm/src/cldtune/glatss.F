      subroutine glatss(khalf,colrad)                                           
      save                                                                      
      dimension colrad(khalf)                                                   
      eps=1.e-12                                                                
c     print 101                                                                 
101   format ('0 i   colat   colrad     wgt', 12x, 'wgtcs',                     
     1 10x, 'iter  res')                                                        
      si=1.e+00                                                                 
      k2=2*khalf                                                                
      rk2=k2                                                                    
      scale=2.e+00/(rk2**2)                                                     
      k1=k2-1                                                                   
      pi=atan(si)*4.e+00                                                        
      dradz=pi/360.e+00                                                         
      rad=0.e+00                                                                
      do 1000 k=1,khalf                                                         
      iter=0                                                                    
      drad=dradz                                                                
1     call poly(k2,rad,p2)                                                      
2     p1 =p2                                                                    
      iter=iter+1                                                               
      rad=rad+drad                                                              
      call poly(k2,rad,p2)                                                      
      if(sign(si,p1).eq.sign(si,p2))go to 2                                     
      if(drad.lt.eps)go to 3                                                    
      rad=rad-drad                                                              
      drad=drad*0.25e+00                                                        
      go to 1                                                                   
3     continue                                                                  
      colrad(k)=rad                                                             
      phi=rad   *180.e+00/pi                                                    
      call poly(k1,rad,p1)                                                      
      x=cos(rad)                                                                
      w=scale*(1.e+00-x*x)/(p1*p1)                                              
      sn=sin(rad)                                                               
      w=w/(sn*sn)                                                               
      rc=1./(sn*sn)                                                             
      call poly(k2,rad,p1)                                                      
      prphi = phi                                                               
      prcol = colrad(k)                                                         
c     print 102, k, prphi, prcol, wgt(k), wgtcs(k), iter, p1                    
102   format(1h ,i2,2x,f6.2,2x,f10.7,2x,e13.7,2x,e13.7,2x,i4,2x,e13.7)          
1000  continue                                                                  
      return                                                                    
      end                                                                       
