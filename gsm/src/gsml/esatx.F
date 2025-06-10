      subroutine esatx(tt, es, desdt)                                           
                                                                                
      implicit none                                                             
      real tt, es, desdt                                                        
                                                                                
      real        airmw, h2omw, lice, runiv, rgas                               
      parameter ( airmw  = 28.97      )                                         
      parameter ( h2omw  = 18.01      )                                         
      parameter ( lice   = 2.834e6    )                                         
      parameter ( runiv  = 8314.3     )                                         
      parameter ( rgas   = runiv/airmw)                                         
                                                                                
      real tmin, tmax, zeroc, tmix                                              
      real tstarr1, tstarr2, tstarr3, tstarr4                                   
      parameter ( zeroc   = 273.16 )                                            
      parameter ( tmin    = -95. )                                              
      parameter ( tstarr1 = -75. )                                              
      parameter ( tstarr2 = -65. )                                              
      parameter ( tstarr3 = -50. )                                              
      parameter ( tstarr4 = -40. )                                              
      parameter ( tmix    = -20. )                                              
      parameter ( tmax    = +60. )                                              
                                                                                
                                                                                
      real        esfac, erfac                                                  
      parameter ( esfac = h2omw/airmw       )                                   
      parameter ( erfac = (1.0-esfac)/esfac )                                   
                                                                                
                                                                                
      real       b0, b1, b2, b3, b4, b5, b6                                     
      parameter (b6= 6.136820929e-11*esfac)                                     
      parameter (b5= 2.034080948e-8 *esfac)                                     
      parameter (b4= 3.031240396e-6 *esfac)                                     
      parameter (b3= 2.650648471e-4 *esfac)                                     
      parameter (b2= 1.428945805e-2 *esfac)                                     
      parameter (b1= 4.436518521e-1 *esfac)                                     
      parameter (b0= 6.107799961e+0 *esfac)                                     
                                                                                
      real bi0, bi1, bi2, bi3, bi4, bi5, bi6                                    
      parameter (bi6= 1.838826904e-10*esfac)                                    
      parameter (bi5= 4.838803174e-8 *esfac)                                    
      parameter (bi4= 5.824720280e-6 *esfac)                                    
      parameter (bi3= 4.176223716e-4 *esfac)                                    
      parameter (bi2= 1.886013408e-2 *esfac)                                    
      parameter (bi1= 5.034698970e-1 *esfac)                                    
      parameter (bi0= 6.109177956e+0 *esfac)                                    
                                                                                
      real s10, s11, s12, s13, s14, s15, s16                                    
      parameter (s16= 0.516000335e-11*esfac)                                    
      parameter (s15= 0.276961083e-8 *esfac)                                    
      parameter (s14= 0.623439266e-6 *esfac)                                    
      parameter (s13= 0.754129933e-4 *esfac)                                    
      parameter (s12= 0.517609116e-2 *esfac)                                    
      parameter (s11= 0.191372282e+0 *esfac)                                    
      parameter (s10= 0.298152339e+1 *esfac)                                    
                                                                                
      real s20, s21, s22, s23, s24, s25, s26                                    
      parameter (s26= 0.314296723e-10*esfac)                                    
      parameter (s25= 0.132243858e-7 *esfac)                                    
      parameter (s24= 0.236279781e-5 *esfac)                                    
      parameter (s23= 0.230325039e-3 *esfac)                                    
      parameter (s22= 0.129690326e-1 *esfac)                                    
      parameter (s21= 0.401390832e+0 *esfac)                                    
      parameter (s20= 0.535098336e+1 *esfac)                                    
                                                                                
      real       c1, c2, c3, c4, c5, c6                                         
      parameter (c1= b1   )                                                     
      parameter (c2= b2*2.)                                                     
      parameter (c3= b3*3.)                                                     
      parameter (c4= b4*4.)                                                     
      parameter (c5= b5*5.)                                                     
      parameter (c6= b6*6.)                                                     
                                                                                
      real       ci1, ci2, ci3, ci4, ci5, ci6                                   
      parameter (ci1= bi1   )                                                   
      parameter (ci2= bi2*2.)                                                   
      parameter (ci3= bi3*3.)                                                   
      parameter (ci4= bi4*4.)                                                   
      parameter (ci5= bi5*5.)                                                   
      parameter (ci6= bi6*6.)                                                   
                                                                                
      real       d11, d12, d13, d14, d15, d16                                   
      parameter (d11= s11   )                                                   
      parameter (d12= s12*2.)                                                   
      parameter (d13= s13*3.)                                                   
      parameter (d14= s14*4.)                                                   
      parameter (d15= s15*5.)                                                   
      parameter (d16= s16*6.)                                                   
                                                                                
      real       d21, d22, d23, d24, d25, d26                                   
      parameter (d21= s21   )                                                   
      parameter (d22= s22*2.)                                                   
      parameter (d23= s23*3.)                                                   
      parameter (d24= s24*4.)                                                   
      parameter (d25= s25*5.)                                                   
      parameter (d26= s26*6.)                                                   
                                                                                
      real t, w, qx, dqx                                                        
      t = amax1(amin1(tt-zeroc,tmax),tmin)                                      
      dqx = 0.                                                                  
      qx = 0.                                                                   
      if(t.gt.0.) then                                                          
       qx = (t*(t*(t*(t*(t*(t*b6+b5)+b4)+b3)+b2)+b1)+b0)                        
       dqx = (t*(t*(t*(t*(t*c6+c5)+c4)+c3)+c2)+c1)                              
      elseif(t.lt.tstarr1) then                                                 
       qx = (t*(t*(t*(t*(t*(t*s16+s15)+s14)+s13)+s12)+s11)+s10)                 
       dqx = (t*(t*(t*(t*(t*d16+d15)+d14)+d13)+d12)+d11)                        
      elseif(t.lt.tstarr2) then                                                 
       w = (tstarr2 - t)/(tstarr2-tstarr1)                                      
       qx =     w *(t*(t*(t*(t*(t*(t*s16+s15)+s14)+s13)+s12)+s11)+s10)          
     *    + (1.-w)*(t*(t*(t*(t*(t*(t*s26+s25)+s24)+s23)+s22)+s21)+s20)          
       dqx =     w *(t*(t*(t*(t*(t*d16+d15)+d14)+d13)+d12)+d11)                 
     *     + (1.-w)*(t*(t*(t*(t*(t*d26+d25)+d24)+d23)+d22)+d21)                 
      elseif(t.lt.tstarr3) then                                                 
       qx = (t*(t*(t*(t*(t*(t*s26+s25)+s24)+s23)+s22)+s21)+s20)                 
       dqx = (t*(t*(t*(t*(t*d26+d25)+d24)+d23)+d22)+d21)                        
      elseif(t.lt.tstarr4) then                                                 
       w = (tstarr4 - t)/(tstarr4-tstarr3)                                      
       qx =     w *(t*(t*(t*(t*(t*(t*s26+s25)+s24)+s23)+s22)+s21)+s20)          
     *    + (1.-w)*(t*(t*(t*(t*(t*(t*bi6+bi5)+bi4)+bi3)+bi2)+bi1)+bi0)          
       dqx =     w *(t*(t*(t*(t*(t*d26+d25)+d24)+d23)+d22)+d21)                 
     *     + (1.-w)*(t*(t*(t*(t*(t*ci6+ci5)+ci4)+ci3)+ci2)+ci1)                 
      elseif(t.lt.tmix) then                                                    
       qx = (t*(t*(t*(t*(t*(t*bi6+bi5)+bi4)+bi3)+bi2)+bi1)+bi0)                 
       dqx = (t*(t*(t*(t*(t*ci6+ci5)+ci4)+ci3)+ci2)+ci1)                        
      else                                                                      
       w = (tmix - t)/tmix                                                      
       qx =     w *(t*(t*(t*(t*(t*(t*b6+b5)+b4)+b3)+b2)+b1)+b0)                 
     &    + (1.-w)*(t*(t*(t*(t*(t*(t*bi6+bi5)+bi4)+bi3)+bi2)+bi1)+bi0)          
       dqx =     w *(t*(t*(t*(t*(t*c6+c5)+c4)+c3)+c2)+c1)                       
     &      + (1.-w)*(t*(t*(t*(t*(t*ci6+ci5)+ci4)+ci3)+ci2)+ci1)                
      endif                                                                     
      es    = qx                                                                
      desdt = dqx                                                               
                                                                                
      return                                                                    
      end                                                                       
