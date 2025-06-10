      subroutine prof5(eq,ee,ud)                                        
!
!***********************************************************************
!*****    gaussian type mixing profile....******************************
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  this subroutine integrates the area under the curve in the gaussian  
!  distribution...the numerical approximation to the integral is taken from
!  "handbook of mathematical functions with formulas, graphs and mathematics tables"
!  ed. by abramowitz and stegun, natl bureau of standards applied
!  mathematics series.  june, 1964., may, 1968.                         
!                                     jack kain                         
!                                     7/6/89                            
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
      real,         intent(in   )   :: eq
      real,         intent(inout)   :: ee,ud
      real ::       sqrt2p,a1,a2,a3,p,sigma,fe,x,y,ey,e45,t1,t2,c1,c2

      data sqrt2p,a1,a2,a3,p,sigma,fe/2.506628,0.4361836,-0.1201676, 
     1     0.9372980,0.33267,0.166666667,0.202765151/   
      x=(eq-0.5)/sigma                                                  
      y=6.*eq-3.                                                        
      ey=exp(y*y/(-2))                                                  
      e45=exp(-4.5)                                                     
      t2=1./(1.+p*abs(y))                                               
      t1=0.500498                                                       
      c1=a1*t1+a2*t1*t1+a3*t1*t1*t1                                     
      c2=a1*t2+a2*t2*t2+a3*t2*t2*t2                                     
      if(y.ge.0.)then                                                   
        ee=sigma*(0.5*(sqrt2p-e45*c1-ey*c2)+sigma*(e45-ey))-e45*eq*eq/2.
        ud=sigma*(0.5*(ey*c2-e45*c1)+sigma*(e45-ey))-e45*(0.5+eq*eq/2.- 
     1     eq)                                                          
      else                                                              
        ee=sigma*(0.5*(ey*c2-e45*c1)+sigma*(e45-ey))-e45*eq*eq/2.       
        ud=sigma*(0.5*(sqrt2p-e45*c1-ey*c2)+sigma*(e45-ey))-e45*(0.5+eq*
     1     eq/2.-eq)                                                    
      endif                                                             
      ee=ee/fe                                                          
      ud=ud/fe                                                          

      end subroutine prof5
