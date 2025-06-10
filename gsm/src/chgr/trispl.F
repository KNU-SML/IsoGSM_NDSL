      subroutine trispl(n,p,q,y)                                                
      save                                                                      
c-----------------------------------------------------                          
c    subroutine calculates weights for cubic spline interpolation               
c     n...number ofc     n...number of observations                             
c     y...function y(p) to be interpolated                                      
c     p...coor (independent variables)                                          
c     q,ovh,sh...contain weights                                                
c----------------------------------------------------                           
c                                                                               
      common/spl/ ovh(100),sh(100),theta
      common/spli/ iflag,jflag
      dimension c(100),b(100),p(n),q(n),d(100),y(n),                            
     1 tchs(100),tb(100)                                                        
      nu=n-1                                                                    
      if (iflag.ne.0) go to 99                                                  
      do 9 i=2,n                                                                
      im=i-1                                                                    
      him=p(i)-p(im)                                                            
      if(him .eq. 0.) him=0.001                                                 
      ovh(im)=1/him                                                             
      th=theta*him                                                              
       sh(im)=1/sinh(th)                                                        
 9    tchs(im)=(theta*cosh(th))*sh(im)-ovh(im)                                  
      do 8 i=2,nu                                                               
      im=i-1                                                                    
      c(im)=ovh(im)-theta*sh(im)                                                
 8    b(i)=tchs(im)+tchs(i)                                                     
      c(n-1)=ovh(n-1)-theta*sh(n-1)                                             
      b(1)=tchs(1)                                                              
      b(n)=tchs(n-1)                                                            
      tb(1)=c(1)/b(1)                                                           
      do 7 i=2,nu                                                               
 7    tb(i)=c(i)/(b(i)-c(i-1)*tb(i-1))                                          
      d(1)=0.                                                                   
      d(n)=0.                                                                   
      iflag=1                                                                   
 99   yp=(y(2)-y(1))*ovh(1)                                                     
      do 6 i=2,nu                                                               
      ynow=(y(i+1)-y(i))*ovh(i)                                                 
      d(i)=ynow-yp                                                              
 6    yp=ynow                                                                   
      q(1)=d(1)/b(1)                                                            
      do 5 i=2,n                                                                
      im=i-1                                                                    
 5     q(i)=(d(i)-c(im)*q(im))/(b(i)-c(im)*tb(im))                              
      do 4 i=1,nu                                                               
      ii=n-i                                                                    
 4    q(ii)=q(ii)-tb(ii)*q(ii+1)                                                
      return                                                                    
      end                                                                       
