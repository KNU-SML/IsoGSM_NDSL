      subroutine coeint(rat,ir)                                                 
c **********************************************************************        
c                                                                               
c                                                                               
c            the transmission function between p1 and p2 is assumed to          
c       the  functional form                                                    
c                     tau(p1,p2)= 1.0-sqrt(c*log(1.0+x*path)),                  
c               where                                                           
c                     path(p1,p2)=((p1-p2)**2)*(p1+p2+core)/                    
c                                 (eta*(p1+p2+core)+(p1-p2))                    
c                                                                               
c                                                                               
c        the parameters c and x are functions of p2, and are to be deter        
c        while core is a prespecified number.eta is a function of the th        
c        product (cx);it is obtaited iteratively. the derivation of all         
c        values will be explained in a forthcoming paper.                       
c            subroutine coeint determines c(i) and x(i) by using the act        
c        values of tau(p(i-2),p(i)) and tau(p(i-1),p(i)) and the previou        
c        iteration value of eta.                                                
c             define:                                                           
c                patha=path(p(i),p(i-2),core,eta)                               
c                pathb=path(p(i),p(i-1),core,eta);                              
c        then                                                                   
c                r=(1-tau(p(i),p(i-2)))/(1-tau(p(i),p(i-1)))                    
c                 = sqrt(log(1+x*patha)/log(1+x*pathb)),                        
c        so that                                                                
c                r**2= log(1+x*patha)/log(1+x*pathb).                           
c        this equation can be solved by newton s method for x and then t        
c        result used to find c. this is repeated for each value of i gre        
c        than 2 to give the arrays x(i) and c(i).                               
c             newton s method for solving the equation                          
c                 f(x)=0                                                        
c        makes use of the loop xnew= xold-f(xold)/fprime(xold).                 
c        this is iterated 20 times, which is probably excessive.                
c        the first guess for eta is 3.2e-4*exp(-p(i)/1000),which has            
c        been found to be fairly realistic by experiment; we iterate 5 t        
c        (again,probably excessively) to obtain the values for c,x,eta t        
c        used for interpolation.                                                
c           there are several possible pitfalls:                                
c              1) in the course of iteration, x may reach a value which         
c                 1+x*patha negative; in this case the iteration is stop        
c                 and an error message is printed out.                          
c              2) even if (1) does not occur, it is still possible that         
c                 be negative and large enough to make 1+x*path(p(i),0,c        
c                 negative. this is checked for in a final loop, and if         
c                 a warning is printed out.                                     
c                                                                               
c  *********************************************************************        
c....                                                                           
c     implicit double precision (a-h,o-z)                                       
      common/press/pa(109)                                                      
      real rat,sinv                                                             
      real pa,core,transa,path,uexp,sexp,eta,sexpv                              
      real pa2                                                                  
      common/tran/ transa(109,109)                                              
      common/coefs/xa(109),ca(109),eta(109),sexpv(109),core,uexp,sexp           
      dimension path0(109),etap(109),xap(109),cap(109)                          
      dimension sinv(4)                                                         
      data sinv/2.74992,2.12731,4.38111,0.0832926/                              
cnov89   dimension sinv(3)                                                      
cnov89   data sinv/2.74992,2.12731,4.38111/                                     
co222  old code used 2.7528 rather than 2.74992 ---k.a.c. october 1988          
co222   when 2.7528 was used,we exactly reproduced the mrf co2 arrays           
      core=5.000                                                                
      uexp=0.90                                                                 
      p0=0.7                                                                    
      do 902 i=1,109                                                            
      pa2=pa(i)*pa(i)                                                           
      sexpv(i)=.505+2.0e-5*pa(i)+.035*(pa2-.25)/(pa2+.25)                       
902   continue                                                                  
      do 900 i=1,109                                                            
      eta(i)=3.2e-4*exp(-pa(i)/500.)                                            
      etap(i)=eta(i)                                                            
900   continue                                                                  
      do 1200 np=1,10                                                           
      do 1000 i=3,109                                                           
      sexp=sexpv(i)                                                             
      r=(1.0-transa(i,i-2))/(1.0-transa(i,i-1))                                 
      rexp=r**(uexp/sexp)                                                       
      patha=(path(pa(i),pa(i-2),core,eta(i)))**uexp                             
      pathb=(path(pa(i),pa(i-1),core,eta(i)))**uexp                             
      xx=2.0*(pathb*rexp-patha)/(pathb*pathb*rexp-patha*patha)                  
      do 1010 ll=1,20                                                           
      f1=log(1.0+xx*patha)                                                      
      f2=log(1.0+xx*pathb)                                                      
      f=f1/f2-rexp                                                              
      fprime=(f2*patha/(1.0+xx*patha)-f1*pathb/(1.0+xx*pathb))/                 
     1    (f2*f2)                                                               
      xx=xx-f/fprime                                                            
      check=1.0+xx*patha                                                        
      if (check) 1020,1020,1025                                                 
 1020 continue                                                                  
      write (6,360) i,ll,check                                                  
360   format (' error,i=',i3,'ll=',i3,'check=',f20.10)                          
      stop                                                                      
1025  continue                                                                  
 1010 continue                                                                  
      ca(i)=(1.0-transa(i,i-2))**(uexp/sexp)/                                   
     1 (log(1.0+xx*patha)+1.0e-20)                                              
      xa(i)=xx                                                                  
1000  continue                                                                  
      xa(2)=xa(3)                                                               
      xa(1)=xa(3)                                                               
      ca(2)=ca(3)                                                               
      ca(1)=ca(3)                                                               
      do 1100 i=3,109                                                           
      path0(i)=(path(pa(i),0.,core,eta(i)))**uexp                               
      path0(i)=1.0+xa(i)*path0(i)                                               
      if (path0(i).lt.0.) write (6,361) i,path0(i),xa(i)                        
1100  continue                                                                  
      do 1035 i=1,109                                                           
      sexp=sexpv(i)                                                             
      etap(i)=eta(i)                                                            
      eta(i)=(sinv(ir)/rat)**(1./sexp)*                                         
     1  (ca(i)*xa(i))**(1./uexp)                                                
1035  continue                                                                  
c                                                                               
c     the eta formulation is detailed in schwarzkopf and fels(1985).            
c        the quantity sinv=(g*deltanu)/(rco2*d*s)                               
c      in cgs units,with d,the diffusicity factor=2, and                        
c      s,the sum of co2 line strengths over the 15um co2 band                   
c       also,the denominator is multiplied by                                   
c      1000 to permit use of mb units for pressure.                             
c        s is actually weighted by b(250) at 10 cm-1 wide intervals,in          
c      order to be consistent with the methods used to obtain the lbl           
c      1-band consolidated trancmission functions.                              
c      for the 490-850 interval (deltanu=360,ir=1) sinv=2.74992.                
c      (slightly different from 2.7528 used in earlier versions)                
c      for the 490-670 interval (ir=2) sinv=2.12731                             
c      for the 670-850 interval (ir=3) sinv=4.38111                             
c      for the 2270-2380 interval (ir=4) sinv=0.0832926                         
c      sinv has been obtained using the 1982 afgl catalog for co2               
c        rat is the actual co2 mixing ratio in units of 330 ppmv,               
c      letting use of this formulation for any co2 concentration.               
c                                                                               
c     write (6,366) (np,i,ca(i),xa(i),eta(i),sexpv(i),i=1,109)                  
c366   format (2i4,4e20.12)                                                     
1200  continue                                                                  
 361  format (' **warning:** 1+xa*path(pa(i),0) is negative,i= ',i3,/           
     1 20x,'path0(i)=',f16.6,' xa(i)=',f16.6)                                   
      return                                                                    
      end                                                                       
