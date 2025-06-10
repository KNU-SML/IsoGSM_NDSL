#define LEVS levs
      subroutine co2ptz(sgtemp,t41,t42,t43,t44,sglvnu,siglnu,lread)             
#include <paramodel.h>
      parameter (nl=levs_,nlp=nl+1,nlp2=nl+2)                                   
c                                                                               
c **         this program calculates temperatures ,h2o mixing ratios            
c **         and o3 mixing ratios by using an analytical                        
c **         function which approximates                                        
c **         the us standard (1976).  this is                                   
c **         calculated in function 'antemp', which is called by the            
c **         main program.  the form of the analytical function was             
c **         suggested to me in 1971 by richard s. lindzen.                     
c ******************************************************************            
c         code to save stemp,gtemp on data set,bracketed by co222  **           
c             ....k. campana march 88,october 88                                
      dimension sgtemp(nlp,2),t41(nlp2,2),t42(nlp),                             
     1          t43(nlp2,2),t44(nlp)                                            
      dimension sglvnu(nlp),siglnu(nl)                                          
c ******************************************************************            
c                                                                               
c*****this version is only usable for 1976 us std atm and obtains               
c     quantities for co2 interpolation and insertion into opera-                
c     tional radiation codes                                                    
c                                                                               
      character*20 profil                                                       
      dimension press(nlp),temp(nlp),alt(nlp),wmix(nlp),o3mix(nlp)              
      dimension wmxint(nlp,4),wmxout(nlp2),omxint(nlp,4),omxout(nlp2)           
      dimension pd(nlp2),gtemp(nlp)                                             
      dimension prs(nlp),temps(nlp),prsint(nlp),tmpint(nlp,4),a(nlp,4)          
      dimension prout(nlp2),tmpout(nlp2),tmpflx(nlp2),tmpmid(nlp2)              
c                                                                               
c                                                                               
      data profil/                                                              
     6   'us standard 1976'/                                                    
      data psmax/1013.250/                                                      
c                                                                               
c **         ntype is an integer variable which has the following               
c **        values:    0 =sigma levels are used;   1= skyhi l40 levels          
c **        are used;   2 = skyhi l80 levels are used. default: 0               
c                                                                               
      ntype=0                                                                   
co222 read (*,*) ntype                                                          
    5 nlev=nl                                                                   
      delzap=0.5                                                                
      r=8.31432                                                                 
      g0=9.80665                                                                
      zmass=28.9644                                                             
      aa=6356.766                                                               
         alt(1)=0.0                                                             
         temp(1)=antemp(6,0.0)                                                  
c*******determine the pressures (press)                                         
      pstar=psmax                                                               
c                                                                               
co222 if (ntype.eq.1) call skyp(pstar,pd,gtemp)                                 
co222 if (ntype.eq.2) call sky80p(pstar,pd,gtemp)                               
co222 if (ntype.eq.0) call sigp(pstar,pd,gtemp)                                 
ccc----      call sigp(pstar,pd,gtemp)                                          
      call sigp(pstar,pd,gtemp,t41,t42,t43,t44,sglvnu,siglnu,lread)             
      pd(nlp2)=pstar                                                            
      do 40 n=1,nlp                                                             
      prsint(n)=pd(nlp2+1-n)                                                    
 40   continue                                                                  
c    *** calculate temps for several pressures to do quadrature                 
      do 504 nq=1,4                                                             
      do 505 n=2,nlp                                                            
 505  press(n)=prsint(n)+0.25*(nq-1)*(prsint(n-1)-prsint(n))                    
      press(1)=prsint(1)                                                        
c*********************                                                          
      do 100 n=1,nlev                                                           
c                                                                               
c **         establish computatational levels between user levels at            
c **         intervals of approximately 'delzap' km.                            
c                                                                               
      dlogp=7.0*alog(press(n)/press(n+1))                                       
      nint=dlogp/delzap                                                         
      nint=nint+1                                                               
      znint=nint                                                                
      g=g0                                                                      
      dz=r*dlogp/(7.0*zmass*g*znint)                                            
      ht=alt(n)                                                                 
c                                                                               
c **         calculate height at next user level by means of                    
c **                   runge-kutta integration.                                 
c                                                                               
      do 200 m=1,nint                                                           
      rk1=antemp(6,ht)*dz                                                       
      rk2=antemp(6,ht+0.5*rk1)*dz                                               
      rk3=antemp(6,ht+0.5*rk2)*dz                                               
      rk4=antemp(6,ht+rk3)*dz                                                   
      ht=ht+0.16666667*(rk1+rk2+rk2+rk3+rk3+rk4)                                
  200 continue                                                                  
      alt(n+1)=ht                                                               
      temp(n+1)=antemp(6,ht)                                                    
  100 continue                                                                  
      do 506 n=1,nlp                                                            
      tmpint(n,nq)=temp(n)                                                      
      a(n,nq)=alt(n)                                                            
506   continue                                                                  
      do 523 n=1,nlp                                                            
      call mixrat(6,press(n),tmpint(n,nq),wmxint(n,nq))                         
      wmxint(n,nq)=wmxint(n,nq)*1.0e-6                                          
      call ozone(6,press(n),tmpint(n,nq),omxint(n,nq))                          
523   continue                                                                  
504   continue                                                                  
c    ***apply simpson s rule                                                    
      do 507 n=2,nlp                                                            
      temp(n)=1./12.*(tmpint(n-1,1)+tmpint(n,1)+4.*tmpint(n,2)+                 
     1 2.*tmpint(n,3)+4.*tmpint(n,4))                                           
      wmix(n)=1./12.*(wmxint(n-1,1)+wmxint(n,1)+4.*wmxint(n,2)+                 
     1 2.*wmxint(n,3)+4.*wmxint(n,4))                                           
      o3mix(n)=1./12.*(omxint(n-1,1)+omxint(n,1)+4.*omxint(n,2)+                
     1 2.*omxint(n,3)+4.*omxint(n,4))                                           
507   continue                                                                  
c***output for line-by-line calcs                                               
      tmpout(1)=tmpint(nlp,1)                                                   
      wmxout(1)=wmxint(nlp,1)                                                   
      omxout(1)=omxint(nlp,1)                                                   
      tmpmid(1)=tmpint(nlp,1)                                                   
      do 520 i=2,nlp2                                                           
      tmpout(i)=temp(nlp2+1-i)                                                  
      wmxout(i)=wmix(nlp2+1-i)                                                  
      omxout(i)=o3mix(nlp2+1-i)                                                 
      tmpmid(i)=tmpint(nlp2+1-i,3)                                              
520   continue                                                                  
      do 521 i=1,nlp2                                                           
      prout(i)=pd(i)                                                            
521   continue                                                                  
      do 5221 i=2,nlp2                                                          
      tmpflx(i)=tmpint(nlp2+1-i,1)                                              
5221  continue                                                                  
      tmpflx(1)=tmpflx(2)                                                       
c                                                                               
c **        calculate water mixing ratio using luther program                   
c                                                                               
      do 300 n=1,nlp                                                            
      call mixrat(6,prsint(n),tmpint(n,1),wmix(n))                              
      wmix(n)=wmix(n)*1.0e-06                                                   
  300 continue                                                                  
c                                                                               
c   calculate ozone mixing ratio using schwarzkopf program                      
c                                                                               
      do 400 n=1,nlp                                                            
      call ozone(6,prsint(n),tmpint(n,1),o3mix(n))                              
400   continue                                                                  
c                                                                               
      write (6,101) profil                                                      
  101 format (1x,a20)                                                           
      write (6,201)                                                             
  201 format(5x,'   height    temperature     pressure     r(h2o)     ',        
     1' r(o3)')                                                                 
      write (6,202)  a(1,1),tmpint(1,1),prsint(1),wmix(1),o3mix(1)              
      do 210 n=2,nlp                                                            
      write (6,203)  tmpout(nlp2+1-n),wmxout(nlp2+1-n),omxout(nlp2+1-n)         
      write (6,202)  a(n,1),tmpint(n,1),prsint(n),wmix(n),o3mix(n)              
210   continue                                                                  
  203 format (1x,14x,f14.6,14x,2e14.6)                                          
  202 format(1x,2f14.6,e14.6,2e14.6)                                            
c   determine what formats to use***                                            
      nrep=nlp/5                                                                
      nrem=nlp-5*nrep                                                           
      if (nrem.eq.0) then                                                       
        nrep=nrep-1                                                             
        nrem=5                                                                  
      endif                                                                     
      nrep2=nl/5                                                                
      nrem2=nl-5*nrep2                                                          
      if (nrem2.eq.0) then                                                      
        nrep2=nrep2-1                                                           
        nrem2=5                                                                 
      endif                                                                     
      nrep3=nlp/4                                                               
      nrem3=nlp-4*nrep3                                                         
      if (nrem3.eq.0) then                                                      
        nrep3=nrep3-1                                                           
        nrem3=4                                                                 
      endif                                                                     
co222   *****************************************************                   
cc    rewind 66                                                                 
co222   *****************************************************                   
c***output temperatures                                                         
      do 800 iout=1,2                                                           
         if (iout.eq.1)  then                                                   
           write (16,701)                                                       
         else                                                                   
           write (16,702)                                                       
         endif                                                                  
701   format (6x,'data dtemp /')                                                
702   format (6x,'data stemp /')                                                
      nf=0                                                                      
co222   *****************************************************                   
c         save stemp                                                            
cc      if(iout.eq.2) write(66) (tmpint(nlp2-n,1),n=1,nlp)                      
      do 901 n=1,nlp                                                            
        sgtemp(n,1) = tmpint(nlp2-n,1)                                          
  901 continue                                                                  
co222   *****************************************************                   
      if (nrep.ne.0) then                                                       
         do 801 nr=1,nrep                                                       
           ns=5*(nr-1)+1                                                        
           nf=ns+4                                                              
           write (16,656) (tmpint(nlp2-n,1),n=ns,nf)                            
801      continue                                                               
      endif                                                                     
      nf=nf+1                                                                   
      if (nrem.eq.1) then                                                       
          write (16,616) (tmpint(nlp2-n,1),n=nf,nlp)                            
      endif                                                                     
      if (nrem.eq.2) then                                                       
          write (16,626) (tmpint(nlp2-n,1),n=nf,nlp)                            
      endif                                                                     
      if (nrem.eq.3) then                                                       
          write (16,636) (tmpint(nlp2-n,1),n=nf,nlp)                            
      endif                                                                     
      if (nrem.eq.4) then                                                       
          write (16,646) (tmpint(nlp2-n,1),n=nf,nlp)                            
      endif                                                                     
      if (nrem.eq.5) then                                                       
          write (16,606) (tmpint(nlp2-n,1),n=nf,nlp)                            
      endif                                                                     
616   format (5x,'*',1x,f12.6,'/')                                              
626   format (5x,'*',1x,f12.6,',',f12.6,'/')                                    
636   format (5x,'*',1x,f12.6,',',f12.6,',',f12.6,'/')                          
646   format (5x,'*',1x,f12.6,',',f12.6,',',f12.6,',',                          
     1                  f12.6,'/')                                              
656   format (5x,'*',1x,f12.6,',',f12.6,',',f12.6,',',                          
     1                  f12.6,',',f12.6,',')                                    
606   format (5x,'*',1x,f12.6,',',f12.6,',',f12.6,',',                          
     1                  f12.6,',',f12.6,'/')                                    
800   continue                                                                  
c***output gtemp                                                                
co222   *****************************************************                   
c         save gtemp                                                            
cc             write(66) (gtemp(n),n=1,nlp)                                     
      do 902 n=1,nlp                                                            
        sgtemp(n,2) = gtemp(n)                                                  
  902 continue                                                                  
co222   *****************************************************                   
      write (16,706)                                                            
706   format (6x,'data gtemp /')                                                
      nf=0                                                                      
      if (nrep3.ne.0) then                                                      
         do 805 nr=1,nrep3                                                      
           ns=4*(nr-1)+1                                                        
           nf=ns+3                                                              
           write (16,648) (gtemp(n),n=ns,nf)                                    
805      continue                                                               
      endif                                                                     
      nf=nf+1                                                                   
      if (nrem3.eq.1) then                                                      
          write (16,618) (gtemp(n),n=nf,nlp)                                    
      endif                                                                     
      if (nrem3.eq.2) then                                                      
          write (16,628) (gtemp(n),n=nf,nlp)                                    
      endif                                                                     
      if (nrem3.eq.3) then                                                      
          write (16,638) (gtemp(n),n=nf,nlp)                                    
      endif                                                                     
      if (nrem3.eq.4) then                                                      
          write (16,608) (gtemp(n),n=nf,nlp)                                    
      endif                                                                     
c***output wmix in 10**2 gm/gm                                                  
      do 630 n=2,nlp                                                            
        wmix(n)=wmix(n)*1.0e2                                                   
630   continue                                                                  
      nf=0                                                                      
      write (16,703)                                                            
703   format (6x,'data rr /')                                                   
      if (nrep2.ne.0) then                                                      
         do 802 nr=1,nrep2                                                      
           ns=5*(nr-1)+1                                                        
           nf=ns+4                                                              
           write (16,657) (wmix(nlp2-n),n=ns,nf)                                
802      continue                                                               
      endif                                                                     
      nf=nf+1                                                                   
      if (nrem2.eq.1) then                                                      
          write (16,617) (wmix(nlp2-n),n=nf,nl)                                 
      endif                                                                     
      if (nrem2.eq.2) then                                                      
          write (16,627) (wmix(nlp2-n),n=nf,nl)                                 
      endif                                                                     
      if (nrem2.eq.3) then                                                      
          write (16,637) (wmix(nlp2-n),n=nf,nl)                                 
      endif                                                                     
      if (nrem2.eq.4) then                                                      
          write (16,647) (wmix(nlp2-n),n=nf,nl)                                 
      endif                                                                     
      if (nrem2.eq.5) then                                                      
          write (16,607) (wmix(nlp2-n),n=nf,nl)                                 
      endif                                                                     
617   format (5x,'*',1x,e12.6,'/')                                              
627   format (5x,'*',1x,e12.6,',',e12.6,'/')                                    
637   format (5x,'*',1x,e12.6,',',e12.6,',',e12.6,'/')                          
647   format (5x,'*',1x,e12.6,',',e12.6,',',e12.6,',',                          
     1                  e12.6,'/')                                              
657   format (5x,'*',1x,e12.6,',',e12.6,',',e12.6,',',                          
     1                  e12.6,',',e12.6,',')                                    
607   format (5x,'*',1x,e12.6,',',e12.6,',',e12.6,',',                          
     1                  e12.6,',',e12.6,'/')                                    
c***output o3mix in 10**3 gm/gm                                                 
      do 631 n=2,nlp                                                            
         o3mix(n)=o3mix(n)*1.0e3                                                
631   continue                                                                  
      nf=0                                                                      
      write (16,704)                                                            
704   format (6x,'data qqo3 /')                                                 
      if (nrep2.ne.0) then                                                      
         do 803 nr=1,nrep2                                                      
           ns=5*(nr-1)+1                                                        
           nf=ns+4                                                              
           write (16,657) (o3mix(nlp2-n),n=ns,nf)                               
803      continue                                                               
      endif                                                                     
      nf=nf+1                                                                   
      if (nrem2.eq.1) then                                                      
          write (16,617) (o3mix(nlp2-n),n=nf,nl)                                
      endif                                                                     
      if (nrem2.eq.2) then                                                      
          write (16,627) (o3mix(nlp2-n),n=nf,nl)                                
      endif                                                                     
      if (nrem2.eq.3) then                                                      
          write (16,637) (o3mix(nlp2-n),n=nf,nl)                                
      endif                                                                     
      if (nrem2.eq.4) then                                                      
          write (16,647) (o3mix(nlp2-n),n=nf,nl)                                
      endif                                                                     
      if (nrem2.eq.5) then                                                      
          write (16,607) (o3mix(nlp2-n),n=nf,nl)                                
      endif                                                                     
c***output prout in cgs                                                         
      do 629 n=2,nlp2                                                           
        prout(n)=prout(n)*1.0e3                                                 
629   continue                                                                  
      write (16,705)                                                            
705   format (6x,'data ppress /')                                               
      nf=0                                                                      
      if (nrep3.ne.0) then                                                      
         do 804 nr=1,nrep3                                                      
           ns=4*(nr-1)+2                                                        
           nf=ns+3                                                              
           write (16,648) (prout(n),n=ns,nf)                                    
804      continue                                                               
      endif                                                                     
      nf=nf+1                                                                   
      if (nrem3.eq.1) then                                                      
          write (16,618) (prout(n),n=nf,nlp2)                                   
      endif                                                                     
      if (nrem3.eq.2) then                                                      
          write (16,628) (prout(n),n=nf,nlp2)                                   
      endif                                                                     
      if (nrem3.eq.3) then                                                      
          write (16,638) (prout(n),n=nf,nlp2)                                   
      endif                                                                     
      if (nrem3.eq.4) then                                                      
          write (16,608) (prout(n),n=nf,nlp2)                                   
      endif                                                                     
618   format (5x,'*',1x,e15.9,'/')                                              
628   format (5x,'*',1x,e15.9,',',e15.9,'/')                                    
638   format (5x,'*',1x,e15.9,',',e15.9,',',e15.9,'/')                          
648   format (5x,'*',1x,e15.9,',',e15.9,',',e15.9,',',e15.9,',')                
608   format (5x,'*',1x,e15.9,',',e15.9,',',e15.9,',',e15.9,'/')                
      return                                                                    
      end                                                                       
