#define IGRD12 igrd12
#define IGRD12P igrd12p
#define IGRD1P igrd1p
#define IGRD igrd
#define JGRD12P jgrd12p
#define JGRD1P jgrd1p
#define LEVS levs
#define LONFP lonfp
      subroutine table                                                          
#ifdef CRAY_THREAD
cfpp$ noconcur r                                                                
#endif
c     subroutine table computes table entries used in the longwave radia        
c     program. also calculated are indices used in strip-mining and for         
c     some pre-computable functions.                                            
c         inputs:                                                               
c         outputs:                                                              
c       em1,em1wde,table1,table2,table3         tabcom                          
c       em3,source,dsrce,ind,indx2,kmaxv        tabcom                          
c       kmaxvm,                                 tabcom                          
c       ao3rnd,bo3rnd,ab15                      bandta                          
c       ab15wd,skc1r,sko3r,sko2d                bdwide                          
c                                                                               
#include <paramodel.h>
#include <hcon.h>
#include <rdparm.h>
#include <rnddta.h>
#include <tabcom.h>
c                                                                               
      dimension sum(28,180),pertsm(28,180),sum3(28,180),                        
     1 sumwde(28,180),srcwd(28,nblx),src1nb(28,nblw),dbdtnb(28,nblw)            
      dimension zmass(181),zroot(181),sc(28),dsc(28),xtemv(28),                 
     1 tfour(28),fortcu(28),x(28),x1(28),x2(180),srcs(28),                      
     2 sum4(28),sum6(28),sum7(28),sum8(28),sum4wd(28),                          
     3 r1(28),r2(28),s2(28),t3(28),r1wd(28)                                     
      dimension expo(180),fac(180)                                              
      dimension cnusb(30),dnusb(30)                                             
      dimension alfanb(nblw),arotnb(nblw)                                       
      dimension anb(nblw),bnb(nblw),centnb(nblw),delnb(nblw),                   
     1          betanb(nblw)                                                    
      common/tbltmp/ delcm(nbly)                                                
c****************************************                                       
c***compute local quantities and ao3,bo3,ab15                                   
c....for narrow-bands...                                                        
      do 101 n=1,nblw                                                           
      anb(n)=arndm(n)                                                           
      bnb(n)=brndm(n)                                                           
      centnb(n)=haf*(bandlo(n)+bandhi(n))                                       
      delnb(n)=bandhi(n)-bandlo(n)                                              
      betanb(n)=betad(n)                                                        
101   continue                                                                  
      ab15(1)=anb(57)*bnb(57)                                                   
      ab15(2)=anb(58)*bnb(58)                                                   
c....for wide bands...                                                          
      ab15wd=awide*bwide                                                        
c                                                                               
c***compute indices: ind,indx2,kmaxv                                            
      do 111 i=1,imax                                                           
      ind(i)=i                                                                  
111   continue                                                                  
      icnt=0                                                                    
      do 113 i1=1,l                                                             
        i2e=lp1-i1                                                              
        do 115 i2=1,i2e                                                         
          icnt=icnt+1                                                           
          indx2(icnt)=lp1*(i2-1)+lp2*i1                                         
115     continue                                                                
113   continue                                                                  
      kmaxv(1)=1                                                                
      do 117 i=2,l                                                              
      kmaxv(i)=kmaxv(i-1)+(lp2-i)                                               
117   continue                                                                  
      kmaxvm=kmaxv(l)                                                           
c***compute ratios of cont. coeffs                                              
      skc1r=betawd/betinw                                                       
      sko3r=betad(61)/betinw                                                    
      sko2d=one/betinw                                                          
c                                                                               
c****begin table computations here***                                           
c***compute temps, masses for table entries                                     
c---note: the dimensioning and initialization of xtemv and other arrays         
c   with dimension of 28 imply a restriction of model temperatures from         
c   100k to 370k.                                                               
c---the dimensioning of zmass,zroot and other arrays with dimension of          
c   180 imply a restriction of model h2o amounts such that optical paths        
c   are between 10**-16 and 10**2, in cgs units.                                
      zmass(1)=h1m16                                                            
      do 201 j=1,180                                                            
      jp=j+1                                                                    
      zroot(j)=sqrt(zmass(j))                                                   
      zmass(jp)=zmass(j)*h1p25892                                               
201   continue                                                                  
      do 203 i=1,28                                                             
      xtemv(i)=hninety+ten*i                                                    
      tfour(i)=xtemv(i)*xtemv(i)*xtemv(i)*xtemv(i)                              
      fortcu(i)=four*xtemv(i)*xtemv(i)*xtemv(i)                                 
203   continue                                                                  
c******the computation of source,dsrce is  needed only                          
c   for the combined wide-band case.to obtain them,the source                   
c   must be computed for each of the (nblx) wide bands(=srcwd)                  
c   then combined (using iband) into source.                                    
      do 205 n=1,nbly                                                           
      do 205 i=1,28                                                             
      source(i,n)=zero                                                          
205   continue                                                                  
      do 207 n=1,nblx                                                           
      do 207 i=1,28                                                             
      srcwd(i,n)=zero                                                           
207   continue                                                                  
c---begin freq. loop (on n)                                                     
      do 211 n=1,nblx                                                           
        if (n.le.46) then                                                       
c***the 160-1200 band cases                                                     
          cent=centnb(n+16)                                                     
          del=delnb(n+16)                                                       
          bdlo=bandlo(n+16)                                                     
          bdhi=bandhi(n+16)                                                     
        endif                                                                   
        if (n.eq.nblx) then                                                     
c***the 2270-2380 band case                                                     
          cent=centnb(nblw)                                                     
          del=delnb(nblw)                                                       
          bdlo=bandlo(nblw)                                                     
          bdhi=bandhi(nblw)                                                     
        endif                                                                   
c***for purposes of accuracy, all evaluations of planck fctns are made          
c  on 10 cm-1 intervals, then summed into the (nblx) wide bands.                
      nsubds=(del-h1m3)/10+1                                                    
      do 213 nsb=1,nsubds                                                       
      if (nsb.ne.nsubds) then                                                   
        cnusb(nsb)=ten*(nsb-1)+bdlo+five                                        
        dnusb(nsb)=ten                                                          
      else                                                                      
        cnusb(nsb)=haf*(ten*(nsb-1)+bdlo+bdhi)                                  
        dnusb(nsb)=bdhi-(ten*(nsb-1)+bdlo)                                      
      endif                                                                     
      c1=(h37412m5)*cnusb(nsb)**3                                               
c---begin temp. loop (on i)                                                     
      do 215 i=1,28                                                             
      x(i)=h1p4387*cnusb(nsb)/xtemv(i)                                          
      x1(i)=exp(x(i))                                                           
      srcs(i)=c1/(x1(i)-one)                                                    
      srcwd(i,n)=srcwd(i,n)+srcs(i)*dnusb(nsb)                                  
215   continue                                                                  
213   continue                                                                  
211   continue                                                                  
c***the following loops create the combined wide band quantities source         
c   and dsrce                                                                   
      do 221 n=1,40                                                             
      do 221 i=1,28                                                             
      source(i,iband(n))=source(i,iband(n))+srcwd(i,n)                          
221   continue                                                                  
      do 223 n=9,nbly                                                           
      do 223 i=1,28                                                             
      source(i,n)=srcwd(i,n+32)                                                 
223   continue                                                                  
      do 225 n=1,nbly                                                           
      do 225 i=1,27                                                             
      dsrce(i,n)=(source(i+1,n)-source(i,n))*hp1                                
225   continue                                                                  
      do 231 n=1,nblw                                                           
      alfanb(n)=bnb(n)*anb(n)                                                   
      arotnb(n)=sqrt(alfanb(n))                                                 
231   continue                                                                  
c***first compute planck fctns (src1nb) and derivatives (dbdtnb) for            
c   use in table evaluations. these are different from source,dsrce             
c   because different frequency pts are used in evaluation, the freq.           
c   ranges are different, and the derivative algorithm is different.            
c                                                                               
      do 301 n=1,nblw                                                           
      cent=centnb(n)                                                            
      del=delnb(n)                                                              
c---note: at present, the ia loop is only used for ia=2. the loop struct        
c   is kept so that in the future, we may use a quadrature scheme for           
c   the planck fctn evaluation, rather than use the mid-band frequency.         
      do 303 ia=1,3                                                             
      anu=cent+haf*(ia-2)*del                                                   
      c1=(h37412m5)*anu*anu*anu+h1m20                                           
c---temperature loop---      
#ifdef SX6
!CDIR NOVECTOR      
#endif                                             
      do 305 i=1,28                                                             
         x(i)=h1p4387*anu/xtemv(i)                                              
         x1(i)=exp(x(i))                                                        
         sc(i)=c1/((x1(i)-one)+h1m20)                                           
         dsc(i)=sc(i)*sc(i)*x(i)*x1(i)/(xtemv(i)*c1)                            
305      continue                                                               
      if (ia.eq.2) then                                                         
         do 307 i=1,28                                                          
         src1nb(i,n)=del*sc(i)                                                  
         dbdtnb(i,n)=del*dsc(i)                                                 
307      continue                                                               
      endif                                                                     
303   continue                                                                  
301   continue                                                                  
c***next compute r1,r2,s2,and t3- coefficients used for e3 function             
c   when the optical path is less than 10-4. in this case, we assume a          
c   different dependence on (zmass).                                            
c---also obtain r1wd, which is r1 summed over the 160-560 cm-1 range            
      do 311 i=1,28                                                             
      sum4(i)=zero                                                              
      sum6(i)=zero                                                              
      sum7(i)=zero                                                              
      sum8(i)=zero                                                              
      sum4wd(i)=zero                                                            
311   continue                                                                  
      do 313 n=1,nblw                                                           
      cent=centnb(n)                                                            
c***perform summations for freq. ranges of 0-560,1200-2200 cm-1 for sum4        
c   sum6,sum7,sum8                                                              
      if (cent.lt.560. .or. cent.gt.1200..and.cent.le.2200.) then               
         do 315 i=1,28                                                          
         sum4(i)=sum4(i)+src1nb(i,n)                                            
         sum6(i)=sum6(i)+dbdtnb(i,n)                                            
         sum7(i)=sum7(i)+dbdtnb(i,n)*arotnb(n)                                  
         sum8(i)=sum8(i)+dbdtnb(i,n)*alfanb(n)                                  
315      continue                                                               
      endif                                                                     
c***perform summations over 160-560 cm-1 freq range for e1 calcs (sum4wd        
      if (cent.gt.160. .and. cent.lt.560.) then                                 
         do 316 i=1,28                                                          
         sum4wd(i)=sum4wd(i)+src1nb(i,n)                                        
316      continue                                                               
      endif                                                                     
313   continue                                                                  
      do 317 i=1,28                                                             
      r1(i)=sum4(i)/tfour(i)                                                    
      r2(i)=sum6(i)/fortcu(i)                                                   
      s2(i)=sum7(i)/fortcu(i)                                                   
      t3(i)=sum8(i)/fortcu(i)                                                   
      r1wd(i)=sum4wd(i)/tfour(i)                                                
317   continue                                                                  
      do 401 j=1,180                                                            
      do 401 i=1,28                                                             
      sum(i,j)=zero                                                             
      pertsm(i,j)=zero                                                          
      sum3(i,j)=zero                                                            
      sumwde(i,j)=zero                                                          
401   continue                                                                  
c---frequency loop begins---                                                    
      do 411 n=1,nblw                                                           
      cent=centnb(n)                                                            
c***perform calculations for freq. ranges of 0-560,1200-2200 cm-1               
      if (cent.lt.560. .or. cent.gt.1200..and.cent.le.2200.) then               
         do 413 j=1,180                                                         
         x2(j)=arotnb(n)*zroot(j)                                               
         expo(j)=exp(-x2(j))                                                    
413      continue                                                               
         do 415 j=1,180                                                         
         if (x2(j).ge.hundred) then                                             
              expo(j)=zero                                                      
         endif                                                                  
415      continue                                                               
         do 417 j=121,180                                                       
         fac(j)=zmass(j)*(one-(one+x2(j))*expo(j))/(x2(j)*x2(j))                
417      continue                                                               
         do 419 j=1,180                                                         
         do 419 i=1,28                                                          
         sum(i,j)=sum(i,j)+src1nb(i,n)*expo(j)                                  
         pertsm(i,j)=pertsm(i,j)+dbdtnb(i,n)*expo(j)                            
419      continue                                                               
         do 421 j=121,180                                                       
         do 421 i=1,28                                                          
         sum3(i,j)=sum3(i,j)+dbdtnb(i,n)*fac(j)                                 
421      continue                                                               
      endif                                                                     
c---compute sum over 160-560 cm-1 range for use in e1 calcs (sumwde)            
      if (cent.gt.160. .and. cent.lt.560.) then                                 
         do 420 j=1,180                                                         
         do 420 i=1,28                                                          
         sumwde(i,j)=sumwde(i,j)+src1nb(i,n)*expo(j)                            
420      continue                                                               
      endif                                                                     
411   continue                                                                  
      do 431 j=1,180                                                            
      do 431 i=1,28                                                             
      em1(i,j)=sum(i,j)/tfour(i)                                                
      table1(i,j)=pertsm(i,j)/fortcu(i)                                         
431   continue                                                                  
      do 433 j=121,180                                                          
      do 433 i=1,28                                                             
      em3(i,j)=sum3(i,j)/fortcu(i)                                              
433   continue                                                                  
      do 441 j=1,179                                                            
      do 441 i=1,28                                                             
      table2(i,j)=(table1(i,j+1)-table1(i,j))*ten                               
441   continue                                                                  
      do 443 j=1,180                                                            
      do 443 i=1,27                                                             
      table3(i,j)=(table1(i+1,j)-table1(i,j))*hp1                               
443   continue                                                                  
      do 445 i=1,28                                                             
      table2(i,180)=zero                                                        
445   continue                                                                  
      do 447 j=1,180                                                            
      table3(28,j)=zero                                                         
447   continue                                                                  
      do 449 j=1,2                                                              
      do 449 i=1,28                                                             
      em1(i,j)=r1(i)                                                            
449   continue                                                                  
      do 451 j=1,120                                                            
      do 451 i=1,28                                                             
      em3(i,j)=r2(i)/two-s2(i)*sqrt(zmass(j))/three+t3(i)*zmass(j)/eight        
451   continue                                                                  
      do 453 j=121,180                                                          
      do 453 i=1,28                                                             
      em3(i,j)=em3(i,j)/zmass(j)                                                
453   continue                                                                  
c***now compute e1 tables for 160-560 cm-1 bands only.                          
c   we use r1wd and sumwde obtained above.                                      
      do 501 j=1,180                                                            
      do 501 i=1,28                                                             
      em1wde(i,j)=sumwde(i,j)/tfour(i)                                          
501   continue                                                                  
      do 503 j=1,2                                                              
      do 503 i=1,28                                                             
      em1wde(i,j)=r1wd(i)                                                       
503   continue                                                                  
      return                                                                    
      end                                                                       
