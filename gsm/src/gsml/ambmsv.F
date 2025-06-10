#define CP cp
#define LEVS levs
#define RD rd
#define RERTH rerth
      subroutine ambmsv(kmx,si,sl,tov,am,bm,sv,gv,cm)                           
c$$$  subprogram documentation block                                            
c                .      .    .                                       .          
c subprogram:    ambmsv      computes constant matrices.                        
c   prgmmr: joseph sela      org: w/nmc23    date: 88-04-06                     
c                                                                               
c abstract: computes the 3 matrices am and bm and sv.                           
c   the matrices represent the linearized gravity wave terms                    
c   of the equations involved in the semi-implicit time integration.            
c   the matrices are dependent only on the vertical structure.                  
c   am is the divergence equations linear dependence on temperature.           
c   bm is the temperature equations linear dependence on divergence.           
c   sv is the continuity equations linear dependence on divergence.            
c                                                                               
c program history log:                                                          
c   88-04-06  joseph sela                                                       
c   93-02-23  mark iredell   compact vertical formulation                       
c                                                                               
c usage:    call ambmsv(km,sl,tov,am,bm,sv)                                     
c   input argument list:                                                        
c     km       - integer number of vertical levels.                             
c     sl       - real (km) sigma level values.                                  
c     tov      - real (km) reference temperatures.                              
c                                                                               
c   output argument list:                                                       
c     am       - real (km,km) such that dd(k)/dt = ... + am(k,j)*t(j)           
c     bm       - real (km,km) such that dt(k)/dt = ... + bm(k,j)*d(j)           
c     sv       - real (km) such that dq/dt = ... + sv(j)*d(j)                   
c                                                                               
c attributes:                                                                   
c   language: cray fortran.                                                     
c                                                                               
c$$$                                                                            
#include <paramodel.h>
      dimension si(kmx+1),sl(kmx),tov(kmx)                                      
      dimension am(kmx,kmx),bm(kmx,kmx),sv(kmx),gv(kmx),cm(kmx,kmx)             
      parameter(rd=rd_,cp=cp_,rerth=rerth_)                                     
      parameter(rocp=rd/cp,raa=rd/(rerth**2))                                   
c local array                                                                   
      dimension cd(levs_,levs_+1),ci(levs_+1,levs_),                            
     1          cq(levs_+1,levs_),cql(levs_,levs_)                              
      dimension rnu(levs_),rmu(levs_),ti(2:levs_),dt(levs_+1,levs_)             
      km=kmx                                                                    
      call bnmc(km,si,sl,cd,ci,cq,cql)                                          
      do 10 j=1,km                                                              
        sv(j)=cq(km+1,j)                                                        
        gv(j)=raa*tov(j)                                                        
        dt(1,j)=0.                                                              
        dt(km+1,j)=0.                                                           
10    continue                                                                  
      do 20 k=1,km-1                                                            
        ti(k+1)=0.                                                              
        del=si(k+1)-si(k)                                                       
        rsl=log(si(k+1)/si(k))/del                                              
        rnu(k)=(1.-rsl*si(k))/del                                               
        rmu(k)=(rsl*si(k+1)-1.)/del                                             
20    continue                                                                  
      rnu(km)=0.                                                                
      rmu(km)=1./si(km)                                                         
      do 40 k=2,km                                                              
        do 30 j=1,km                                                            
          ti(k)=ti(k)+ci(k,j)*tov(j)                                            
          dt(k,j)=(1-si(k))*cq(km+1,j)-cq(k,j)                                  
30      continue                                                                
40    continue                                                                  
      do 70 j=1,km                                                              
        do 50 k=1,km                                                            
          am(k,j)=-raa*cql(k,j)                                                 
          bm(k,j)=(rocp-1)*tov(k)*cq(km+1,j)                                    
     &           +rocp*tov(k)*(rnu(k)*dt(k+1,j)+rmu(k)*dt(k,j))                 
50      continue                                                                
        bm(j,j)=bm(j,j)-tov(j)                                                  
        do 60 k=1,km                                                            
        do 60 i=2,km                                                            
          bm(k,j)=bm(k,j)-cd(k,i)*ti(i)*dt(i,j)                                 
60      continue                                                                
70    continue                                                                  
#ifdef ORIGIN_THREAD
c$doacross share(km,gv,sv,am,bm,cm),local(j,k,i)                                
#endif
#ifdef CRAY_THREAD
cmic$ do all shared(km,gv,sv,am,bm,cm) private(j,k,i)                           
#endif
#ifdef OPENMP
!$omp parallel do private(j,k,i)
#endif
      do 100 j=1,km                                                             
        do 80 k=1,km                                                            
          cm(k,j)=gv(k)*sv(j)                                                   
80      continue                                                                
        do 90 k=1,km                                                            
        do 90 i=1,km                                                            
          cm(k,j)=cm(k,j)+am(k,i)*bm(i,j)                                       
90      continue                                                                
100   continue                                                                  
      return                                                                    
      end                                                                       
